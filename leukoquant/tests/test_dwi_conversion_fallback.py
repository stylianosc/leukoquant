"""Tests for DWI DICOM conversion fallback cascade and input validation."""

import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import nibabel as nib
import numpy as np
import pytest

from leukoquant.utils.dwi_utils import (
    _convert_dicom_with_fallback,
    _find_converter_candidates,
    _record_input_failure,
    InputDataError,
    prepare_dwi_input,
)


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


def _create_mock_nifti_4d(filepath: Path, shape=(10, 10, 10, 7)):
    data = np.ones(shape, dtype=np.float32)
    img = nib.Nifti1Image(data, affine=np.eye(4))
    nib.save(img, str(filepath))


def _create_mock_nifti_3d(filepath: Path, shape=(10, 10, 10)):
    data = np.ones(shape, dtype=np.float32)
    img = nib.Nifti1Image(data, affine=np.eye(4))
    nib.save(img, str(filepath))


class TestFindConverterCandidates:
    def test_finds_available_converters_in_order(self, temp_dir):
        # Create dummy binaries
        dcm2niix = temp_dir / "dcm2niix"
        dcm2niix.write_text("#!/bin/sh\nexit 0\n")
        dcm2niix.chmod(0o755)

        dcm2niix_nonspatial = temp_dir / "dcm2niix_nonspatial"
        dcm2niix_nonspatial.write_text("#!/bin/sh\nexit 0\n")
        dcm2niix_nonspatial.chmod(0o755)

        mri_convert = temp_dir / "mri_convert"
        mri_convert.write_text("#!/bin/sh\nexit 0\n")
        mri_convert.chmod(0o755)

        with patch("shutil.which") as mock_which:
            mock_which.side_effect = lambda cmd: {
                "dcm2niix": str(dcm2niix),
                "dcm2niix_nonspatial": str(dcm2niix_nonspatial),
                "mri_convert": str(mri_convert),
            }.get(cmd)

            candidates = _find_converter_candidates(str(dcm2niix))
            names = [c[0] for c in candidates]
            assert names == ["dcm2niix", "dcm2niix_nonspatial", "mri_convert"]


class TestConvertDicomWithFallback:
    def test_primary_succeeds_first_time(self, temp_dir):
        src_dir = temp_dir / "dicom_src"
        src_dir.mkdir()
        out_dir = temp_dir / "out"
        out_dir.mkdir()

        def fake_run(cmd, capture_output=True, text=True):
            # Simulate successful dcm2niix
            dwi_path = out_dir / "dwi.nii.gz"
            _create_mock_nifti_4d(dwi_path, (10, 10, 10, 7))
            np.savetxt(out_dir / "dwi.bval", [0, 1000, 1000, 1000, 1000, 1000, 1000])
            np.savetxt(out_dir / "dwi.bvec", np.zeros((3, 7)))
            proc = MagicMock()
            proc.returncode = 0
            return proc

        with patch("leukoquant.utils.dwi_utils._find_converter_candidates") as mock_cands, \
             patch("subprocess.run", side_effect=fake_run) as mock_sub:
            mock_cands.return_value = [("dcm2niix", "/fake/dcm2niix")]
            dwi, bvec, bval = _convert_dicom_with_fallback(src_dir, out_dir)

            assert dwi.name == "dwi.nii.gz"
            assert bvec is not None
            assert bval is not None
            assert mock_sub.call_count == 1

    def test_fallback_to_nonspatial_on_primary_failure(self, temp_dir):
        src_dir = temp_dir / "dicom_src"
        src_dir.mkdir()
        out_dir = temp_dir / "out"
        out_dir.mkdir()

        call_log = []

        def fake_run(cmd, capture_output=True, text=True):
            call_log.append(cmd[0])
            proc = MagicMock()
            if "dcm2niix_nonspatial" in cmd[0]:
                dwi_path = out_dir / "dwi.nii.gz"
                _create_mock_nifti_4d(dwi_path, (10, 10, 10, 7))
                np.savetxt(out_dir / "dwi.bval", [0, 1000, 1000, 1000, 1000, 1000, 1000])
                np.savetxt(out_dir / "dwi.bvec", np.zeros((3, 7)))
                proc.returncode = 0
            else:
                proc.returncode = 1
                proc.stderr = "InstanceNumberOrderIsNotSpatial error"
            return proc

        with patch("leukoquant.utils.dwi_utils._find_converter_candidates") as mock_cands, \
             patch("subprocess.run", side_effect=fake_run):
            mock_cands.return_value = [
                ("dcm2niix", "/fake/dcm2niix"),
                ("dcm2niix_nonspatial", "/fake/dcm2niix_nonspatial"),
            ]
            dwi, bvec, bval = _convert_dicom_with_fallback(src_dir, out_dir)

            assert dwi.name == "dwi.nii.gz"
            assert call_log == ["/fake/dcm2niix", "/fake/dcm2niix_nonspatial"]

    def test_rejects_flattened_3d_when_bvals_are_multivolume(self, temp_dir):
        src_dir = temp_dir / "dicom_src"
        src_dir.mkdir()
        out_dir = temp_dir / "out"
        out_dir.mkdir()

        def fake_run(cmd, capture_output=True, text=True):
            dwi_path = out_dir / "dwi.nii.gz"
            # Output is incorrectly 3D
            _create_mock_nifti_3d(dwi_path, (10, 10, 10))
            np.savetxt(out_dir / "dwi.bval", [0, 1000, 1000, 1000, 1000, 1000, 1000])
            np.savetxt(out_dir / "dwi.bvec", np.zeros((3, 7)))
            proc = MagicMock()
            proc.returncode = 0
            return proc

        with patch("leukoquant.utils.dwi_utils._find_converter_candidates") as mock_cands, \
             patch("subprocess.run", side_effect=fake_run):
            mock_cands.return_value = [("dcm2niix", "/fake/dcm2niix")]
            with pytest.raises(InputDataError, match="Corrupted conversion: image is 3D"):
                _convert_dicom_with_fallback(src_dir, out_dir)


class TestInputChecksAndAtomicFailureMarker:
    """The marker is written to the session folder the caller names, and only
    for unusable input (InputDataError), never for other failures."""

    @staticmethod
    def _dwi(temp_dir, n_vols, bvals):
        dwi_file = temp_dir / "raw_dwi.nii.gz"
        _create_mock_nifti_4d(dwi_file, (10, 10, 10, n_vols))
        bval_file = temp_dir / "raw_dwi.bval"
        bvec_file = temp_dir / "raw_dwi.bvec"
        np.savetxt(bval_file, bvals)
        np.savetxt(bvec_file, np.ones((3, len(bvals))))
        return dwi_file, bvec_file, bval_file

    def test_record_input_failure_writes_marker_in_named_session_dir(self, temp_dir):
        session_dir = temp_dir / "subj-001"
        _record_input_failure(session_dir, "Corrupted DWI missing slices", source="/raw/dwi.nii.gz")

        marker = session_dir / "input_checks" / "dwi_failed.txt"
        assert marker.read_text().splitlines() == ["Corrupted DWI missing slices", "input: /raw/dwi.nii.gz"]
        assert not list((session_dir / "input_checks").glob(".tmp_*"))

    def test_too_few_directions_records_marker(self, temp_dir):
        session_dir = temp_dir / "subj-002"
        out_dir = temp_dir / "scratch" / "scan_0"
        out_dir.mkdir(parents=True)
        dwi, bvec, bval = self._dwi(temp_dir, 4, [0, 1000, 1000, 1000])

        with pytest.raises(InputDataError, match="minimum 6 diffusion-weighted"):
            prepare_dwi_input(str(dwi), bvecs=str(bvec), bvals=str(bval),
                              output_dir=str(out_dir), session_dir=str(session_dir))

        assert "minimum 6 diffusion-weighted" in (session_dir / "input_checks" / "dwi_failed.txt").read_text()

    def test_b0_volumes_do_not_count_as_directions(self, temp_dir):
        session_dir = temp_dir / "subj-004"
        out_dir = temp_dir / "scratch" / "scan_0"
        out_dir.mkdir(parents=True)
        # 8 volumes but only 3 are diffusion-weighted
        dwi, bvec, bval = self._dwi(temp_dir, 8, [0, 0, 0, 0, 0, 1000, 1000, 1000])

        with pytest.raises(InputDataError):
            prepare_dwi_input(str(dwi), bvecs=str(bvec), bvals=str(bval),
                              output_dir=str(out_dir), session_dir=str(session_dir))
        assert (session_dir / "input_checks" / "dwi_failed.txt").exists()

    def test_bval_volume_mismatch_records_marker(self, temp_dir):
        session_dir = temp_dir / "subj-003"
        out_dir = temp_dir / "scratch" / "scan_0"
        out_dir.mkdir(parents=True)
        dwi, bvec, bval = self._dwi(temp_dir, 7, [0] + [1000] * 7)  # 7 volumes, 8 b-values

        with pytest.raises(InputDataError, match="Volume count mismatch"):
            prepare_dwi_input(str(dwi), bvecs=str(bvec), bvals=str(bval),
                              output_dir=str(out_dir), session_dir=str(session_dir))
        assert "Volume count mismatch" in (session_dir / "input_checks" / "dwi_failed.txt").read_text()

    def test_other_errors_leave_no_marker(self, temp_dir):
        """A transient/environment failure must never get a healthy session skipped."""
        session_dir = temp_dir / "subj-005"
        out_dir = temp_dir / "scratch" / "scan_0"
        out_dir.mkdir(parents=True)
        dwi, bvec, bval = self._dwi(temp_dir, 7, [0] + [1000] * 6)

        with patch("leukoquant.utils.dwi_utils._prepare_dwi_input_impl", side_effect=OSError("disk full")):
            with pytest.raises(OSError):
                prepare_dwi_input(str(dwi), bvecs=str(bvec), bvals=str(bval),
                                  output_dir=str(out_dir), session_dir=str(session_dir))
        assert not (session_dir / "input_checks").exists()

    def test_no_marker_without_session_dir(self, temp_dir):
        out_dir = temp_dir / "scratch" / "scan_0"
        out_dir.mkdir(parents=True)
        dwi, bvec, bval = self._dwi(temp_dir, 4, [0, 1000, 1000, 1000])

        with pytest.raises(InputDataError):
            prepare_dwi_input(str(dwi), bvecs=str(bvec), bvals=str(bval), output_dir=str(out_dir))
        assert not list(temp_dir.rglob("input_checks"))

    def test_converter_without_gradient_tables_is_rejected(self, temp_dir):
        """mri_convert writes the image only; that must not count as a conversion."""
        src_dir = temp_dir / "dicom_src"
        src_dir.mkdir()
        out_dir = temp_dir / "out"
        out_dir.mkdir()

        def fake_run(cmd, capture_output=True, text=True, env=None):
            _create_mock_nifti_4d(out_dir / "data.nii.gz", (10, 10, 10, 7))
            proc = MagicMock()
            proc.returncode = 0
            return proc

        with patch("leukoquant.utils.dwi_utils._find_converter_candidates",
                   return_value=[("mri_convert", "/fake/mri_convert")]), \
             patch("subprocess.run", side_effect=fake_run):
            (src_dir / "a.dcm").write_bytes(b"x")
            with pytest.raises(InputDataError, match="no bval/bvec"):
                _convert_dicom_with_fallback(src_dir, out_dir)


pydicom = pytest.importorskip("pydicom")


def _write_dicom(path: Path, slice_location=None, ipp_z=None):
    from pydicom.dataset import FileDataset, FileMetaDataset
    from pydicom.uid import ExplicitVRLittleEndian, generate_uid

    meta = FileMetaDataset()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian
    meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.4"
    meta.MediaStorageSOPInstanceUID = generate_uid()
    ds = FileDataset(str(path), {}, file_meta=meta, preamble=b"\0" * 128)
    if slice_location is not None:
        ds.SliceLocation = slice_location
    if ipp_z is not None:
        ds.ImagePositionPatient = [0.0, 0.0, ipp_z]
    ds.save_as(str(path), enforce_file_format=True)


class TestSliceFilter:
    def test_slice_location_prefers_slicelocation_then_ipp(self, temp_dir):
        from leukoquant.utils.dwi_utils import _extract_dicom_slice_location

        _write_dicom(temp_dir / "a.dcm", slice_location=12.5, ipp_z=99.0)
        _write_dicom(temp_dir / "b.dcm", ipp_z=-7.25)
        _write_dicom(temp_dir / "c.dcm")
        (temp_dir / "d.dcm").write_bytes(b"not a dicom")

        assert _extract_dicom_slice_location(temp_dir / "a.dcm") == 12.5
        assert _extract_dicom_slice_location(temp_dir / "b.dcm") == -7.25
        assert _extract_dicom_slice_location(temp_dir / "c.dcm") is None
        assert _extract_dicom_slice_location(temp_dir / "d.dcm") is None

    def test_trailing_partial_volume_is_trimmed(self, temp_dir):
        from leukoquant.utils.dwi_utils import _filter_incomplete_slices

        src = temp_dir / "dicoms"
        src.mkdir()
        n = 0
        for vol in range(6):            # 6 complete volumes of 5 slices
            for z in range(5):
                _write_dicom(src / f"v{vol}_z{z}.dcm", slice_location=float(z))
                n += 1
        for z in (5, 6):                # aborted acquisition: 2 stray slices
            _write_dicom(src / f"partial_z{z}.dcm", slice_location=float(z))

        filtered = _filter_incomplete_slices(src)
        assert filtered != src
        assert len(list(filtered.iterdir())) == n

    def test_complete_series_is_left_alone(self, temp_dir):
        from leukoquant.utils.dwi_utils import _filter_incomplete_slices

        src = temp_dir / "dicoms"
        src.mkdir()
        for vol in range(6):
            for z in range(5):
                _write_dicom(src / f"v{vol}_z{z}.dcm", slice_location=float(z))
        assert _filter_incomplete_slices(src) == src
