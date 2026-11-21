"""Smoke tests for metrics_calc module.

These tests verify basic functionality without requiring full NIfTI fixtures.
"""

import argparse
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import nibabel as nib
import pytest
from leukoquant.utils.metrics_calc import (
    _load,
    _load_inputs,
    _voxel_vol_mm3,
    _assert_shape_compatibility,
    _load_qc_report,
    _qc_pass_for_tract,
)


@pytest.fixture
def temp_dir():
    """Create temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def simple_nifti(temp_dir):
    """Create a simple test NIfTI file."""
    data = np.random.rand(10, 10, 10)
    affine = np.diag([2.0, 2.0, 2.0, 1.0])  # 2mm isotropic
    img = nib.Nifti1Image(data, affine=affine)
    path = temp_dir / "test.nii.gz"
    nib.save(img, path)
    return path


class TestLoadFunction:
    """Test _load function."""

    def test_load_valid_nifti(self, simple_nifti):
        """Test loading a valid NIfTI file."""
        img = _load(str(simple_nifti))
        assert img is not None
        assert img.shape == (10, 10, 10)

    def test_load_nonexistent_raises(self):
        """Test that loading nonexistent file raises error."""
        with pytest.raises(Exception):  # FileNotFoundError or nibabel error
            _load("/nonexistent/file.nii.gz")


class TestVoxelVolFunction:
    """Test _voxel_vol_mm3 function."""

    def test_isotropic_voxel(self, simple_nifti):
        """Test voxel volume calculation for isotropic voxels."""
        img = _load(str(simple_nifti))
        voxel_vol = _voxel_vol_mm3(img)
        # 2mm x 2mm x 2mm = 8 mm³
        assert abs(voxel_vol - 8.0) < 0.01

    def test_anisotropic_voxel(self, temp_dir):
        """Test voxel volume calculation for anisotropic voxels."""
        data = np.zeros((5, 5, 5))
        affine = np.diag([1.0, 2.0, 4.0, 1.0])  # 1mm x 2mm x 4mm
        img = nib.Nifti1Image(data, affine=affine)
        path = temp_dir / "aniso.nii.gz"
        nib.save(img, path)

        loaded = _load(str(path))
        voxel_vol = _voxel_vol_mm3(loaded)
        # 1 x 2 x 4 = 8 mm³
        assert abs(voxel_vol - 8.0) < 0.01


class TestAssertShapeCompatibility:
    """Test _assert_shape_compatibility function."""

    def test_compatible_shapes(self):
        """Test arrays with matching shapes."""
        arrays = {
            "arr1": np.zeros((5, 5, 5)),
            "arr2": np.ones((5, 5, 5)),
        }
        # Should not raise
        _assert_shape_compatibility((5, 5, 5), arrays, "test_group")

    def test_incompatible_shapes_raises(self):
        """Test arrays with mismatched shapes raise ValueError."""
        arrays = {
            "arr1": np.zeros((5, 5, 5)),
            "arr2": np.ones((6, 6, 6)),
        }
        with pytest.raises(ValueError, match="Shape mismatch detected"):
            _assert_shape_compatibility((5, 5, 5), arrays, "test_group")

    def test_empty_arrays_dict(self):
        """Test with empty arrays dictionary."""
        arrays = {}
        # Should not raise
        _assert_shape_compatibility((5, 5, 5), arrays, "test_group")

    def test_single_mismatched_raises(self):
        """Test that single mismatched array raises."""
        arrays = {
            "correct": np.zeros((5, 5, 5)),
            "wrong": np.ones((4, 4, 4)),
        }
        with pytest.raises(ValueError, match="Shape mismatch detected"):
            _assert_shape_compatibility((5, 5, 5), arrays, "test_group")


class TestLoadQcReport:
    """Test _load_qc_report and _qc_pass_for_tract functions."""

    def test_load_qc_report_maps_dir_and_base_name(self, temp_dir):
        qc_csv = temp_dir / "qc_report.csv"
        df = pd.DataFrame([
            {"file": "/path/to/dpath/lh.cst_avg16_syn_bbr/path.pd.trk", "qc_pass": False, "tract_name": "lh.cst_avg16_syn_bbr"},
            {"file": "/path/to/dpath/acomm_avg16_syn_bbr/path.pd.trk", "qc_pass": True, "tract_name": "acomm_avg16_syn_bbr"},
        ])
        df.to_csv(qc_csv, index=False)

        qc_map = _load_qc_report(str(qc_csv))
        assert qc_map["lh.cst_avg16_syn_bbr"] is False
        assert qc_map["lh.cst"] is False
        assert qc_map["acomm_avg16_syn_bbr"] is True
        assert qc_map["acomm"] is True

    def test_qc_pass_for_tract_matches(self):
        qc_map = {"lh.cst": False, "acomm": True}
        assert _qc_pass_for_tract("lh.cst_avg16_syn_bbr", "/dummy/tracts.nii.gz", qc_map) is False
        assert _qc_pass_for_tract("acomm_avg16_syn_bbr", "/dummy/tracts.nii.gz", qc_map) is True


class TestTractLoadingFallback:
    """Test tract loading with atlas fallback for missing tracts."""

    def test_missing_tract_in_tractography_atlas_mode(self, temp_dir):
        import argparse
        from leukoquant.utils.metrics_calc import _load_inputs

        affine = np.eye(4)
        t1_img = nib.Nifti1Image(np.ones((10, 10, 10)), affine=affine)
        t1_path = temp_dir / "t1.nii.gz"
        nib.save(t1_img, t1_path)

        # 1 tract in tractography: acomm
        tract_dir = temp_dir / "tractography" / "acomm_avg16_syn_bbr"
        tract_dir.mkdir(parents=True)
        acomm_img = nib.Nifti1Image(np.full((10, 10, 10), 0.8), affine=affine)
        nib.save(acomm_img, tract_dir / "tracts.nii.gz")

        # 2 tracts in atlas: acomm and lh.cst (lh.cst omitted by tractography due to FOV)
        atlas_dir = temp_dir / "tract_atlas"
        atlas_dir.mkdir(parents=True)
        nib.save(acomm_img, atlas_dir / "acomm.nii.gz")
        cst_img = nib.Nifti1Image(np.full((10, 10, 10), 0.5), affine=affine)
        nib.save(cst_img, atlas_dir / "lh.cst.nii.gz")

        qc_csv = temp_dir / "qc_report.csv"
        pd.DataFrame([
            {"file": str(tract_dir / "path.pd.trk"), "qc_pass": True, "tract_name": "acomm_avg16_syn_bbr"},
            {"file": "missing/lh.cst_avg16_syn_bbr/path.pd.trk", "qc_pass": False, "tract_name": "lh.cst_avg16_syn_bbr"},
        ]).to_csv(qc_csv, index=False)

        args = argparse.Namespace(
            subject="test_subj",
            t1_path=str(t1_path),
            lesion_path="",
            tracts_path=str(temp_dir / "tractography" / "*" / "tracts.nii.gz"),
            tract_atlas_path=str(atlas_dir / "*.nii.gz"),
            qc_report=str(qc_csv),
            tract_mode="tractography-atlas",
            skip_skullstrip=True,
            t1_brain_path="",
            skeleton_path="",
            maps="",
        )

        bundle = _load_inputs(args)

        assert "acomm" in bundle.tract_arrays
        assert "lh.cst" in bundle.tract_arrays
        assert bundle.tract_source["acomm"] == "tractography"
        assert bundle.tract_source["lh.cst"] == "atlas"

    def test_atlas_mode_loads_all_atlas_tracts(self, temp_dir):
        import argparse
        from leukoquant.utils.metrics_calc import _load_inputs

        affine = np.eye(4)
        t1_img = nib.Nifti1Image(np.ones((10, 10, 10)), affine=affine)
        t1_path = temp_dir / "t1.nii.gz"
        nib.save(t1_img, t1_path)

        atlas_dir = temp_dir / "tract_atlas"
        atlas_dir.mkdir(parents=True)
        img = nib.Nifti1Image(np.full((10, 10, 10), 0.5), affine=affine)
        nib.save(img, atlas_dir / "tract_a.nii.gz")
        nib.save(img, atlas_dir / "tract_b.nii.gz")

        args = argparse.Namespace(
            subject="test_subj",
            t1_path=str(t1_path),
            lesion_path="",
            tracts_path="",
            tract_atlas_path=str(atlas_dir / "*.nii.gz"),
            qc_report="",
            tract_mode="atlas",
            skip_skullstrip=True,
            t1_brain_path="",
            skeleton_path="",
            maps="",
        )

        bundle = _load_inputs(args)

        assert "tract_a" in bundle.tract_arrays
        assert "tract_b" in bundle.tract_arrays
        assert bundle.tract_source["tract_a"] == "atlas"
        assert bundle.tract_source["tract_b"] == "atlas"

    def test_atlas_mode_raises_when_no_atlas_files(self, tmp_path):
        affine = np.eye(4)
        t1_path = tmp_path / "t1.nii.gz"
        nib.save(nib.Nifti1Image(np.zeros((10, 10, 10)), affine=affine), t1_path)

        args = argparse.Namespace(
            subject="test_subj",
            t1_path=str(t1_path),
            lesion_path="",
            tracts_path="",
            tract_atlas_path=str(tmp_path / "nonexistent" / "*.nii.gz"),
            qc_report="",
            tract_mode="atlas",
            skip_skullstrip=True,
            t1_brain_path="",
            skeleton_path="",
            maps="",
        )

        with pytest.raises(FileNotFoundError, match="Tract mode is 'atlas', but no tract atlas files were found"):
            _load_inputs(args)


