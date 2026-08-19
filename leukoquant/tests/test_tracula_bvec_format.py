"""Regression test for the trac-preproc bvec-format bug.

TRACULA's trac-preproc reads bvecs as one row per DWI volume (N rows x 3
cols) and bvals as N rows of one value each. But tracula_workflow.smk
prepares its DWI/bvec/bval inputs via dwi_utils.prepare_dwi_input() -- the
same function dti/noddi use -- which always normalizes to FSL format
(3 rows x N cols / 1 row x N cols) instead. Left uncorrected, trac-preproc
reads the FSL bvecs file's 3 axis-rows as if they were 3 gradient vectors
regardless of the real direction count, and fails with "Found N b-values
but 3 gradient vectors" for any subject with more than 3 real directions --
confirmed in production logs across OASIS3/ADNI3/EPAD/DPUK.

write_dmrirc() (tracula_utils.py) now converts back to TRACULA format via
normalize_bvec_bval_to_tracula() right before referencing the files in the
dmrirc template. These tests verify the conversion function directly, and
that write_dmrirc() actually applies it to every scan's bvec/bval pair.
"""
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from leukoquant.utils.dwi_utils import (
    normalize_bvec_bval_to_fsl,
    normalize_bvec_bval_to_tracula,
)
from leukoquant.utils.tracula_utils import write_dmrirc


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


def _write_tracula_native_bvec_bval(temp_dir, n_directions=24):
    """N rows x 3 cols bvecs, N rows x 1 col bvals -- TRACULA's own convention."""
    rng = np.random.default_rng(0)
    bvecs = rng.normal(size=(n_directions, 3))
    bvecs /= np.linalg.norm(bvecs, axis=1, keepdims=True)
    bvals = np.full((n_directions, 1), 1000.0)
    bvals[0] = 0  # one b0

    bvecs_path = temp_dir / "test.bvecs"
    bvals_path = temp_dir / "test.bvals"
    np.savetxt(bvecs_path, bvecs, fmt="%.8f")
    np.savetxt(bvals_path, bvals, fmt="%g")
    return str(bvecs_path), str(bvals_path), bvecs


class TestNormalizeBvecBvalToTracula:
    def test_converts_fsl_format_to_tracula_format(self, temp_dir):
        """The exact real-world failure: dwi_utils.prepare_dwi_input()'s FSL
        normalization produces a 3-row bvecs file for any direction count,
        which a naive line-count reader (trac-preproc) misreads as "3
        gradient vectors" -- reproduced here with 24 real directions, the
        same count seen in production logs."""
        n = 24
        bvecs_path, bvals_path, original_bvecs = _write_tracula_native_bvec_bval(
            temp_dir, n_directions=n
        )

        # Step 1: what dti/noddi's prep (and, today, tracula_workflow.smk's
        # DWI-prep loop) already does -- normalize to FSL format.
        normalize_bvec_bval_to_fsl(bvecs_path, bvals_path)
        fsl_bvecs = np.loadtxt(bvecs_path)
        assert fsl_bvecs.shape == (3, n), (
            "sanity check: FSL format is 3 rows x N cols -- a line-count "
            "reader would see exactly 3 'vectors' here, reproducing the bug"
        )

        # Step 2: the fix -- write_dmrirc() now converts back before use.
        normalize_bvec_bval_to_tracula(bvecs_path, bvals_path)
        tracula_bvecs = np.loadtxt(bvecs_path)
        tracula_bvals = np.loadtxt(bvals_path)

        assert tracula_bvecs.shape == (n, 3), (
            f"expected TRACULA-native ({n}, 3) shape, got {tracula_bvecs.shape} "
            "-- trac-preproc would still misread the direction count"
        )
        assert tracula_bvals.reshape(-1).shape == (n,)
        assert np.allclose(tracula_bvecs, original_bvecs, atol=1e-6), (
            "converted content should round-trip back to the original directions"
        )

    def test_idempotent_on_already_tracula_format(self, temp_dir):
        bvecs_path, bvals_path, original_bvecs = _write_tracula_native_bvec_bval(temp_dir)
        normalize_bvec_bval_to_tracula(bvecs_path, bvals_path)
        result = np.loadtxt(bvecs_path)
        assert result.shape == (24, 3)
        assert np.allclose(result, original_bvecs, atol=1e-6)

    def test_ambiguous_3x3_left_unchanged(self, temp_dir):
        bvecs_path, bvals_path, _ = _write_tracula_native_bvec_bval(temp_dir, n_directions=3)
        before = np.loadtxt(bvecs_path)
        normalize_bvec_bval_to_tracula(bvecs_path, bvals_path)
        after = np.loadtxt(bvecs_path)
        assert np.allclose(before, after), "3x3 is ambiguous -- must be left unchanged, not guessed at"


class TestWriteDmrircAppliesTraculaFormat:
    def test_write_dmrirc_converts_fsl_format_bvecs_before_use(self, temp_dir):
        """End-to-end: write_dmrirc() must leave the referenced bvec/bval
        files in TRACULA-native format on disk, regardless of what format
        they arrived in (matching what tracula_workflow.smk's upstream
        FSL-normalizing DWI-prep step actually hands it today)."""
        n = 24
        bvecs_path, bvals_path, original_bvecs = _write_tracula_native_bvec_bval(
            temp_dir, n_directions=n
        )
        # Simulate the files arriving already FSL-normalized, exactly as
        # dwi_utils.prepare_dwi_input() leaves them.
        normalize_bvec_bval_to_fsl(bvecs_path, bvals_path)
        assert np.loadtxt(bvecs_path).shape == (3, n)

        dwi_path = temp_dir / "dwi.nii.gz"
        import nibabel as nib
        nib.save(nib.Nifti1Image(np.zeros((4, 4, 4, n)), np.eye(4)), dwi_path)

        with patch("leukoquant.utils.tracula_utils._extract_dwi_metadata") as mock_meta:
            mock_meta.return_value = {"pe_dir": "", "epi_factor": "", "echo_spacing": ""}
            write_dmrirc(
                subject="sub-01",
                subjects_dir=str(temp_dir),
                outdir=str(temp_dir / "out"),
                dwi=str(dwi_path),
                bvecs=bvecs_path,
                bvals=bvals_path,
                doeddy=1,
            )

        final_bvecs = np.loadtxt(bvecs_path)
        final_bvals = np.loadtxt(bvals_path)
        assert final_bvecs.shape == (n, 3), (
            f"write_dmrirc() left bvecs in shape {final_bvecs.shape} -- trac-preproc "
            f"would misread this as {final_bvecs.shape[0]} gradient vectors"
        )
        assert final_bvals.reshape(-1).shape == (n,)
        assert np.allclose(final_bvecs, original_bvecs, atol=1e-6)
