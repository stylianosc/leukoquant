"""A single b0 volume (e.g. ADNI b0_AP / b0_PA) is accepted by the DWI preparation only when asked for.

TRACULA's reverse-polarity correction (dob0 = 2) needs these scans next to the diffusion series, but every other
consumer of prepare_dwi_input needs real diffusion data, so the default must keep rejecting them.
"""
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from leukoquant.utils.dwi_utils import InputDataError, prepare_dwi_input


def _out(folder: Path) -> Path:
    """The preparation writes into an existing output folder."""
    out = folder / "out"
    out.mkdir(exist_ok=True)
    return out


def _b0_scan(folder: Path) -> Path:
    dwi = folder / "b0.nii.gz"
    nib.save(nib.Nifti1Image(np.ones((6, 6, 6), dtype=np.float32), np.eye(4)), str(dwi))
    (folder / "b0.bval").write_text("0\n")
    (folder / "b0.bvec").write_text("0\n0\n0\n")
    return dwi


def test_single_b0_volume_is_rejected_by_default(tmp_path):
    with pytest.raises(InputDataError):
        prepare_dwi_input(str(_b0_scan(tmp_path)), output_dir=str(_out(tmp_path)))


def test_single_b0_volume_is_accepted_and_made_4d_when_allowed(tmp_path):
    result = prepare_dwi_input(str(_b0_scan(tmp_path)), output_dir=str(_out(tmp_path)), allow_b0_only=True)
    assert nib.load(result["dwi"]).shape == (6, 6, 6, 1)
    assert np.loadtxt(result["bvals"]).reshape(-1).tolist() == [0.0]


def test_allow_b0_only_does_not_excuse_a_diffusion_scan_with_too_few_directions(tmp_path):
    dwi = tmp_path / "few.nii.gz"
    nib.save(nib.Nifti1Image(np.ones((6, 6, 6, 3), dtype=np.float32), np.eye(4)), str(dwi))
    (tmp_path / "few.bval").write_text("0 1000 1000\n")
    (tmp_path / "few.bvec").write_text("0 1 0\n0 0 1\n0 0 0\n")
    with pytest.raises(InputDataError):
        prepare_dwi_input(str(dwi), output_dir=str(_out(tmp_path)), allow_b0_only=True)


def test_b0_only_failure_leaves_no_marker_when_allowed(tmp_path):
    session = tmp_path / "session"
    prepare_dwi_input(str(_b0_scan(tmp_path)), output_dir=str(_out(tmp_path)), session_dir=str(session), allow_b0_only=True)
    assert not (session / "input_checks" / "dwi_failed.txt").exists()


# ---- TRACULA format of a one-volume scan ---------------------------------------------------------------------
from leukoquant.utils.dwi_utils import normalize_bvec_bval_to_tracula


@pytest.mark.parametrize("content", ["0\n0\n0\n", "0 0 0\n"], ids=["fsl_3x1", "tracula_1x3"])
def test_single_volume_bvec_is_written_as_one_row_of_three(tmp_path, content):
    # np.loadtxt gives shape (3,) for both layouts, so the old code skipped the file and left the 3 axis-rows in
    # place; trac-preproc then reads them as 3 gradient vectors for 1 b-value.
    bvec, bval = tmp_path / "bvecs", tmp_path / "bvals"
    bvec.write_text(content)
    bval.write_text("0\n")
    normalize_bvec_bval_to_tracula(str(bvec), str(bval))
    assert np.loadtxt(bvec, ndmin=2).shape == (1, 3)
    assert np.loadtxt(bval, ndmin=2).shape == (1, 1)
