"""dob0 (reverse-polarity correction) in the dmrirc: switched on only for --topup auto and a real b0 pair."""
import re
from unittest.mock import patch

import nibabel as nib
import numpy as np
import pytest

from leukoquant.utils.tracula_utils import has_reverse_polarity_pair, write_dmrirc


# ---- pure decision: which scans make a reverse-polarity pair ---------------------------------------------------
def test_pair_of_single_b0_scans_with_opposite_directions():
    assert has_reverse_polarity_pair(["PA", "AP", "PA"], [[0, 1000, 1000], [0], [0]])


@pytest.mark.parametrize("a, b", [("PA", "AP"), ("AP", "PA"), ("RL", "LR"), ("IS", "SI")])
def test_all_opposite_direction_names(a, b):
    assert has_reverse_polarity_pair(["PA", a, b], [[0, 1000], [0], [0]])


def test_same_direction_is_no_pair():
    assert not has_reverse_polarity_pair(["PA", "PA", "PA"], [[0, 1000], [0], [0]])


def test_opposite_direction_on_a_diffusion_scan_is_no_pair():
    # only b0-only scans count; a diffusion scan with another direction is not a reverse-polarity acquisition
    assert not has_reverse_polarity_pair(["PA", "AP"], [[0, 1000], [0, 1000]])


def test_missing_direction_is_no_pair():
    assert not has_reverse_polarity_pair(["PA", "", "AP"], [[0, 1000], [0], [0]])
    assert not has_reverse_polarity_pair(["PA", "AP"], [[0, 1000], [0]])


def test_one_scan_is_no_pair():
    assert not has_reverse_polarity_pair(["PA"], [[0, 1000]])


# ---- written dmrirc ---------------------------------------------------------------------------------------------
def _scan(folder, name, bvals):
    """A tiny DWI with matching FSL-format bvals/bvecs; returns (dwi, bvecs, bvals)."""
    n = len(bvals)
    dwi = folder / f"{name}.nii.gz"
    nib.save(nib.Nifti1Image(np.zeros((4, 4, 4, n), dtype=np.float32), np.eye(4)), dwi)
    bval = folder / f"{name}.bval"
    bvec = folder / f"{name}.bvec"
    bval.write_text(" ".join(str(b) for b in bvals) + "\n")
    rng = np.random.default_rng(0)
    vecs = rng.normal(size=(3, n))
    vecs[:, np.array(bvals) == 0] = 0
    bvec.write_text("\n".join(" ".join(f"{x:.4f}" for x in row) for row in vecs) + "\n")
    return str(dwi), str(bvec), str(bval)


def _dmrirc(tmp_path, topup, pe_dirs, bvals_per_scan):
    scans = [_scan(tmp_path, f"data_{i}", b) for i, b in enumerate(bvals_per_scan)]
    meta = [{"pe_dir": pe, "epi_factor": "81", "echo_spacing": "0.7"} for pe in pe_dirs]
    with patch("leukoquant.utils.tracula_utils._extract_dwi_metadata", side_effect=meta):
        path, _ = write_dmrirc(
            subject="sub-01", subjects_dir=str(tmp_path), outdir=str(tmp_path / "out"),
            dwi=[s[0] for s in scans], bvecs=[s[1] for s in scans], bvals=[s[2] for s in scans],
            doeddy=1, topup=topup,
        )
    return path.read_text()


def _dob0(text):
    return int(re.search(r"^set dob0 = (\d)", text, re.M).group(1))


MAIN = [0, 0, 1000, 1000, 1000, 1000, 1000, 1000, 2000, 2000, 2000, 2000]


def test_topup_auto_with_b0_pair_sets_dob0_2(tmp_path):
    text = _dmrirc(tmp_path, "auto", ["PA", "PA", "AP"], [MAIN, [0], [0]])
    assert _dob0(text) == 2
    assert "set pedir = ( PA PA AP )" in text


def test_topup_none_never_sets_dob0(tmp_path):
    assert _dob0(_dmrirc(tmp_path, "none", ["PA", "PA", "AP"], [MAIN, [0], [0]])) == 0


def test_topup_auto_without_a_pair_keeps_dob0_0(tmp_path):
    assert _dob0(_dmrirc(tmp_path, "auto", ["PA"], [MAIN])) == 0


def test_default_is_no_topup(tmp_path):
    scans = [_scan(tmp_path, f"data_{i}", b) for i, b in enumerate([MAIN, [0], [0]])]
    meta = [{"pe_dir": pe, "epi_factor": "81", "echo_spacing": "0.7"} for pe in ("PA", "PA", "AP")]
    with patch("leukoquant.utils.tracula_utils._extract_dwi_metadata", side_effect=meta):
        path, _ = write_dmrirc(
            subject="sub-01", subjects_dir=str(tmp_path), outdir=str(tmp_path / "out"),
            dwi=[s[0] for s in scans], bvecs=[s[1] for s in scans], bvals=[s[2] for s in scans], doeddy=1,
        )
    assert _dob0(path.read_text()) == 0


def test_invalid_topup_value_is_rejected(tmp_path):
    with pytest.raises(ValueError):
        _dmrirc(tmp_path, "sometimes", ["PA"], [MAIN])
