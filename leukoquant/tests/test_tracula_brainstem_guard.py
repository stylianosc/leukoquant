"""Brainstem field-of-view recovery: coverage measure and tract exclusion."""
import nibabel as nib
import numpy as np
import pytest

from leukoquant.utils import tracula_utils as tu


def _save(path, data, affine):
    nib.save(nib.Nifti1Image(data.astype(np.int16), affine), str(path))


@pytest.fixture
def aseg(tmp_path):
    """1 mm 20^3 aseg whose brainstem (label 16) fills z = 0..9."""
    data = np.zeros((20, 20, 20))
    data[5:15, 5:15, 0:10] = tu.BRAINSTEM_LABEL
    p = tmp_path / "aseg.nii.gz"
    _save(p, data, np.eye(4))
    return p


def test_full_coverage(aseg, tmp_path):
    dwi = tmp_path / "dwi.nii.gz"
    _save(dwi, np.zeros((20, 20, 20, 3)), np.eye(4))
    assert tu.brainstem_coverage(aseg, dwi) == 1.0


def test_lower_brainstem_cut_off(aseg, tmp_path):
    """DWI grid starting at z=4 loses z=0..3 of the brainstem's z=0..9 -> 60% left."""
    dwi = tmp_path / "dwi.nii.gz"
    affine = np.eye(4)
    affine[2, 3] = 4
    _save(dwi, np.zeros((20, 20, 16, 3)), affine)
    assert tu.brainstem_coverage(aseg, dwi) == pytest.approx(0.6)


def test_coarser_dwi_grid(aseg, tmp_path):
    """Different voxel size: coverage is measured in scanner space, not voxel counts."""
    dwi = tmp_path / "dwi.nii.gz"
    _save(dwi, np.zeros((10, 10, 10, 3)), np.diag([2.0, 2.0, 2.0, 1.0]))
    assert tu.brainstem_coverage(aseg, dwi) == 1.0


def test_no_brainstem_in_aseg_means_nothing_to_protect(tmp_path):
    aseg = tmp_path / "aseg.nii.gz"
    _save(aseg, np.zeros((5, 5, 5)), np.eye(4))
    dwi = tmp_path / "dwi.nii.gz"
    _save(dwi, np.zeros((5, 5, 5, 3)), np.eye(4))
    assert tu.brainstem_coverage(aseg, dwi) == 1.0


def test_read_default_pathlist(tmp_path):
    f = tmp_path / "pathlist.txt"
    f.write_text("acomm\t0\t7\nlh.cst\t1\t8\nrh.cst\t2\t8\nmcp\t3\t6\n")
    assert tu.read_default_pathlist(f) == ["acomm", "lh.cst", "rh.cst", "mcp"]


def _pathlist(tmp_path):
    pl = tmp_path / "pathlist.txt"
    pl.write_text("acomm\t0\t7\nlh.cst\t1\t8\nrh.cst\t2\t8\nmcp\t3\t6\n")
    return pl


def test_exclude_tracts_appends_reduced_pathlist(tmp_path):
    rc = tmp_path / "dmrirc"
    rc.write_text("set dob0 = 0\n")

    kept = tu.exclude_tracts_from_dmrirc(str(rc), str(_pathlist(tmp_path)), "lh.cst")

    assert kept == ["acomm", "rh.cst", "mcp"]
    lines = rc.read_text().splitlines()
    assert lines[0] == "set dob0 = 0"          # existing settings untouched
    assert lines[-1] == "set pathlist = ( acomm rh.cst mcp )"


def test_exclusions_accumulate_across_retries(tmp_path):
    rc = tmp_path / "dmrirc"
    rc.write_text("")
    pl = _pathlist(tmp_path)
    tu.exclude_tracts_from_dmrirc(str(rc), str(pl), "lh.cst")
    kept = tu.exclude_tracts_from_dmrirc(str(rc), str(pl), "rh.cst")
    assert kept == ["acomm", "mcp"]
    assert rc.read_text().strip().splitlines()[-1] == "set pathlist = ( acomm mcp )"


def test_excluding_an_already_excluded_tract_raises(tmp_path):
    """Lets the workflow stop instead of retrying the same failure forever."""
    rc = tmp_path / "dmrirc"
    rc.write_text("")
    pl = _pathlist(tmp_path)
    tu.exclude_tracts_from_dmrirc(str(rc), str(pl), "lh.cst")
    with pytest.raises(ValueError):
        tu.exclude_tracts_from_dmrirc(str(rc), str(pl), "lh.cst")


def test_failing_tract_is_the_last_one_loaded(tmp_path):
    log = tmp_path / "tracula_all.log"
    log.write_text(
        "Loading streamlines from /usr/local/freesurfer/trctrain/hcp/mgh_1001/syn/mcp.bbr.prep.trk\n"
        "WARN: Could not find satisfactory control point fit - try 1\n"
        "Loading streamlines from /usr/local/freesurfer/trctrain/hcp/mgh_1001/syn/lh.cst.bbr.prep.trk\n"
        "Loading streamlines from /usr/local/freesurfer/trctrain/hcp/mgh_1006/syn/lh.cst.bbr.prep.trk\n"
        "WARN: Could not find satisfactory control point fit - try 17025\n")
    assert tu.failing_tract(str(log)) == "lh.cst"


def test_failing_tract_empty_when_nothing_was_loaded(tmp_path):
    log = tmp_path / "tracula_all.log"
    log.write_text("nothing relevant\n")
    assert tu.failing_tract(str(log)) == ""


def test_cli_entrypoints_used_by_the_workflow(tmp_path):
    import subprocess, sys
    log = tmp_path / "log.txt"
    log.write_text("Loading streamlines from /x/syn/rh.cst.bbr.prep.trk\n")
    rc = tmp_path / "dmrirc"
    rc.write_text("")
    pl = _pathlist(tmp_path)

    tract = subprocess.run([sys.executable, tu.__file__, "failing-tract", str(log)],
                           capture_output=True, text=True, check=True).stdout.strip()
    assert tract == "rh.cst"
    subprocess.run([sys.executable, tu.__file__, "exclude-tracts", str(rc), str(pl), tract],
                   capture_output=True, text=True, check=True)
    again = subprocess.run([sys.executable, tu.__file__, "exclude-tracts", str(rc), str(pl), tract],
                           capture_output=True, text=True)
    assert again.returncode == 3      # already excluded -> the workflow stops


def test_fov_features_full_and_cut(aseg, tmp_path):
    full = tmp_path / "full.nii.gz"
    _save(full, np.zeros((20, 20, 20, 3)), np.eye(4))
    f = tu.brainstem_fov_features(aseg, full)
    assert f["coverage"] == 1.0 and f["inferior_coverage"] == 1.0
    assert f["inferior_margin_mm"] <= 0 and f["brainstem_height_mm"] == pytest.approx(9.0)

    cut = tmp_path / "cut.nii.gz"            # grid starts at z=4: loses z=0..3 of 0..9
    affine = np.eye(4)
    affine[2, 3] = 4
    _save(cut, np.zeros((20, 20, 16, 3)), affine)
    f = tu.brainstem_fov_features(aseg, cut)
    assert f["coverage"] == pytest.approx(0.6)
    assert f["inferior_coverage"] == pytest.approx(0.0)   # lowest third (z=0..3) entirely cut
    assert f["inferior_margin_mm"] == pytest.approx(4.0)
