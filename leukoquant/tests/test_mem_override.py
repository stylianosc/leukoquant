"""LEUKOQUANT_TRACULA_MEM_GB: the TRACULA memory request is 8 GB unless a valid override is set."""
from pathlib import Path

from leukoquant.utils import mem_override as mo

GB = 1024


def test_default_is_8_gb_when_unset_or_empty():
    assert mo.tracula_mem_mb({}) == 8 * GB
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": ""}) == 8 * GB
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": "  "}) == 8 * GB


def test_valid_override_sets_whole_gb():
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": "16"}) == 16 * GB
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": " 12 "}) == 12 * GB
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": "4"}) == 4 * GB       # lower bound
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": "24"}) == 24 * GB     # upper bound (user cap)


def test_invalid_override_is_ignored_with_a_warning(capsys):
    for bad in ("3", "25", "64", "-8", "8.5", "lots", "16GB"):
        assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": bad}) == 8 * GB, bad
    assert "ignoring LEUKOQUANT_TRACULA_MEM_GB" in capsys.readouterr().out


# ---- LEUKOQUANT_TRACULA_SCRATCH_GB: node-local scratch request of the TRACULA job -------------------------------------
# bedpostx_post of a 256x256 DWI merges per-slice sample files next to the originals, so its scratch need is several
# times that of a 116x116 DWI. The request is a scheduling number only (SGE tscratch); the default stays 3 GB.

def test_scratch_default_is_3_gb_when_unset_or_empty():
    assert mo.tracula_scratch_mb(None, {}) == 3 * GB
    assert mo.tracula_scratch_mb(None, {"LEUKOQUANT_TRACULA_SCRATCH_GB": ""}) == 3 * GB
    assert mo.tracula_scratch_mb(None, {"LEUKOQUANT_TRACULA_SCRATCH_GB": " "}) == 3 * GB


def test_scratch_override_accepts_whole_gb_in_range():
    assert mo.tracula_scratch_mb(None, {"LEUKOQUANT_TRACULA_SCRATCH_GB": "15"}) == 15 * GB
    assert mo.tracula_scratch_mb(None, {"LEUKOQUANT_TRACULA_SCRATCH_GB": " 10 "}) == 10 * GB
    assert mo.tracula_scratch_mb(None, {"LEUKOQUANT_TRACULA_SCRATCH_GB": "3"}) == 3 * GB      # lower bound
    assert mo.tracula_scratch_mb(None, {"LEUKOQUANT_TRACULA_SCRATCH_GB": "40"}) == 40 * GB    # upper bound


def test_scratch_override_rejects_bad_values_and_keeps_default(capsys):
    for bad in ["2", "41", "0", "-5", "abc", "10.5"]:
        assert mo.tracula_scratch_mb(None, {"LEUKOQUANT_TRACULA_SCRATCH_GB": bad}) == 3 * GB, bad
    assert "ignoring LEUKOQUANT_TRACULA_SCRATCH_GB" in capsys.readouterr().out


def test_scratch_override_does_not_change_memory_request():
    env = {"LEUKOQUANT_TRACULA_SCRATCH_GB": "20"}
    assert mo.tracula_mem_mb(env) == 8 * GB


def test_workflow_requests_scratch_per_subject_from_dwi():
    smk = (Path(__file__).resolve().parents[1] / "workflow" / "workflows" / "tracula_workflow.smk").read_text()
    assert "scratch_size=lambda wildcards: tracula_scratch_mb(dwi_map[wildcards.subject])" in smk
    assert "scratch_size=3*1024" not in smk


# ---- scratch estimated from the DWI size (megavoxel-volumes = voxels x volumes / 1e6) ----------------------------------
# Measured on 1828 completed ADNI-3 runs: 116x116x80x55 (59 Mvv) 1.1-1.9 GB, 256x256x80x36 (189 Mvv) 3.2-4.4 GB,
# 256x256x80x54 (283 Mvv) 3.7-4.9 GB. Request = (1.0 + 0.020 x Mvv) x 1.25, whole GB, never below 3 GB, at most 40 GB.
import gzip
import struct
import zipfile

import pytest


def _write_nifti_gz(path, dims):
    hdr = bytearray(348)
    struct.pack_into("<i", hdr, 0, 348)
    struct.pack_into("<8h", hdr, 40, len(dims), *dims, *([1] * (7 - len(dims))))
    with gzip.open(path, "wb") as f:
        f.write(bytes(hdr))


def test_estimate_from_nifti_header(tmp_path):
    p = tmp_path / "dwi.nii.gz"
    _write_nifti_gz(p, (256, 256, 80, 54))
    assert mo.dwi_megavoxel_volumes([str(p)]) == pytest.approx(256 * 256 * 80 * 54 / 1e6)


def test_estimate_from_zip_uses_uncompressed_16bit_size(tmp_path):
    p = tmp_path / "I1.zip"
    with zipfile.ZipFile(p, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("a.dcm", b"\0" * 2_000_000)
        z.writestr("b.dcm", b"\0" * 2_000_000)
    assert mo.dwi_megavoxel_volumes([str(p)]) == pytest.approx(2.0)


def test_estimate_unreadable_or_unknown_is_none(tmp_path):
    assert mo.dwi_megavoxel_volumes([]) is None
    assert mo.dwi_megavoxel_volumes([str(tmp_path / "missing.zip")]) is None
    assert mo.dwi_megavoxel_volumes([str(tmp_path / "x.txt")]) is None


def test_scratch_is_the_fitted_line_plus_margin_rounded_up(tmp_path):
    import math
    from leukoquant.utils import tracula_resource_model as rm
    big = tmp_path / "b.nii.gz"
    _write_nifti_gz(big, (256, 256, 80, 54))
    mvv = 256 * 256 * 80 * 54 / 1e6
    fit = rm.load_model()["scratch_gb"]
    # the fit is on scratch AT CLEANUP; the peak during the run is higher, hence the provisional peak factor
    expect = math.ceil(mo.TRACULA_SCRATCH_PEAK_FACTOR * (fit["intercept"] + fit["slope"] * mvv + fit["margin"]))
    assert mo.tracula_scratch_mb([str(big)], {}) == expect * GB


def test_peak_factor_is_a_documented_provisional_two():
    assert mo.TRACULA_SCRATCH_PEAK_FACTOR == 2.0
    assert "provisional" in mo.__doc__.lower()


def test_scratch_grows_with_dwi_size_and_has_no_size_classes(tmp_path):
    sizes = []
    for k, dims in enumerate([(116, 116, 80, 55), (128, 128, 80, 100), (256, 256, 80, 36), (256, 256, 80, 54), (256, 256, 80, 130)]):
        p = tmp_path / f"d{k}.nii.gz"
        _write_nifti_gz(p, dims)
        sizes.append(mo.tracula_scratch_mb([str(p)], {}))
    assert sizes == sorted(sizes)
    assert sizes[-1] > sizes[0]


def test_scratch_has_headroom_over_the_measured_256_runs(tmp_path):
    # completed ADNI-3 runs: 256x256x80x36 held up to 4.35 GB, 256x256x80x54 up to 4.88 GB when the run ended; the peak
    # during the run is unknown (higher), so at least twice the largest cleanup value is requested (about 10 GB)
    for vols, used_gb in ((36, 4.35), (54, 4.88)):
        p = tmp_path / f"v{vols}.nii.gz"
        _write_nifti_gz(p, (256, 256, 80, vols))
        assert mo.tracula_scratch_mb([str(p)], {}) >= 2 * used_gb * GB * 0.99, vols


def test_scratch_override_wins_over_the_estimate(tmp_path):
    big = tmp_path / "b.nii.gz"
    _write_nifti_gz(big, (256, 256, 80, 54))
    assert mo.tracula_scratch_mb([str(big)], {"LEUKOQUANT_TRACULA_SCRATCH_GB": "15"}) == 15 * GB


def test_scratch_is_clamped_only_at_the_ends(tmp_path):
    huge = tmp_path / "h.nii.gz"
    _write_nifti_gz(huge, (512, 512, 200, 300))
    assert mo.tracula_scratch_mb([str(huge)], {}) == 40 * GB          # upper end
    assert mo.tracula_scratch_mb(None, {}) == 3 * GB                  # nothing known: lower end
    assert mo.tracula_scratch_mb([str(tmp_path / "missing.zip")], {}) == 3 * GB


# ---- memory request from the DWI size: fitted line, clamped to 8-16 GB, override up to the 24 GB cap ------------------

def test_memory_without_dwi_is_the_8_gb_floor_as_before():
    assert mo.tracula_mem_mb({}) == 8 * GB
    assert mo.tracula_mem_mb({}, dwi_paths=None) == 8 * GB


def test_memory_is_the_peak_stage_prediction_clamped_to_8_16_gb(tmp_path):
    from leukoquant.utils import tracula_resource_model as rm
    model = rm.load_model()
    for dims in [(116, 116, 80, 55), (256, 256, 80, 54), (256, 256, 80, 300), (512, 512, 200, 300)]:
        p = tmp_path / ("x%d_%d.nii.gz" % (dims[0], dims[3]))
        _write_nifti_gz(p, dims)
        mvv = dims[0] * dims[1] * dims[2] * dims[3] / 1e6
        peak = max(model[k]["intercept"] + model[k]["slope"] * mvv + model[k]["margin"] for k in ("prep_rss_gb", "path_rss_gb"))
        import math
        expect = min(max(math.ceil(peak), 8), 16)
        assert mo.tracula_mem_mb({}, dwi_paths=[str(p)]) == expect * GB, dims


def test_memory_override_wins_and_never_exceeds_24_gb(tmp_path):
    p = tmp_path / "x.nii.gz"
    _write_nifti_gz(p, (256, 256, 80, 54))
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": "24"}, dwi_paths=[str(p)]) == 24 * GB
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": "30"}, dwi_paths=[str(p)]) <= 16 * GB


def test_workflow_requests_memory_per_subject_from_dwi():
    smk = (Path(__file__).resolve().parents[1] / "workflow" / "workflows" / "tracula_workflow.smk").read_text()
    assert "mem_mb=lambda wildcards: tracula_mem_mb(dwi_paths=dwi_map[wildcards.subject])" in smk


def test_memory_request_exceeds_the_direct_large_image_evidence(tmp_path):
    """The only direct evidence about 256x256 memory: trac-paths failed with bad_alloc at the 8 GB cap for 256x256x80x36
    sessions (peak RSS 7996 MB = the cap, so the need is above 8 GB), and one 256x256x80x54 run completed with a 10.9 GB
    peak. A model that requests 8 GB for the first or less than 10.9 GB for the second repeats a known failure."""
    for vols, floor_mb in ((36, 8 * GB + 1), (54, 10893)):      # the need is above 8 GB; a 10.9 GB peak was measured
        p = tmp_path / f"e{vols}.nii.gz"
        _write_nifti_gz(p, (256, 256, 80, vols))
        request = mo.tracula_mem_mb({}, dwi_paths=[str(p)])
        assert floor_mb <= request <= 16 * GB, (vols, request)


def test_scratch_factor_only_applies_to_a_cleanup_based_fit(monkeypatch):
    """A fit on true sampled peaks ("scratch_basis": "peak") must not be multiplied by the provisional factor."""
    from leukoquant.utils import mem_override as mo
    fit = {"intercept": 0.0, "slope": 0.0, "margin": 5.0}
    monkeypatch.setattr(mo, "dwi_megavoxel_volumes", lambda p: 100.0)
    monkeypatch.setattr(mo, "load_model", lambda: {"scratch_gb": fit})
    assert mo.tracula_scratch_mb(["x"], environ={}) == 10 * 1024
    monkeypatch.setattr(mo, "load_model", lambda: {"scratch_gb": fit, "scratch_basis": "peak"})
    assert mo.tracula_scratch_mb(["x"], environ={}) == 5 * 1024
