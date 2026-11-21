"""Tests for the per-dataset memory request of the metrics rule (leukoquant/utils/metrics_utils.py)."""
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from leukoquant.utils import metrics_utils as mu

GB = 1024
WORKFLOW = Path(__file__).resolve().parents[1] / "workflow" / "workflows" / "metrics_workflow.smk"


def write_header_only_t1(path, shape, zooms):
    """A NIfTI file that holds only a header: enough for nibabel to report the shape and zooms."""
    header = nib.Nifti1Header()
    header.set_data_shape(shape)
    header.set_zooms(zooms)
    header.set_data_dtype(np.int8)
    with open(path, "wb") as fh:
        header.write_to(fh)
        fh.write(b"\0\0\0\0")  # empty header extension, so the data would start at the usual offset
    return str(path)


# The four sample cases measured by the skull-strip memory test.
CASES = [
    # name, shape, zooms, padded millions of 1 mm voxels, expected request in GB
    ("DPUK 256^3 at 1 mm", (256, 256, 256), (1.0, 1.0, 1.0), 16.78, 12),
    ("OASIS-3 176x240x256 at 1.2/1.05/1.05", (176, 240, 256), (1.2, 1.05, 1.05), 20.97, 14),
    ("ADNI-3 196x256x256 at 1.2/1.05/1.05", (196, 256, 256), (1.2, 1.05, 1.05), 26.21, 14),
    ("EPAD 196x256x256 at 1.2/1.05/1.05", (196, 256, 256), (1.2, 1.05, 1.05), 26.21, 14),
]


@pytest.mark.parametrize("name,shape,zooms,mvox,gb", CASES, ids=[c[0] for c in CASES])
def test_sample_cases(tmp_path, name, shape, zooms, mvox, gb):
    assert mu.padded_mvoxels(shape, zooms) == pytest.approx(mvox, abs=0.01)
    t1 = write_header_only_t1(tmp_path / "t1.nii", shape, zooms)
    assert mu.compute_metrics_mem_mb([t1]) == gb * GB


def test_small_image_is_held_at_the_floor(tmp_path):
    t1 = write_header_only_t1(tmp_path / "small.nii", (100, 100, 100), (1.0, 1.0, 1.0))
    assert mu.compute_metrics_mem_mb(t1) == mu.METRICS_MEM_FLOOR_MB  # a single path is accepted too


def test_huge_image_is_held_at_the_cap(tmp_path):
    t1 = write_header_only_t1(tmp_path / "huge.nii", (512, 512, 512), (2.0, 2.0, 2.0))
    assert mu.compute_metrics_mem_mb([t1]) == mu.METRICS_MEM_CAP_MB


def test_request_is_sized_for_the_largest_t1_in_whole_gb(tmp_path):
    small = write_header_only_t1(tmp_path / "a.nii", (176, 240, 256), (1.2, 1.05, 1.05))
    large = write_header_only_t1(tmp_path / "b.nii", (196, 256, 256), (1.2, 1.05, 1.05))
    request = mu.compute_metrics_mem_mb([small, large, small])
    assert request == 14 * GB  # the larger T1 (26.2 M padded voxels): the verified 14 GB
    assert request % 1024 == 0


def test_size_beyond_the_verified_range_uses_the_model_plus_margin(tmp_path):
    """Above the largest size seen (26.3 M) nothing is verified: the model's upper bound plus the margin applies."""
    huge = write_header_only_t1(tmp_path / "huge.nii", (320, 320, 320), (1.0, 1.0, 1.0))  # 32.8 M padded voxels
    assert mu.compute_metrics_mem_mb([huge]) == 19 * GB  # 4.0 + 0.40 * 32.8 = 17.1 -> 17, + 2 GB


def test_verified_request_boundaries(tmp_path):
    """Below 20 M the floor applies, from 20 M the verified 14 GB (12 GB failed for the 21 M sessions)."""
    below = write_header_only_t1(tmp_path / "a.nii", (256, 256, 256), (1.0, 1.0, 1.0))      # 16.8 M
    above = write_header_only_t1(tmp_path / "b.nii", (176, 240, 256), (1.2, 1.05, 1.05))    # 21.0 M
    assert mu.compute_metrics_mem_mb([below]) == 12 * GB
    assert mu.compute_metrics_mem_mb([above]) == 14 * GB


def test_missing_and_unreadable_paths_are_skipped_with_a_warning(tmp_path, capsys):
    good = write_header_only_t1(tmp_path / "good.nii", (196, 256, 256), (1.2, 1.05, 1.05))
    garbage = tmp_path / "garbage.nii"
    garbage.write_bytes(b"not a nifti file")
    request = mu.compute_metrics_mem_mb([str(tmp_path / "missing.nii"), str(garbage), good, None, ""])
    assert request == 14 * GB
    assert "skipped" in capsys.readouterr().out


def test_no_readable_header_falls_back_to_the_largest_known_request(tmp_path, capsys):
    """A path problem must not silently under-request memory: the fallback covers the largest known T1."""
    assert mu.compute_metrics_mem_mb([str(tmp_path / "missing.nii")]) == mu.METRICS_MEM_FALLBACK_MB
    assert mu.compute_metrics_mem_mb([]) == mu.METRICS_MEM_FALLBACK_MB
    assert mu.compute_metrics_mem_mb(None) == mu.METRICS_MEM_FALLBACK_MB
    assert mu.METRICS_MEM_FALLBACK_MB >= 16 * 1024
    assert "fallback" in capsys.readouterr().out

def test_resolve_t1_candidates(tmp_path):
    session = tmp_path / "sub-01"
    session.mkdir()
    t1 = write_header_only_t1(session / "T1.nii", (176, 240, 256), (1.2, 1.05, 1.05))
    assert mu.resolve_t1_candidates(t1, "sub-01") == [t1]                                # plain file path
    assert mu.resolve_t1_candidates(f"{tmp_path}:T1*.nii", "sub-01") == [t1]             # base:glob form
    assert mu.resolve_t1_candidates(f"{tmp_path}:T1*.nii", "sub-02") == []               # no such subject
    assert mu.resolve_t1_candidates(str(tmp_path / "nope.nii"), "sub-01") == []
    assert mu.resolve_t1_candidates(None, "sub-01") == []
    assert mu.resolve_t1_candidates(str(tmp_path / "{subject}" / "T1.nii"), "sub-01") == [t1]  # placeholder


def test_workflow_uses_the_computed_request():
    text = WORKFLOW.read_text()
    assert "mem_mb=METRICS_MEM_MB" in text
    assert "mem_mb=16 * 1024" not in text
    assert "compute_metrics_mem_mb" in text


def test_override_replaces_the_computed_request(capsys):
    assert mu.apply_metrics_mem_override(16 * GB, {"LEUKOQUANT_METRICS_MEM_GB": "12"}) == 12 * GB
    assert "override 12 GB (computed 16 GB)" in capsys.readouterr().out


def test_no_override_keeps_the_computed_request():
    assert mu.apply_metrics_mem_override(16 * GB, {}) == 16 * GB
    assert mu.apply_metrics_mem_override(16 * GB, {"LEUKOQUANT_METRICS_MEM_GB": ""}) == 16 * GB


@pytest.mark.parametrize("bad", ["abc", "7", "33", "12.5", "-4", "0"])
def test_invalid_override_is_ignored_with_a_warning(bad, capsys):
    assert mu.apply_metrics_mem_override(16 * GB, {"LEUKOQUANT_METRICS_MEM_GB": bad}) == 16 * GB
    assert "ignoring" in capsys.readouterr().out


def test_override_limits_are_inclusive():
    assert mu.apply_metrics_mem_override(16 * GB, {"LEUKOQUANT_METRICS_MEM_GB": "8"}) == 8 * GB
    assert mu.apply_metrics_mem_override(16 * GB, {"LEUKOQUANT_METRICS_MEM_GB": "32"}) == 32 * GB

