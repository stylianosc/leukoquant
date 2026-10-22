"""The metrics step's separate output root must be mounted in the container (metrics_output_dir).

The 2026-10-18 pilot failed on every metrics job with ``mkdir: cannot create directory '/metrics_output':
Read-only file system`` because the folder was declared for the rule but never added to the bind list.
"""
from pathlib import Path

from leukoquant.utils.bind_utils import METRICS_OUTPUT_CONTAINER_DIR, metrics_output_bind


def test_separate_root_is_mounted_and_created(tmp_path):
    out = tmp_path / "outputs"
    metrics = tmp_path / "metrics_v2_1"
    assert metrics_output_bind(str(metrics), str(out)) == [f"{metrics}:{METRICS_OUTPUT_CONTAINER_DIR}"]
    assert metrics.is_dir()  # a bind source has to exist


def test_no_entry_without_a_separate_root(tmp_path):
    assert metrics_output_bind(None, str(tmp_path)) == []
    assert metrics_output_bind("", str(tmp_path)) == []


def test_no_entry_when_it_is_the_output_dir(tmp_path):
    assert metrics_output_bind(str(tmp_path), str(tmp_path)) == []
    assert metrics_output_bind(str(tmp_path / "."), str(tmp_path)) == []


def test_container_path_matches_the_workflow():
    smk = (Path(__file__).resolve().parents[1] / "workflow" / "workflows" / "process_all_workflow.smk").read_text()
    assert f'"{METRICS_OUTPUT_CONTAINER_DIR}"' in smk
