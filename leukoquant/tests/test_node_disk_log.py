"""TRACULA logs the node and its free disk space around each stage, so a kill for lack of node scratch can be proven.

Background: bedpostx_post of some 256x256 DWI runs was SIGKILLed with 3-5 GB of scratch in use and no record of the free
space on the node. log_node_disk_state writes one block per call: label, node name, df of the scratch and tmp
directories and the size of the job's scratch folder."""
import subprocess
from pathlib import Path

BASH_UTILS = Path(__file__).resolve().parents[1] / "utils" / "bash_utils.sh"
SMK = Path(__file__).resolve().parents[1] / "workflow" / "workflows" / "tracula_workflow.smk"


def run(snippet, env_extra=None):
    import os
    env = os.environ.copy()
    env.update(env_extra or {})
    return subprocess.run(["bash", "-c", f"source '{BASH_UTILS}'; {snippet}"], capture_output=True, text=True, env=env)


def test_block_has_label_node_df_and_scratch_size(tmp_path):
    (tmp_path / "f").write_bytes(b"x" * 2_000_000)
    res = run(f"log_node_disk_state 'before bedpostx' '{tmp_path}'", {"TMPDIR": str(tmp_path)})
    assert res.returncode == 0
    out = res.stdout
    assert "[node] before bedpostx" in out
    assert "host:" in out
    assert "df " in out or "Filesystem" in out
    assert str(tmp_path) in out
    assert "scratch folder size" in out


def test_never_fails_with_missing_paths():
    res = run("log_node_disk_state 'x' /nonexistent/scratch", {"TMPDIR": "/nonexistent/tmp"})
    assert res.returncode == 0
    assert "[node] x" in res.stdout


def test_workflow_logs_before_each_stage_and_after_bedpostx():
    smk = SMK.read_text()
    for label in ("start", "after prep", "before bedpostx", "after bedpostx", "after path"):
        assert f'log_node_disk_state "{label}"' in smk, label


def test_node_state_is_also_logged_at_job_end():
    assert 'log_node_disk_state "at job end" "$scratch_path"; scratch_cleanup "$scratch_path"' in SMK.read_text()
