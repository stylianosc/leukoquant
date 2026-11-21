"""A failed bedpostx stage must leave its own logs (esp. bedpostx_post.e*) in the output logs folder.

Background: bedpostx_post of 021-s-7055 was SIGKILLed after 28 h; its error file lived in node scratch, which is
deleted at job end, so the reason was lost. rescue_bedpostx_logs copies dmri.bedpostX/logs out before the cleanup."""
import subprocess
from pathlib import Path

BASH_UTILS = Path(__file__).resolve().parents[1] / "utils" / "bash_utils.sh"
SMK = Path(__file__).resolve().parents[1] / "workflow" / "workflows" / "tracula_workflow.smk"


def rescue(bedpostx_dir, dest):
    return subprocess.run(["bash", "-c", f"source '{BASH_UTILS}'; rescue_bedpostx_logs '{bedpostx_dir}' '{dest}'"],
                          capture_output=True, text=True)


def test_logs_are_copied_with_error_files(tmp_path):
    logs = tmp_path / "dmri.bedpostX" / "logs"
    (logs / "monitor").mkdir(parents=True)
    (logs / "bedpostx_post.e31473").write_text("Killed\n")
    (logs / "bedpostx_post.o31473").write_text("")
    dest = tmp_path / "out_logs"
    dest.mkdir()
    res = rescue(tmp_path / "dmri.bedpostX", dest)
    assert res.returncode == 0
    assert (dest / "bedpostx_logs" / "bedpostx_post.e31473").read_text() == "Killed\n"


def test_missing_bedpostx_dir_is_harmless(tmp_path):
    dest = tmp_path / "out_logs"
    dest.mkdir()
    res = rescue(tmp_path / "nope.bedpostX", dest)
    assert res.returncode == 0
    assert not (dest / "bedpostx_logs").exists()


def test_rerun_replaces_the_previous_copy(tmp_path):
    logs = tmp_path / "dmri.bedpostX" / "logs"
    logs.mkdir(parents=True)
    dest = tmp_path / "out_logs"
    (dest / "bedpostx_logs").mkdir(parents=True)
    (dest / "bedpostx_logs" / "old.e1").write_text("old")
    (logs / "new.e2").write_text("new")
    assert rescue(tmp_path / "dmri.bedpostX", dest).returncode == 0
    assert (dest / "bedpostx_logs" / "new.e2").exists()
    assert not (dest / "bedpostx_logs" / "old.e1").exists()


def test_workflow_rescues_logs_after_failed_bedpostx():
    smk = SMK.read_text()
    assert "rescue_bedpostx_logs" in smk
