"""The --topup option (none|auto) must reach DWI selection and the dmrirc writer in both TRACULA entry points."""
from pathlib import Path

import pytest
from click.testing import CliRunner

from leukoquant.cli import commands
from leukoquant.cli.main import main
from leukoquant.core import tracula_processor
from leukoquant.utils.duplicate_dwi import Selection

WORKFLOW = Path(__file__).resolve().parents[1] / "workflow" / "workflows" / "tracula_workflow.smk"


@pytest.mark.parametrize("command", ["process-tracula", "process-all"])
def test_cli_lists_topup_option_with_both_choices(command):
    result = CliRunner().invoke(main, [command, "--help"])
    assert result.exit_code == 0
    assert "--topup" in result.output and "none" in result.output and "auto" in result.output


@pytest.mark.parametrize("command", ["process-tracula", "process-all"])
def test_cli_rejects_unknown_topup_value(command):
    result = CliRunner().invoke(main, [command, "--topup", "sometimes"])
    assert result.exit_code != 0


def test_discovery_forwards_topup_to_the_selection(monkeypatch):
    seen = {}

    def fake_select(candidates, *args, **kw):
        seen.update(kw)
        return Selection([])

    monkeypatch.setattr(tracula_processor, "resolve_subject_pattern", lambda *a, **k: ["x.zip"])
    monkeypatch.setattr(tracula_processor, "select_scans", fake_select)
    with pytest.raises(FileNotFoundError):
        tracula_processor._discover_dwi_with_sidecars("*.zip", ["s1"], None, None, topup="auto")
    assert seen["topup"] == "auto"


def test_process_tracula_passes_topup_on(monkeypatch):
    seen = {}
    monkeypatch.setattr(commands, "apply_tracula", lambda **kw: seen.update(kw) or {"success": True})
    commands.process_tracula(subject="s", dwi="d", output_dir="o", topup="auto")
    assert seen["topup"] == "auto"


def test_process_all_passes_topup_on(monkeypatch):
    seen = {}
    monkeypatch.setattr(commands, "apply_process_all", lambda **kw: seen.update(kw) or {"success": True})
    commands.process_all(subject="s", t1="t", flair="f", dwi="d", output_dir="o", topup="auto")
    assert seen["topup"] == "auto"


def test_workflow_hands_topup_to_the_dmrirc_writer_instead_of_a_fixed_dob0():
    text = WORKFLOW.read_text()
    assert "--topup {params.topup}" in text
    assert "--dob0 0" not in text
    assert 'config.get("topup", "none")' in text


def test_workflow_lets_the_dwi_preparation_accept_b0_only_scans_only_with_topup_auto():
    text = WORKFLOW.read_text()
    assert '[ "{params.topup}" = "auto" ]' in text
    assert '--allow-b0-only' in text
    # the flag is handed to the per-scan preparation call of dwi_utils
    prep_call = text[text.index("/leukoquant/leukoquant/utils/dwi_utils.py"):][:400]
    assert "$b0_only_flag" in prep_call


def test_workflow_points_trac_preproc_at_this_containers_topup_config():
    # The container keeps b02b0.cnf in $FSLDIR/etc/flirtsch; trac-preproc looks in $FSLDIR/src/topup/flirtsch.
    text = WORKFLOW.read_text()
    assert "src/topup/flirtsch/b02b0.cnf" in text and "etc/flirtsch/b02b0.cnf" in text
    assert text.index("etc/flirtsch/b02b0.cnf") > text.index('cp "$trac_preproc_file" "$local_bin/trac-preproc"')
