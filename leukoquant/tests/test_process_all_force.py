"""process-all --force: stage names -> Snakemake rule names.

Snakemake silently ignores a --forcerun name that matches no rule, so a wrong
entry in PROCESS_ALL_FORCE_STAGES would make --force quietly do nothing. These
tests pin the resolution behaviour and check every mapped rule name against
the rule names process_all_workflow.smk actually builds (module prefixes and
the atlas converter nested inside tracula included).
"""
import importlib
import re
from pathlib import Path

import pytest
from click.testing import CliRunner

# The module, not the click group that leukoquant.cli re-exports as "main".
cli_main = importlib.import_module("leukoquant.cli.main")
from leukoquant.cli.commands import (PROCESS_ALL_FORCE_STAGES,
                                     process_all_force_rules)

WORKFLOWS = Path(__file__).resolve().parents[1] / "workflow" / "workflows"


# --------------------------------------------------------------------------
# Resolution
# --------------------------------------------------------------------------
def test_nothing_forced():
    assert process_all_force_rules(None, "gif,freesurfer") is None
    assert process_all_force_rules("", "gif,freesurfer") is None


def test_metrics_expands_per_parcellation():
    assert process_all_force_rules("metrics", "gif,freesurfer") == [
        "metrics_extract_metrics_gif", "metrics_extract_metrics_freesurfer"]


def test_several_stages_in_order_without_duplicates():
    assert process_all_force_rules("zscore, METRICS,zscore", "gif") == [
        "z_score_z_score", "metrics_extract_metrics_gif"]


def test_unknown_stage_is_an_error():
    with pytest.raises(ValueError, match="unknown stage"):
        process_all_force_rules("zscores", "gif")


def test_zscore_with_skip_zscore_is_an_error():
    with pytest.raises(ValueError, match="skip-zscore"):
        process_all_force_rules("zscore", "gif", skip_zscore=True)


def test_atlas_needs_a_non_freesurfer_parcellation():
    with pytest.raises(ValueError, match="non-FreeSurfer"):
        process_all_force_rules("atlas", "freesurfer")
    assert process_all_force_rules("atlas", "gif,freesurfer") == ["tracula_other_convert_atlas"]


# --------------------------------------------------------------------------
# Every mapped rule exists in process-all's DAG
# --------------------------------------------------------------------------
_RULE = re.compile(r"^\s*rule\s+(\w+)\s*:", re.M)
_NAMED_RULE = re.compile(r"^\s*name\s*:\s*f?\"([^\"]+)\"", re.M)
_MODULE = re.compile(r"^\s*module\s+(\w+)\s*:\s*\n\s*snakefile\s*:\s*\"([^\"]+)\"", re.M)
_USE_ALL = re.compile(r"^\s*use rule \* from (\w+) as (\w+)_\*", re.M)


def _rule_names(smk: Path, parc: str) -> set:
    """Rule names a workflow file defines, including rules imported from its
    modules under their prefix. Loop-generated names (name: f"..._{_parc}...")
    are expanded for `parc`; any other f-string field (e.g. metrics' mode tag,
    empty under process-all) is dropped."""
    text = smk.read_text()
    names = set(_RULE.findall(text))
    for n in _NAMED_RULE.findall(text):
        names.add(re.sub(r"\{[^}]*\}", "", n.replace("{_parc}", parc)))
    modules = {m: WORKFLOWS / f for m, f in _MODULE.findall(text)}
    for module, prefix in _USE_ALL.findall(text):
        names |= {f"{prefix}_{n}" for n in _rule_names(modules[module], parc)}
    return names


@pytest.mark.parametrize("stage", sorted(PROCESS_ALL_FORCE_STAGES))
def test_mapped_rules_exist_in_process_all(stage):
    for parc in ("gif", "freesurfer"):
        dag_rules = _rule_names(WORKFLOWS / "process_all_workflow.smk", parc)
        for rule in PROCESS_ALL_FORCE_STAGES[stage]:
            name = rule.format(parc=parc)
            if stage == "atlas" and parc == "freesurfer":
                continue  # atlas conversion only runs for non-FreeSurfer parcellations
            assert name in dag_rules, f"--force {stage}: no rule {name!r} in process-all"


# --------------------------------------------------------------------------
# CLI: bare --force keeps meaning "metrics"; --force STAGES passes the list
# --------------------------------------------------------------------------
@pytest.mark.parametrize("args, expected", [
    ([], None),
    (["--force"], "metrics"),
    (["--force", "--verbose"], "metrics"),
    (["--force", "zscore,metrics"], "zscore,metrics"),
])
def test_cli_force_value(monkeypatch, args, expected):
    seen = {}

    def fake_process_all(*a):
        seen["force"] = a[17]  # positional order of commands.process_all
        return {"success": True}

    monkeypatch.setattr(cli_main, "process_all", fake_process_all)
    monkeypatch.setattr(cli_main, "check_sge_plugin", lambda *_: None)
    result = CliRunner().invoke(cli_main.main, ["process-all", *args])
    assert result.exit_code == 0, result.output
    assert seen["force"] == expected


def test_healthy_cohort_is_always_queued_but_gets_no_metrics_targets():
    """Healthy reference sessions are built (normal and forced runs) but only the non-healthy subjects are metrics
    targets, so a forced metrics rerun redoes exactly the requested sessions."""
    from leukoquant.core.process_all_processor import queue_with_healthy

    requested = ["sub-1", "sub-2"]
    healthy = ["sub-h1", "sub-h2"]
    assert queue_with_healthy(requested, healthy) == requested + healthy
    assert queue_with_healthy(requested + ["sub-h1"], healthy) == requested + ["sub-h1", "sub-h2"]
    assert requested == ["sub-1", "sub-2"]  # the caller's list is not modified
    workflow = (WORKFLOWS / "process_all_workflow.smk").read_text()
    block = workflow[workflow.index("config_metrics = {"):]
    assert '"subjects": target_subjects' in block[:700]
