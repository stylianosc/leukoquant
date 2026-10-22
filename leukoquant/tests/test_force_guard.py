"""The safety check for ``process-all --force`` (leukoquant/utils/force_guard.py).

A forced run must not rerun a stage upstream of the forced one: Snakemake deletes the outputs of every job
it reruns, which destroyed finished TRACULA outputs on 2026-10-05 and 2026-10-18. These tests pin the decision
on small, hand-written plans, and check that every stage the CLI can force is known to the guard.
"""
import pytest

from leukoquant.cli.commands import PROCESS_ALL_FORCE_STAGES, process_all_force_rules
from leukoquant.utils.force_guard import (STAGE_ORDER, describe_unexpected, parse_job_stats,
                                          stage_of_rule, unexpected_upstream_rules)

# What Snakemake prints for a dry run (tabulate output), trimmed to the part the parser reads.
PLAN_TEXT = """\
Building DAG of jobs...
Job stats:
job                                 count
--------------------------------  -------
all                                     1
metrics_extract_metrics_freesurfer     10
metrics_extract_metrics_gif            10
tracula_trac_all_gif                   55
total                                  76

Reasons:
    (check the reasons)
"""


class TestParseJobStats:
    def test_reads_every_rule_and_skips_total(self):
        counts = parse_job_stats(PLAN_TEXT)
        assert counts == {"all": 1, "metrics_extract_metrics_freesurfer": 10,
                          "metrics_extract_metrics_gif": 10, "tracula_trac_all_gif": 55}

    def test_empty_when_nothing_is_planned(self):
        assert parse_job_stats("Nothing to be done (all requested files are present and up to date).") == {}

    def test_last_table_wins(self):
        two_tables = PLAN_TEXT + "\nJob stats:\njob   count\n-----  -----\nall        1\ntotal      1\n"
        assert parse_job_stats(two_tables) == {"all": 1}


class TestStageOfRule:
    @pytest.mark.parametrize("rule, stage", [
        ("recon_recon_all", "recon"),
        ("gif_run_gif", "gif"),
        ("bamos_run_bamos_correction", "bamos"),
        ("tracula_other_convert_atlas", "atlas"),
        ("tracula_trac_all_freesurfer", "tracula"),
        ("dti_dti", "dti"),
        ("noddi_noddi", "noddi"),
        ("tract_qc_tract_qc_gif", "tractqc"),
        ("z_score_z_score", "zscore"),
        ("metrics_extract_metrics_gif", "metrics"),
        ("all", ""),
    ])
    def test_rule_to_stage(self, rule, stage):
        assert stage_of_rule(rule) == stage

    def test_every_forceable_rule_belongs_to_its_own_stage(self):
        # PROCESS_ALL_FORCE_STAGES owns the rule names; the guard must map each back to the same stage.
        for stage, rules in PROCESS_ALL_FORCE_STAGES.items():
            for rule in rules:
                for name in (rule.format(parc="gif"), rule.format(parc="freesurfer")):
                    assert stage_of_rule(name) == stage, (stage, name)

    def test_every_forceable_stage_is_ordered(self):
        assert set(STAGE_ORDER) == set(PROCESS_ALL_FORCE_STAGES)


class TestUnexpectedUpstreamRules:
    def test_the_2026_10_18_plan_is_refused(self):
        """Forcing metrics while TRACULA, DTI, NODDI and tract QC are also planned (the ADNI-3 incident)."""
        plan = {"metrics_extract_metrics_gif": 10, "tracula_trac_all_gif": 55, "tracula_trac_all_freesurfer": 55,
                "dti_dti": 4, "noddi_noddi": 4, "tract_qc_tract_qc_gif": 4, "all": 1}
        assert unexpected_upstream_rules(plan, ["metrics"]) == {
            "tracula_trac_all_gif": 55, "tracula_trac_all_freesurfer": 55, "dti_dti": 4, "noddi_noddi": 4,
            "tract_qc_tract_qc_gif": 4}

    def test_metrics_only_plan_is_allowed(self):
        assert unexpected_upstream_rules({"metrics_extract_metrics_gif": 10, "all": 1}, ["metrics"]) == {}

    def test_downstream_of_the_forced_stage_is_allowed(self):
        plan = {"z_score_z_score": 3, "metrics_extract_metrics_gif": 10}
        assert unexpected_upstream_rules(plan, ["zscore"]) == {}

    def test_earlier_forced_stage_allows_later_stages(self):
        plan = {"tracula_trac_all_gif": 5, "dti_dti": 5, "metrics_extract_metrics_gif": 5}
        assert unexpected_upstream_rules(plan, ["tracula"]) == {}

    def test_stage_before_the_earliest_forced_one_is_refused(self):
        plan = {"tracula_trac_all_gif": 5, "metrics_extract_metrics_gif": 5}
        assert unexpected_upstream_rules(plan, ["zscore", "metrics"]) == {"tracula_trac_all_gif": 5}

    def test_zero_counts_and_aggregation_rules_are_ignored(self):
        assert unexpected_upstream_rules({"tracula_trac_all_gif": 0, "all": 1}, ["metrics"]) == {}

    def test_no_forced_stage_means_no_check(self):
        assert unexpected_upstream_rules({"tracula_trac_all_gif": 5}, []) == {}


def test_forced_rules_resolve_to_known_stages():
    """What the processor passes to the guard (rule names -> stages) covers every stage the CLI accepts."""
    for stage in PROCESS_ALL_FORCE_STAGES:
        parcellation = "gif" if stage == "atlas" else "gif,freesurfer"
        rules = process_all_force_rules(stage, parcellation)
        assert {stage_of_rule(r) for r in rules} == {stage}


def test_message_names_rules_and_the_override():
    message = describe_unexpected({"tracula_trac_all_gif": 55}, ["metrics"])
    assert "tracula_trac_all_gif: 55 job(s)" in message
    assert "--allow-upstream-rerun" in message
    assert "Nothing was submitted" in message


class TestOnlyForced:
    """--only-forced: Snakemake may use only the forced rules (and the aggregation rule)."""

    def test_allowed_rules_are_the_forced_rules_plus_all(self):
        from leukoquant.utils.force_guard import only_forced_rules
        rules = ["metrics_extract_metrics_gif", "metrics_extract_metrics_freesurfer"]
        assert only_forced_rules(rules) == frozenset(rules) | {"all"}

    def test_upstream_rules_are_not_allowed(self):
        from leukoquant.utils.force_guard import only_forced_rules
        allowed = only_forced_rules(["metrics_extract_metrics_gif"])
        for upstream in ("tracula_trac_all_gif", "dti_dti", "noddi_noddi", "z_score_z_score", "tract_qc_tract_qc_gif"):
            assert upstream not in allowed

    def test_cli_exposes_the_option(self):
        import importlib
        from click.testing import CliRunner
        cli_main = importlib.import_module("leukoquant.cli.main")
        result = CliRunner().invoke(cli_main.main, ["process-all", "--help"])
        assert "--only-forced" in result.output and "--dry-run" in result.output
