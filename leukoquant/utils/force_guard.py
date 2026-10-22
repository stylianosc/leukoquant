"""Safety check for ``process-all --force``.

A forced run is meant to recompute the forced stages (and whatever depends on them). If Snakemake also
schedules a stage *upstream* of them (TRACULA, DTI, ...), it deletes that stage's existing outputs before
rerunning it. That happened on 2026-10-05 (DPUK) and 2026-10-18 (ADNI-3: TRACULA ``dmri`` outputs of 55
sessions), so a forced run is first planned as a dry run, and refused when the plan contains such a stage.

This module only reads a dry-run plan and decides; it does not run Snakemake.
"""
import re
from typing import Dict, Iterable, List

# Every stage of process-all in dependency order (earliest first). A stage listed after a forced stage depends
# on it and may legitimately rerun. Keep in step with ``PROCESS_ALL_FORCE_STAGES`` (leukoquant/cli/commands.py).
STAGE_ORDER = ["recon", "gif", "bamos", "atlas", "tracula", "dti", "noddi", "tractqc", "zscore", "metrics"]

# Rule-name prefix of each stage in the process-all DAG (sub-workflow modules are prefixed by module name).
# ``tracula_other_`` is the atlas conversion nested in the tracula module, so it is matched before ``tracula_``.
_STAGE_PREFIXES = [
    ("atlas", "tracula_other_convert_atlas"),
    ("recon", "recon_"),
    ("gif", "gif_"),
    ("bamos", "bamos_"),
    ("tracula", "tracula_"),
    ("dti", "dti_"),
    ("noddi", "noddi_"),
    ("tractqc", "tract_qc_"),
    ("zscore", "z_score_"),
    ("metrics", "metrics_"),
]

# Rules that are bookkeeping or aggregation, never a processing stage.
_NEUTRAL_RULES = {"all"}

_STATS_ROW = re.compile(r"^\s*(\S+)\s+(\d+)(?:\s+\S+)*\s*$")


def stage_of_rule(rule: str) -> str:
    """Name of the stage a Snakemake rule belongs to, or ``""`` when it is not a processing stage."""
    for stage, prefix in _STAGE_PREFIXES:
        if rule.startswith(prefix):
            return stage
    return ""


def parse_job_stats(output: str) -> Dict[str, int]:
    """Jobs per rule from the ``Job stats:`` table Snakemake prints for a dry run (``{}`` when there is none).

    Several tables can occur (one per planning pass); the counts of the last one are returned because it is
    the final plan.
    """
    counts: Dict[str, int] = {}
    in_table = False
    for line in output.splitlines():
        if line.strip().startswith("Job stats:"):
            counts, in_table = {}, True
            continue
        if not in_table:
            continue
        stripped = line.strip()
        if not stripped or stripped.startswith("job ") or set(stripped) <= {"-", " "}:
            continue
        match = _STATS_ROW.match(line)
        if match is None:
            in_table = False
            continue
        if match.group(1) != "total":
            counts[match.group(1)] = int(match.group(2))
    return counts


def unexpected_upstream_rules(job_counts: Dict[str, int], forced_stages: Iterable[str]) -> Dict[str, int]:
    """Planned rules that belong to a stage *upstream* of every forced stage.

    A stage is allowed when it is a forced stage or comes after the earliest forced stage in ``STAGE_ORDER``.
    Rules that match no stage are ignored (aggregation rules).
    """
    forced_positions = [STAGE_ORDER.index(s) for s in forced_stages if s in STAGE_ORDER]
    if not forced_positions:
        return {}
    earliest = min(forced_positions)
    unexpected: Dict[str, int] = {}
    for rule, count in job_counts.items():
        if count <= 0 or rule in _NEUTRAL_RULES:
            continue
        stage = stage_of_rule(rule)
        if stage and STAGE_ORDER.index(stage) < earliest:
            unexpected[rule] = count
    return unexpected


def describe_unexpected(unexpected: Dict[str, int], forced_stages: Iterable[str]) -> str:
    """Error message for a refused run."""
    rows: List[str] = [f"  {rule}: {count} job(s)" for rule, count in sorted(unexpected.items())]
    return (
        "process-all --force "
        + ",".join(forced_stages)
        + " would also rerun stages upstream of the forced one, which deletes their existing outputs:\n"
        + "\n".join(rows)
        + "\nNothing was submitted. Check why those outputs look out of date (a fresh checkout, a different "
        "working folder and moved inputs are typical causes). Pass --allow-upstream-rerun only when rerunning them is intended."
    )


def only_forced_rules(force_rules: Iterable[str]) -> frozenset:
    """Rules a run restricted to the forced stages (``--only-forced``) may use: the forced rules and the aggregation rule.

    Every other rule is left out of the plan, so a stage upstream of the forced one is never scheduled, whatever the
    modification times say. A file such an upstream stage should have produced must exist; if it does not, Snakemake stops
    with a missing-input error instead of rerunning (and deleting the outputs of) that stage."""
    return frozenset(force_rules) | frozenset(_NEUTRAL_RULES)
