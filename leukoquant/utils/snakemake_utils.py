"""Shared helpers for building Snakemake CLI invocations across processors."""
import os
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import yaml

# Per-run Snakemake config files live here, inside each output directory.
RUN_CONFIG_DIR = ".leukoquant_runs"


def write_run_config(record_path, config: Dict[str, object]) -> str:
    """Write this run's Snakemake config and return the path to hand Snakemake.

    Every job Snakemake submits re-reads the config file by path when it
    STARTS (its command carries --configfiles <path>), possibly hours or days
    after submission. Writing each run's config to the same fixed file in the
    output directory therefore let a later run silently change the config of
    an earlier run's still-queued jobs -- including their subject list, so
    those jobs no longer matched their own targets and failed at start. This
    happened on EPAD 2026-09-27: a 94-session rerun rewrote
    process_all_config.yaml and every queued job of the main run failed.

    So the config Snakemake reads is a per-run file,
    <dir>/.leukoquant_runs/<name>_<timestamp>_<random>.yaml, that no later
    run touches. `record_path` (e.g. <dir>/process_all_config.yaml) is still
    written too, as the record of the latest run's settings that people and
    analysis code read. Nothing a running pipeline depends on reads it.
    """
    record_path = Path(record_path)
    run_dir = record_path.parent / RUN_CONFIG_DIR
    run_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    run_path = run_dir / f"{record_path.stem}_{stamp}_{uuid.uuid4().hex[:8]}.yaml"
    with open(run_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False)
    # Record copy: temp file + rename, so a reader never sees a partial file.
    tmp = record_path.with_name(f".{record_path.name}.tmp{os.getpid()}")
    with open(tmp, "w") as f:
        yaml.dump(config, f, default_flow_style=False)
    os.replace(tmp, record_path)
    return str(run_path)


def load_yaml_config(config_path: Optional[str]) -> Dict[str, object]:
    """Load an optional processor config YAML, explicitly passed via --config-yaml.

    Returns {} when config_path is None -- callers merge these values in as
    defaults, with CLI-supplied arguments always taking priority. Raises if
    a path was explicitly given but doesn't exist, since a silently-ignored
    typo would be worse than a clear error.
    """
    if not config_path:
        return {}

    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"Config YAML not found: {config_path}")

    with config_file.open("r") as f:
        data = yaml.safe_load(f) or {}

    if not isinstance(data, dict):
        raise ValueError("Config YAML must contain a top-level mapping")

    return data


def first_truthy(*values: object) -> Optional[str]:
    """Return the first truthy value cast to str, or None.

    Used to let a CLI argument override the equivalent config-YAML value:
    call as first_truthy(cli_value, yaml_cfg.get("key")).
    """
    for v in values:
        if v is not None and v != {} and v != "":
            return str(v)
    return None


def add_rerun_triggers_args(snakemake_cmd: List[str]) -> None:
    """Append --rerun-triggers mtime, in place.

    Snakemake's default rerun_triggers includes CODE (plus PARAMS, INPUT,
    SOFTWARE_ENV) on top of the always-on mtime staleness check -- meaning
    any edit to a workflow file's shell/params, even an unrelated one-line
    bugfix in a rule that has nothing to do with a given subject's actual
    data, makes every subject with recorded metadata for that rule look
    "changed" and get rescheduled, whether real work is needed or not.
    Confirmed as the actual cause of two real incidents in the z-score
    pipeline (2026-08-24/25 and 2026-08-29): a routine same-day edit to
    z_score_workflow.smk correctly-per-Snakemake's-own-logic but
    incorrectly-for-our-purposes re-flagged ~250 already-complete subjects
    as needing a rerun (confirmed via Snakemake's own
    summary(detailed=True): "rule implementation changed" against subjects
    that already had complete, valid output).

    Restricting to mtime alone keeps genuine staleness detection (an input
    file actually changed) working exactly as before -- that check isn't
    gated by rerun_triggers at all, only CODE/PARAMS/INPUT/SOFTWARE_ENV are
    -- while removing the "we edited the pipeline today" false-positive
    path entirely. A genuine methodology change now requires an explicit
    --forcerun, same tradeoff already accepted for the z-score ancient()
    fix from 65054fe1. Applied to every processor, not just process-all,
    since the underlying risk (a routine edit to any rule's workflow file)
    isn't specific to any one tool.
    """
    # "=" form on purpose: --rerun-triggers takes one or more values, so the
    # space-separated form also swallows a following positional target such as
    # "all" ("invalid choice: 'all'"), which broke every standalone processor.
    snakemake_cmd.append("--rerun-triggers=mtime")


def add_forcerun_args(snakemake_cmd: List[str], force_rules: Optional[List[str]]) -> None:
    """Append --forcerun <rule> for each rule in force_rules, in place.

    Snakemake scopes --forcerun to the named rule(s) plus anything downstream
    of them in the DAG -- rules whose outputs are already up to date and not
    downstream of a forced rule are left untouched. This lets e.g. process-all
    be re-invoked to only recompute metrics for subjects whose outputs already
    exist, without redoing completed upstream stages (recon-all, DTI, GIF,
    tractography, ...).
    """
    if not force_rules:
        return
    for rule in force_rules:
        snakemake_cmd.extend(["--forcerun", rule])
