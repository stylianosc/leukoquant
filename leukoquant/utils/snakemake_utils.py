"""Shared helpers for building Snakemake CLI invocations across processors."""
from pathlib import Path
from typing import Dict, List, Optional

import yaml


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
