"""--rerun-triggers must not swallow the positional Snakemake target."""
import argparse

from leukoquant.utils.snakemake_utils import add_rerun_triggers_args


def test_target_after_rerun_triggers_is_not_taken_as_a_trigger():
    cmd = ["snakemake", "-s", "wf.smk"]
    add_rerun_triggers_args(cmd)
    cmd.append("all")

    # Mirror snakemake's option: nargs="+" with a restricted set of choices.
    parser = argparse.ArgumentParser()
    parser.add_argument("-s")
    parser.add_argument("--rerun-triggers", nargs="+", choices=["code", "input", "mtime", "params", "software-env"])
    parser.add_argument("targets", nargs="*")
    args = parser.parse_args(cmd[1:])

    assert args.rerun_triggers == ["mtime"]
    assert args.targets == ["all"]
