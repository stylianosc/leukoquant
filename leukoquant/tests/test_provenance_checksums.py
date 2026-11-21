"""Provenance lines in the metrics log: the node's CPU and a checksum of every intermediate.

Background (EPAD, 2026-11-18): the same metrics job with the same code and inputs gave tract volumes 0.3-0.7% apart in
about one session in five, apparently depending on the node it ran on, while the logs were identical. The log must
show where two runs part ways: which node and CPU ran the job, and a checksum of each transform and each T1-space
input handed to metrics_calc.py, so two runs can be diffed line by line.
"""
import os
import re
import subprocess
from pathlib import Path

UTILS = Path(__file__).resolve().parents[1] / "utils"
BASH_UTILS = UTILS / "bash_utils.sh"
METRICS_SH = UTILS / "metrics_calc.sh"


def _bash(snippet, cwd=None):
    return subprocess.run(["bash", "-c", f"source '{BASH_UTILS}'\n{snippet}"], capture_output=True, text=True,
                          cwd=cwd, env=os.environ.copy(), check=False)


def test_log_host_cpu_names_host_and_cpu():
    res = _bash("log_host_cpu")
    assert res.returncode == 0
    line = res.stdout.strip()
    assert line.startswith("[host] "), line
    assert "cpu=" in line and "simd=" in line


def test_log_checksums_lists_files_sorted_with_md5(tmp_path):
    (tmp_path / "d" / "sub").mkdir(parents=True)
    (tmp_path / "d" / "b.txt").write_text("b")
    (tmp_path / "d" / "sub" / "a.txt").write_text("a")
    (tmp_path / "aff.txt").write_text("1 0 0 0\n")
    res = _bash(f"log_checksums affine '{tmp_path}/aff.txt'; log_checksums inputs '{tmp_path}/d'")
    assert res.returncode == 0
    lines = res.stdout.strip().splitlines()
    assert re.fullmatch(r"\[sha\] affine aff\.txt [0-9a-f]{32}", lines[0]), lines[0]
    # Directory entries are relative to the directory and sorted, so two runs diff cleanly.
    assert [l.split()[2] for l in lines[1:]] == ["b.txt", "sub/a.txt"]
    assert lines[1].split()[3] == "92eb5ffee6ae2fec3ad71c777531578f"  # md5 of "b"


def test_log_checksums_reports_missing_and_never_fails(tmp_path):
    res = _bash(f"set -e; log_checksums cpp '{tmp_path}/nope.nii.gz'; echo after")
    assert res.returncode == 0
    assert "[sha] cpp nope.nii.gz missing" in res.stdout
    assert "[sha] atlas (unset)" in _bash("log_checksums atlas ''").stdout
    assert res.stdout.strip().endswith("after")


def test_metrics_script_logs_provenance_before_the_python_step():
    text = METRICS_SH.read_text()
    py_call = text.index('PY_METRICS_PREP="$SCRIPT_DIR/metrics_calc.py"')
    assert "log_host_cpu" in text[:py_call]
    # The DWI->T1 registration folder (affines, skull-stripped b0), the atlas registration folder (affine, control
    # points, deformation) and every T1-space input of the Python step.
    for target in ('"$REG_DIR"', '"$FA_REGISTERED_DIR"', '"$T1_SPACE_DIR"'):
        assert re.search(r"log_checksums \S+ " + re.escape(target), text[:py_call]), target
