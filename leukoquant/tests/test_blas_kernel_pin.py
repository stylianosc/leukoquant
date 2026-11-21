"""The metrics step must give the same numbers on every cluster node.

Background (EPAD, 2026-11-18): the same metrics job with the same code and inputs gave slightly different tract masks
on the abner nodes (Xeon E5-2650 v2, AVX only) than on the other nodes (AVX2/FMA, some AVX-512). Every checksum before
the Python streamline step (registration affines, control points, deformation fields, skull-stripped images, DTI and
NODDI maps in T1 space) was identical; the first differences were the tractography density maps written by
transform_trk.py. The cause is OpenBLAS: numpy's bundled OpenBLAS picks its matrix kernel at run time from the CPU, a
kernel with FMA rounds differently from one without, and the last bits of the transformed streamline coordinates then
decide which voxel a point near a voxel boundary falls in. Setting OPENBLAS_CORETYPE to one kernel that every node
supports made the density maps of an abner node, an AVX2 node and an AVX-512 node byte-identical.

So the metrics step pins the kernel before it starts any Python.
"""
import os
import re
import subprocess
from pathlib import Path

import pytest

UTILS = Path(__file__).resolve().parents[1] / "utils"
BASH_UTILS = UTILS / "bash_utils.sh"
METRICS_SH = UTILS / "metrics_calc.sh"


def _bash(snippet, env=None):
    return subprocess.run(["bash", "-c", f"source '{BASH_UTILS}'\n{snippet}"], capture_output=True, text=True,
                          env=env or os.environ.copy(), check=False)


def _without_pin():
    env = os.environ.copy()
    env.pop("OPENBLAS_CORETYPE", None)
    return env


def test_pin_blas_kernel_exports_a_kernel_every_node_supports():
    res = _bash('pin_blas_kernel; echo "kernel=$OPENBLAS_CORETYPE"', env=_without_pin())
    assert res.returncode == 0
    # Core2 needs only SSSE3, so the oldest node we know (arbuckle, SSE4.2, no AVX) can run it; SandyBridge killed it with SIGILL.
    assert "kernel=Core2" in res.stdout


def test_pin_blas_kernel_logs_what_it_did():
    res = _bash("pin_blas_kernel", env=_without_pin())
    assert re.fullmatch(r"\[blas\] OPENBLAS_CORETYPE=Core2 \(.+\)", res.stdout.strip()), res.stdout


def test_pin_blas_kernel_respects_a_kernel_set_by_the_caller():
    # Lets a developer test another kernel; the log line then says so.
    env = _without_pin()
    env["OPENBLAS_CORETYPE"] = "Haswell"
    res = _bash('pin_blas_kernel; echo "kernel=$OPENBLAS_CORETYPE"', env=env)
    assert "kernel=Haswell" in res.stdout
    assert "already set" in res.stdout


def test_pin_blas_kernel_never_fails_the_job():
    res = _bash("set -e; pin_blas_kernel; echo after", env=_without_pin())
    assert res.returncode == 0 and res.stdout.strip().endswith("after")


def test_metrics_script_pins_the_kernel_before_any_python():
    lines = METRICS_SH.read_text().splitlines()
    # The pin is a top-level call (exported to every child process), not inside a function or branch.
    pin = next(i for i, l in enumerate(lines) if re.fullmatch(r"pin_blas_kernel", l))
    # First line of code (comments skipped) that starts a Python interpreter, by name or through $PYTHON3.
    first_python = next(i for i, l in enumerate(lines)
                        if not l.lstrip().startswith("#") and re.search(r'(^|[\s"(])(python3|"?\$PYTHON3"?)\s', l)
                        and "PYTHON3=" not in l)
    assert pin < first_python


def test_metrics_script_pins_before_the_provenance_lines_are_written():
    text = METRICS_SH.read_text()
    assert text.index("pin_blas_kernel") < text.index("log_host_cpu")


def _numpy_openblas_architecture():
    """Kernel numpy's OpenBLAS selected in a fresh process started with the pin, or None if not applicable."""
    code = ("from threadpoolctl import threadpool_info as t\n"
            "print(','.join(sorted({str(d.get('architecture')) for d in t() if d.get('internal_api') == 'openblas'})))")
    env = os.environ.copy()
    env["OPENBLAS_CORETYPE"] = "Core2"
    res = subprocess.run(["python3", "-c", code], capture_output=True, text=True, env=env, check=False)
    return res.stdout.strip() if res.returncode == 0 else None


def test_numpy_openblas_honours_the_pin():
    arch = _numpy_openblas_architecture()
    if not arch:
        pytest.skip("threadpoolctl or an OpenBLAS-backed numpy is not available here")
    import platform
    if platform.machine().lower() not in ("x86_64", "amd64"):
        pytest.skip("OPENBLAS_CORETYPE selects x86 kernels only")
    assert arch.lower() == "core2", arch
