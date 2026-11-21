"""The skull-strip masks must not depend on the node's CPU.

Background (ADNI-3, 2026-11-21): mri_synthstrip (a PyTorch network on the CPU) gave a different T1 mask on an AVX2 (Haswell)
node than on AVX-512 (Xeon Gold) and AMD EPYC nodes for the same input (raw mask md5 2d1078ed.. against ee67f545..). The
mask difference moved the DWI-to-T1 affine and with it the DTI/NODDI values sampled in the lesions. The cause is MKL, which
selects its kernels from the CPU at run time. MKL_CBWR=COMPATIBLE restricts it to kernels every x86 CPU produces the same
bits with. With it, the Haswell, Xeon Gold and AMD nodes gave identical b0 and T1 masks and identical reg_aladin affines
for three sessions, two repeats each. ONEDNN_MAX_CPU_ISA, one thread and the OpenBLAS pin did not change the mask.
The variable is read when MKL is loaded, so it has to be exported before mri_synthstrip starts.
"""
import os
import re
import subprocess
from pathlib import Path

UTILS = Path(__file__).resolve().parents[1] / "utils"
BASH_UTILS = UTILS / "bash_utils.sh"
METRICS_SH = UTILS / "metrics_calc.sh"


def _bash(snippet, env):
    return subprocess.run(["bash", "-c", f"source '{BASH_UTILS}'\n{snippet}"], capture_output=True, text=True, env=env, check=False)


def _without_pin():
    env = os.environ.copy()
    env.pop("MKL_CBWR", None)
    return env


def test_pin_mkl_exports_the_compatible_mode():
    res = _bash('pin_mkl_kernel; echo "mode=$MKL_CBWR"', _without_pin())
    assert res.returncode == 0
    assert "mode=COMPATIBLE" in res.stdout


def test_pin_mkl_logs_what_it_did():
    res = _bash("pin_mkl_kernel", _without_pin())
    assert re.fullmatch(r"\[mkl\] MKL_CBWR=COMPATIBLE \(.+\)", res.stdout.strip()), res.stdout


def test_pin_mkl_respects_a_value_set_by_the_caller():
    env = _without_pin()
    env["MKL_CBWR"] = "AUTO"
    res = _bash('pin_mkl_kernel; echo "mode=$MKL_CBWR"', env)
    assert "mode=AUTO" in res.stdout and "already set" in res.stdout


def test_pin_mkl_never_fails_the_job():
    res = _bash("set -e; pin_mkl_kernel; echo after", _without_pin())
    assert res.returncode == 0 and res.stdout.strip().endswith("after")


def test_metrics_script_pins_mkl_before_any_skull_strip_or_python():
    lines = METRICS_SH.read_text().splitlines()
    pin = next(i for i, l in enumerate(lines) if re.fullmatch(r"pin_mkl_kernel", l))
    first_use = next(i for i, l in enumerate(lines)
                     if not l.lstrip().startswith("#") and (re.search(r"\brun_synthstrip\b|\bmri_synthstrip\b", l)
                                                            or re.search(r'(^|[\s"(])(python3|"?\$PYTHON3"?)\s', l) and "PYTHON3=" not in l))
    assert pin < first_use


def test_the_dwi_to_t1_registration_runs_after_the_pin():
    # register_dwi_to_t1.sh skull-strips the b0 with mri_synthstrip; it is called from the metrics script, so the exported
    # variable reaches it, provided the call comes after the pin.
    text = METRICS_SH.read_text()
    assert "register_dwi_to_t1.sh" in text
    assert text.index("pin_mkl_kernel") < text.index("register_dwi_to_t1.sh")
