"""Out-of-memory and other step failures must be loud: printed lines and a non-zero status.

Background: mri_synthstrip ran out of memory under the SGE h_vmem limit in some sessions. Its output went to
/dev/null (metrics_calc.sh) or the failure was only seen in the error log, the rule printed a normal
'[mem] Peak RSS' line, and the stale or empty outputs looked like a result. These tests simulate the failing
tool and a SIGKILLed stage and check what the census can rely on:

  status       the stage's exit status is passed on (non-zero)
  log lines    'OUT OF MEMORY: ...', 'FAILED: ...', '[oom] <stage> ...', '(stage FAILED, exit N)'
"""
import os
import stat
import subprocess
from pathlib import Path

import pytest

UTILS = Path(__file__).resolve().parents[1] / "utils"
BASH_UTILS = UTILS / "bash_utils.sh"

# GNU time (Linux, the cluster containers) reports a SIGKILLed command as status 137; BSD time on macOS does not.
_GNU_TIME = subprocess.run(["/usr/bin/time", "-f", "%M", "true"], capture_output=True).returncode == 0 \
    if os.path.exists("/usr/bin/time") else False


def run(snippet, tmp_path, extra_env=None, path_prefix=None):
    """Run a bash snippet with bash_utils.sh sourced; stderr goes to tmp_path/metrics_error.txt like in the rule."""
    err = tmp_path / "metrics_error.txt"
    err.touch()
    env = os.environ.copy()
    env["LEUKOQUANT_STDERR_FILE"] = str(err)
    if path_prefix:
        env["PATH"] = f"{path_prefix}{os.pathsep}{env['PATH']}"
    if extra_env:
        env.update(extra_env)
    script = f"source '{BASH_UTILS}'\n{snippet}\n"
    with open(err, "a") as err_fh:
        res = subprocess.run(["bash", "-c", script, "bash"], stdout=subprocess.PIPE, stderr=err_fh,
                             text=True, env=env, check=False)
    return res, err.read_text()


def fake_tool(tmp_path, name, body):
    """A stand-in executable placed first on PATH."""
    d = tmp_path / "bin"
    d.mkdir(exist_ok=True)
    f = d / name
    f.write_text("#!/bin/bash\n" + body + "\n")
    f.chmod(f.stat().st_mode | stat.S_IEXEC)
    return str(d)


ALLOC_MSG = "mri_synthstrip: can't allocate memory: 5033164800 bytes"


def test_stage_out_of_memory_message_is_reported(tmp_path):
    """A stage whose stderr says 'can't allocate memory' and exits 1: status kept, loud lines."""
    res, err = run(f"""
        run_with_peak_mem "metrics" bash -c "echo \\"{ALLOC_MSG}\\" >&2; exit 1"
    """, tmp_path)
    assert res.returncode == 1
    out = res.stdout + err
    assert "(stage FAILED, exit 1)" in out
    assert "[mem] Peak RSS (metrics):" in out          # prefix unchanged: existing parsers still match
    assert "[oom] metrics: out of memory" in out
    assert "can't allocate memory" in out
    assert "marker" not in out.lower()
    assert sorted(p.name for p in tmp_path.iterdir()) == ["metrics_error.txt"]   # no extra files are written


@pytest.mark.skipif(not _GNU_TIME, reason="needs GNU time, which reports a killed command as status 137")
def test_stage_killed_by_sigkill_is_reported_as_out_of_memory(tmp_path):
    """Exit status 137 (SIGKILL: memory limit, time limit or the OOM killer) is reported as out of memory."""
    res, err = run("""
        run_with_peak_mem "metrics" bash -c 'kill -9 $$'
    """, tmp_path)
    assert res.returncode == 137
    out = res.stdout + err
    assert "[oom] metrics: out of memory (exit 137" in out
    assert "(stage FAILED, exit 137)" in out


def test_stage_failing_for_another_reason_is_not_called_oom(tmp_path):
    """A plain failure is marked as failed on the Peak line but not as out of memory."""
    res, err = run("""
        run_with_peak_mem "metrics" bash -c 'echo "bad input" >&2; exit 3'
    """, tmp_path)
    assert res.returncode == 3
    out = res.stdout + err
    assert "(stage FAILED, exit 3)" in out
    assert "[oom]" not in out


def test_successful_stage_prints_no_failure_text(tmp_path):
    """A stage that succeeds keeps the plain Peak line and prints no failure or [oom] text."""
    res, err = run("""
        run_with_peak_mem "metrics" bash -c 'echo fine'
    """, tmp_path)
    assert res.returncode == 0
    out = res.stdout + err
    assert "[mem] Peak RSS (metrics):" in out
    assert "(stage FAILED" not in out and "[oom]" not in out and "OUT OF MEMORY" not in out


def test_swallowed_out_of_memory_message_is_flagged_even_with_status_zero(tmp_path):
    """A stage that hides a failure (exit 0) but printed an allocation error is flagged, status unchanged."""
    res, err = run(f"""
        run_with_peak_mem "metrics" bash -c "echo \\"{ALLOC_MSG}\\" >&2; exit 0"
    """, tmp_path)
    assert res.returncode == 0
    assert "[oom] metrics: exit 0 but the log contains an out-of-memory message" in res.stdout + err


def test_earlier_text_in_the_error_log_is_not_blamed_on_a_later_stage(tmp_path):
    """Only the text a stage itself wrote to stderr is scanned."""
    res, err = run(f"""
        echo "{ALLOC_MSG}" >&2
        run_with_peak_mem "metrics" bash -c 'echo fine'
    """, tmp_path)
    assert res.returncode == 0
    assert "[oom]" not in res.stdout + err


def test_run_synthstrip_reports_out_of_memory(tmp_path):
    """A skull strip that cannot allocate: tool output, legacy ERROR line, OUT OF MEMORY line, status 1."""
    path = fake_tool(tmp_path, "mri_synthstrip", f'echo "{ALLOC_MSG}" >&2; exit 1')
    res, err = run('run_synthstrip "T1" in.nii.gz out.nii.gz', tmp_path, path_prefix=path)
    assert res.returncode == 1
    out = res.stdout + err
    assert "ERROR: mri_synthstrip failed on T1. Output:" in err
    assert "can't allocate memory: 5033164800 bytes" in err
    assert "OUT OF MEMORY: mri_synthstrip (T1) ran out of memory (exit 1" in out


def test_run_synthstrip_other_failure_is_labelled_failed(tmp_path):
    path = fake_tool(tmp_path, "mri_synthstrip", 'echo "unreadable image" >&2; exit 2')
    res, err = run('run_synthstrip "b0" in.nii.gz out.nii.gz', tmp_path, path_prefix=path)
    assert res.returncode == 2
    assert "FAILED: mri_synthstrip (b0) (exit 2)" in res.stdout + err
    assert "OUT OF MEMORY" not in res.stdout + err


def test_run_synthstrip_success_is_silent(tmp_path):
    path = fake_tool(tmp_path, "mri_synthstrip", "exit 0")
    res, err = run('run_synthstrip "T1" in.nii.gz out.nii.gz', tmp_path, path_prefix=path)
    assert res.returncode == 0
    assert "OUT OF MEMORY" not in res.stdout + err and "FAILED" not in res.stdout + err


def test_status_137_with_no_text_is_reported_as_out_of_memory(tmp_path):
    """The two ADNI-3 sessions that died without a message: a SIGKILL status alone is enough."""
    res, err = run('report_step_failure "mri_synthstrip (T1 (TIV))" 137 "" || true', tmp_path)
    assert "OUT OF MEMORY: mri_synthstrip (T1 (TIV)) ran out of memory (exit 137" in res.stdout + err


def test_step_script_with_failing_synthstrip_stops_loudly(tmp_path):
    """A miniature of the metrics rule: set -eE, ERR trap, skull strip, then a step that must not run.

    The stage fails, says why, and the step after the skull strip never runs."""
    path = fake_tool(tmp_path, "mri_synthstrip", f'echo "{ALLOC_MSG}" >&2; exit 1')
    step = tmp_path / "step.sh"
    step.write_text(f"""#!/bin/bash
set -eE
source '{BASH_UTILS}'
trap 'report_failed_command "$?" "$LINENO" "$BASH_COMMAND"' ERR
run_synthstrip "T1 (TIV)" in.nii.gz out.nii.gz || exit 1
echo "REACHED THE STEP AFTER THE SKULL STRIP"
""")
    res, err = run(f"run_with_peak_mem metrics bash '{step}'", tmp_path, path_prefix=path)
    assert res.returncode == 1
    out = res.stdout + err
    assert "REACHED THE STEP" not in out
    assert "OUT OF MEMORY: mri_synthstrip (T1 (TIV))" in out
    assert "(stage FAILED, exit 1)" in out


def test_err_trap_names_a_silent_failure(tmp_path):
    """A command that fails under set -e with no output now leaves a line naming it."""
    step = tmp_path / "step.sh"
    step.write_text(f"""#!/bin/bash
set -eE
source '{BASH_UTILS}'
trap 'report_failed_command "$?" "$LINENO" "$BASH_COMMAND"' ERR
false
echo unreachable
""")
    res, err = run(f"bash '{step}'", tmp_path)
    assert res.returncode == 1
    assert "FAILED: step.sh line 5: false (exit 1)" in res.stdout + err


def test_pipeline_scripts_use_the_loud_helpers():
    """The skull-strip call sites no longer discard their output, and the scripts still parse."""
    metrics = (UTILS / "metrics_calc.sh").read_text()
    register = (UTILS / "register_dwi_to_t1.sh").read_text()
    assert "mri_synthstrip -i" not in metrics and "run_synthstrip" in metrics
    assert "mri_synthstrip -i" not in register and register.count("run_synthstrip") >= 2
    assert "trap cleanup EXIT INT TERM" in metrics and "report_failed_command" in metrics
    for script in ("metrics_calc.sh", "register_dwi_to_t1.sh", "bash_utils.sh"):
        assert subprocess.run(["bash", "-n", str(UTILS / script)]).returncode == 0
