"""Unit tests for bash_utils.sh functions.

Tests make_tmp_path, get_job_id, run_with_peak_mem, scratch_cleanup,
and legacy no-ops in isolated bash subprocesses.
"""
import os
import subprocess
from pathlib import Path


BASH_UTILS = Path(__file__).resolve().parents[1] / "utils" / "bash_utils.sh"


def run_bash_func(snippet: str, env: dict = None) -> subprocess.CompletedProcess:
    """Run a bash snippet sourcing bash_utils.sh."""
    full_env = os.environ.copy()
    if env:
        full_env.update(env)
    cmd = f"source '{BASH_UTILS}'\n{snippet}"
    return subprocess.run(
        ["bash", "-c", cmd],
        capture_output=True,
        text=True,
        env=full_env,
        check=False,
    )


def test_make_tmp_path():
    """Test make_tmp_path with various file extensions."""
    res = run_bash_func("""
        make_tmp_path "/path/to/image.nii.gz"
        make_tmp_path "data.csv"
        make_tmp_path "/dir/archive.tar.bz2"
        make_tmp_path "/dir/noext"
    """)
    assert res.returncode == 0
    lines = res.stdout.strip().splitlines()
    assert lines[0] == "/path/to/image_tmp.nii.gz"
    assert lines[1] == "data_tmp.csv"
    assert lines[2] == "/dir/archive_tmp.tar.bz2"
    assert lines[3] == "/dir/noext_tmp"


def test_get_job_id_fallbacks():
    """Test get_job_id with SGE, SLURM, and default local fallback."""
    # Local fallback: blank every scheduler variable, because the test may itself
    # run inside an SGE/SLURM job whose own job ID would otherwise be picked up.
    no_scheduler = {v: "" for v in (
        "AWS_BATCH_JOB_ID", "LSB_JOBID", "SLURM_JOB_ID", "PBS_JOBID", "JOB_ID", "SGE_TASK_ID")}
    res = run_bash_func("get_job_id", env=no_scheduler)
    assert res.returncode == 0
    assert res.stdout.strip().startswith("local_")

    # SGE job id
    res_sge = run_bash_func("get_job_id", env={"JOB_ID": "12345", "SGE_TASK_ID": "2"})
    assert res_sge.returncode == 0
    assert res_sge.stdout.strip() == "12345_2"

    # SLURM job id
    res_slurm = run_bash_func("get_job_id", env={"SLURM_JOB_ID": "98765", "SLURM_ARRAY_TASK_ID": "4"})
    assert res_slurm.returncode == 0
    assert res_slurm.stdout.strip() == "98765_4"


def test_scratch_cleanup(tmp_path):
    """Test scratch_cleanup reports scratch size and removes directory."""
    scratch_dir = tmp_path / "scratch_work"
    scratch_dir.mkdir()
    dummy_file = scratch_dir / "large.bin"
    dummy_file.write_bytes(b"\x00" * (1024 * 1024 * 5))  # 5 MB file

    res = run_bash_func(f"""
        scratch_cleanup "{scratch_dir}"
    """)
    assert res.returncode == 0
    assert "Cleaning up scratch directory:" in res.stdout
    assert "Scratch usage before cleanup:" in res.stdout
    assert "Scratch cleanup complete." in res.stdout
    assert not scratch_dir.exists()


def test_scratch_cleanup_empty_or_missing():
    """Test scratch_cleanup handles empty string or non-existent path cleanly."""
    res = run_bash_func("""
        scratch_cleanup ""
    """)
    assert res.returncode == 0
    assert "Scratch cleanup complete." not in res.stdout


def test_run_with_peak_mem():
    """Test run_with_peak_mem executes command and logs peak RSS."""
    res = run_bash_func("""
        run_with_peak_mem "test_stage" python3 -c "x = bytearray(20000000); x[0] = 1"
    """)
    assert res.returncode == 0
    combined_output = res.stdout + res.stderr
    assert "[mem] Peak RSS (test_stage):" in combined_output


def test_run_with_peak_mem_propagates_exit_code():
    """Test run_with_peak_mem preserves non-zero exit codes from failed commands."""
    res = run_bash_func("""
        run_with_peak_mem "failing_cmd" bash -c "exit 42"
    """)
    assert res.returncode == 42


def test_run_with_peak_mem_under_set_e():
    """Test run_with_peak_mem under set -e correctly halts on error and preserves exit code."""
    res = run_bash_func("""
        set -e
        run_with_peak_mem "ok_step" echo "step 1"
        run_with_peak_mem "bad_step" bash -c "exit 33"
        echo "unreachable"
    """)
    assert res.returncode == 33
    assert "step 1" in res.stdout
    assert "unreachable" not in res.stdout
    assert "[mem] Peak RSS (ok_step):" in (res.stdout + res.stderr)
    assert "[mem] Peak RSS (bad_step):" in (res.stdout + res.stderr)


def test_run_with_peak_mem_argument_preservation():
    """Test run_with_peak_mem preserves spaces and complex arguments."""
    res = run_bash_func("""
        run_with_peak_mem "args_step" python3 -c 'import sys; print(" ".join(sys.argv[1:]))' "arg with spaces" --opt=val
    """)
    assert res.returncode == 0
    assert "arg with spaces --opt=val" in res.stdout


def test_multi_stage_peak_mem_summary():
    """Test multi-stage commands track each stage and report overall total at the end."""
    res = run_bash_func("""
        run_with_peak_mem "stage_one" python3 -c "x = bytearray(20000000); x[0] = 1"
        run_with_peak_mem "stage_two" python3 -c "x = bytearray(50000000); x[0] = 1"
        run_with_peak_mem "stage_three" python3 -c "x = bytearray(10000000); x[0] = 1"
        report_total_peak_mem
    """)
    assert res.returncode == 0
    combined = res.stdout + res.stderr
    assert "[mem] Peak RSS (stage_one):" in combined
    assert "[mem] Peak RSS (stage_two):" in combined
    assert "[mem] Peak RSS (stage_three):" in combined
    assert "[mem] --- Peak RSS Summary (3 stages) ---" in combined
    assert "* stage_one:" in combined
    assert "* stage_two:" in combined
    assert "* stage_three:" in combined
    assert "[mem] Total Peak RSS:" in combined


def test_single_stage_total_peak_mem():
    """Test single stage command reports total peak RSS without multi-stage list header."""
    res = run_bash_func("""
        run_with_peak_mem "single_stage" python3 -c "x = bytearray(25000000); x[0] = 1"
        report_total_peak_mem
    """)
    assert res.returncode == 0
    combined = res.stdout + res.stderr
    assert "[mem] Peak RSS (single_stage):" in combined
    assert "[mem] --- Peak RSS Summary" not in combined
    assert "[mem] Total Peak RSS:" in combined


def test_report_total_peak_mem_idempotent():
    """Test report_total_peak_mem only reports once when called multiple times."""
    res = run_bash_func("""
        run_with_peak_mem "test_step" python3 -c "x = bytearray(15000000); x[0] = 1"
        report_total_peak_mem
        report_total_peak_mem
    """)
    assert res.returncode == 0
    # In stdout alone, Total Peak RSS should appear exactly once
    assert res.stdout.count("[mem] Total Peak RSS:") == 1


def test_scratch_cleanup_invokes_report_total_peak_mem(tmp_path):
    """Test scratch_cleanup reports total peak RSS before deleting scratch."""
    scratch_dir = tmp_path / "scratch_test"
    scratch_dir.mkdir()
    res = run_bash_func(f"""
        run_with_peak_mem "prep_step" python3 -c "x = bytearray(15000000); x[0] = 1"
        scratch_cleanup "{scratch_dir}"
    """)
    assert res.returncode == 0
    combined = res.stdout + res.stderr
    assert "[mem] Total Peak RSS:" in combined
    assert "Scratch cleanup complete." in res.stdout
    assert not scratch_dir.exists()


def test_legacy_monitor_noops():
    """Test start_peak_mem_monitor and stop_peak_mem_monitor are harmless no-ops."""
    res = run_bash_func("""
        start_peak_mem_monitor
        stop_peak_mem_monitor
    """)
    assert res.returncode == 0
