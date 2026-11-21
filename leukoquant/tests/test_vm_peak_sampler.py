"""Per-stage peak address space (VmPeak) and resident (VmHWM) of the stage's processes, sampled from /proc.

Why: the SGE limit (h_vmem) acts on address space, but the pipeline only measured resident memory (GNU time %M), so
nothing showed how large the address space of a TRACULA stage gets. The sampler reads /proc/<pid>/status of the stage's
descendants every few seconds and reports the largest VmPeak and VmHWM at stage end. A process that starts and ends
between two samples is missed (stated in the log line). The proc directory is a variable so these tests use a fake one."""
import os
import subprocess
from pathlib import Path

import pytest

BASH_UTILS = Path(__file__).resolve().parents[1] / "utils" / "bash_utils.sh"
SMK = Path(__file__).resolve().parents[1] / "workflow" / "workflows" / "tracula_workflow.smk"


def make_proc(root, pid, ppid, name, vmpeak_kb, vmhwm_kb):
    d = root / str(pid)
    d.mkdir(parents=True)
    (d / "status").write_text(f"Name:\t{name}\nPid:\t{pid}\nPPid:\t{ppid}\nVmPeak:\t{vmpeak_kb} kB\nVmSize:\t1 kB\nVmHWM:\t{vmhwm_kb} kB\n")


def bash(snippet, proc_root, extra_env=None):
    env = os.environ.copy()
    env["LEUKOQUANT_PROC_DIR"] = str(proc_root)
    env.update(extra_env or {})
    return subprocess.run(["bash", "-c", f"source '{BASH_UTILS}'; {snippet}"], capture_output=True, text=True, env=env, timeout=60)


def test_scan_returns_the_largest_values_among_descendants_only(tmp_path):
    make_proc(tmp_path, 100, 1, "root_shell", 500, 400)         # the root itself is included
    make_proc(tmp_path, 101, 100, "trac-all", 5_000_000, 100_000)
    make_proc(tmp_path, 102, 101, "dmri_paths", 9_000_000, 3_000_000)   # grandchild
    make_proc(tmp_path, 200, 1, "other_users_job", 99_999_999, 99_999_999)   # not a descendant
    res = bash("_vm_scan 100", tmp_path)
    assert res.returncode == 0
    vmpeak_kb, vmpeak_name, vmhwm_kb, vmhwm_name, vmsum_kb = res.stdout.split()
    assert (vmpeak_kb, vmpeak_name) == ("9000000", "dmri_paths")
    assert (vmhwm_kb, vmhwm_name) == ("3000000", "dmri_paths")


def test_scan_excludes_the_given_pid_such_as_the_sampler_itself(tmp_path):
    make_proc(tmp_path, 100, 1, "root_shell", 500, 400)
    make_proc(tmp_path, 150, 100, "sampler", 777_000_000, 1)
    res = bash("_vm_scan 100 150", tmp_path)
    assert res.stdout.split()[0] == "500"


def test_scan_without_proc_prints_nothing_and_succeeds(tmp_path):
    res = bash("_vm_scan 100", tmp_path / "missing")
    assert res.returncode == 0 and res.stdout.strip() == ""


def test_sampler_reports_the_peak_seen_during_the_stage(tmp_path):
    make_proc(tmp_path, 100, 1, "root_shell", 500, 400)
    make_proc(tmp_path, 101, 100, "trac-all", 8_388_608, 2_097_152)     # 8192 MB VmPeak, 2048 MB VmHWM
    res = bash("start_vm_sampler 'trac-all path' 100; sleep 1; stop_vm_sampler", tmp_path,
               {"LEUKOQUANT_VM_SAMPLE_S": "0.2"})
    assert res.returncode == 0
    line = [l for l in res.stdout.splitlines() if l.startswith("[vm]")][0]
    assert "trac-all path" in line
    assert "VmPeak 8192 MB (trac-all)" in line
    assert "VmHWM 2048 MB (trac-all)" in line
    assert "address-space limit" in line


def test_sampler_says_so_when_proc_is_unavailable(tmp_path):
    res = bash("start_vm_sampler 'prep' 100; stop_vm_sampler", tmp_path / "missing")
    assert res.returncode == 0
    assert "[vm] prep: unavailable" in res.stdout


def test_sampler_leaves_no_background_process_behind(tmp_path):
    make_proc(tmp_path, 100, 1, "root_shell", 500, 400)
    res = bash("start_vm_sampler 'x' 100; pid=$VM_SAMPLER_PID; stop_vm_sampler; "
               "sleep 0.5; if kill -0 \"$pid\" 2>/dev/null; then echo ALIVE; else echo GONE; fi", tmp_path,
               {"LEUKOQUANT_VM_SAMPLE_S": "0.1"})
    assert "GONE" in res.stdout


@pytest.mark.skipif(not os.path.exists("/proc/self/status"), reason="needs Linux /proc")
def test_real_proc_sees_a_python_allocation(tmp_path):
    res = bash("start_vm_sampler 'alloc' $$; python3 -c 'import time; a=bytearray(300*1024*1024); time.sleep(2.5)'; stop_vm_sampler",
               "/proc", {"LEUKOQUANT_VM_SAMPLE_S": "0.5"})
    line = [l for l in res.stdout.splitlines() if l.startswith("[vm]")][0]
    assert int(line.split("VmPeak ")[1].split(" MB")[0]) >= 300


def test_run_with_peak_mem_starts_and_stops_the_sampler():
    text = BASH_UTILS.read_text()
    body = text[text.index("run_with_peak_mem() {"):text.index("# Backward compatibility no-ops")]
    assert "start_vm_sampler" in body and "stop_vm_sampler" in body


def test_stage_writes_one_usage_record_with_true_peaks(tmp_path):
    """run_with_peak_mem appends one line to every table in LEUKOQUANT_USAGE_FILES, with the context and the exact RSS."""
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    (scratch / "blob").write_bytes(b"x" * 3_000_000)
    table = tmp_path / "usage.tsv"
    res = bash("run_with_peak_mem demo bash -c 'sleep 2'", tmp_path / "noproc", {
        "LEUKOQUANT_USAGE_FILES": str(table), "LEUKOQUANT_USAGE_SCRATCH": str(scratch),
        "LEUKOQUANT_USAGE_CONTEXT": "dataset=ADNI3 subject=s1 req_mem_mb=8192 req_scratch_mb=3072"})
    assert res.returncode == 0, res.stderr
    lines = table.read_text().splitlines()
    assert lines[0] == "#schema_version=1" and len(lines) == 3
    row = dict(zip(lines[1].split("\t"), lines[2].split("\t")))
    assert row["dataset"] == "ADNI3" and row["stage"] == "demo" and row["exit_status"] == "0" and row["killed"] == "0"
    assert row["req_mem_mb"] == "8192"


@pytest.mark.parametrize("name", ["Thread Pool Worker", "a  b c", ""])
def test_scan_keeps_one_field_per_value_when_a_process_name_has_spaces_or_is_empty(tmp_path, name):
    # Real names with spaces ("Thread Pool Worker") or an empty name used to shift the whitespace-separated fields,
    # and the sampler then died with "syntax error in expression" (ADNI-3 metrics, 2026-11-20).
    make_proc(tmp_path, 100, 1, "root_shell", 500, 400)
    make_proc(tmp_path, 101, 100, name, 9_000_000, 3_000_000)
    res = bash("_vm_scan 100", tmp_path)
    assert res.returncode == 0
    fields = res.stdout.split()
    assert len(fields) == 5, res.stdout            # peak, name, hwm, name, sum: name is always exactly one word
    assert (fields[0], fields[2]) == ("9000000", "3000000")


def test_sampler_survives_repeated_samples_without_a_new_peak(tmp_path):
    # From the second sample on, a sample that is not larger than the stored peak copies the stored name. The read of
    # the stored line used to put the job-wide fields into that name, and the third sample died with "syntax error
    # in expression" (ADNI-3 metrics, 2026-11-20). A name with spaces made no difference.
    make_proc(tmp_path, 100, 1, "root_shell", 500, 400)
    make_proc(tmp_path, 101, 100, "Thread Pool Worker", 8_388_608, 2_097_152)
    res = bash("start_vm_sampler 'stage' 100; sleep 1.2; stop_vm_sampler", tmp_path, {"LEUKOQUANT_VM_SAMPLE_S": "0.2"})
    assert res.returncode == 0
    assert "syntax error" not in res.stderr, res.stderr
    line = [l for l in res.stdout.splitlines() if l.startswith("[vm]")][0]
    assert "VmPeak 8192 MB" in line
