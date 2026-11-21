"""Per-job, per-stage usage records: one machine-readable line appended atomically, plus collection into a central table."""
import multiprocessing
import os
import subprocess
import sys
from pathlib import Path

import pytest

from leukoquant.utils import usage_record as ur

SCRIPT = Path(ur.__file__)


def rec(**kw):
    base = dict(dataset="ADNI3", subject="subj-1", session="", parcellation="gif", rule="trac_all_gif", stage="prep",
                job_id="7793468", task_id="2", host="node-1", dwi_dims="256x256x80x36", dwi_mvv="188.7",
                req_mem_mb="10240", req_scratch_mb="9216", peak_rss_mb="3368", peak_vmpeak_mb="5200",
                peak_vmhwm_mb="3368", peak_vm_sum_mb="6100", peak_scratch_mb="2100", min_free_scratch_mb="40000",
                elapsed_s="42100", exit_status="0", killed="0", sample_interval_s="30", n_samples="1400",
                timestamp_utc="2026-11-17T10:00:00Z")
    base.update(kw)
    return base


def test_schema_has_the_required_columns_and_a_version():
    need = {"dataset", "subject", "session", "parcellation", "rule", "stage", "job_id", "task_id", "host", "dwi_dims",
            "dwi_mvv", "req_mem_mb", "req_scratch_mb", "peak_rss_mb", "peak_vmpeak_mb", "peak_scratch_mb",
            "min_free_scratch_mb", "elapsed_s", "exit_status", "killed"}
    assert need <= set(ur.COLUMNS)
    assert isinstance(ur.SCHEMA_VERSION, int)


def test_append_writes_header_with_schema_version_once(tmp_path):
    f = tmp_path / "t.tsv"
    ur.append_record(f, rec())
    ur.append_record(f, rec(stage="bedp"))
    lines = f.read_text().splitlines()
    assert lines[0] == f"#schema_version={ur.SCHEMA_VERSION}"
    assert lines[1].split("\t") == ur.COLUMNS
    assert len(lines) == 4
    rows = ur.read_table(f)
    assert [r["stage"] for r in rows] == ["prep", "bedp"]
    assert rows[0]["peak_rss_mb"] == "3368"


def test_values_with_tabs_and_newlines_cannot_break_a_line(tmp_path):
    f = tmp_path / "t.tsv"
    ur.append_record(f, rec(host="bad\thost\nname"))
    assert len(f.read_text().splitlines()) == 3
    assert ur.read_table(f)[0]["host"] == "bad host name"


def test_missing_fields_are_empty_and_unknown_fields_are_rejected(tmp_path):
    f = tmp_path / "t.tsv"
    ur.append_record(f, {"stage": "prep"})
    assert ur.read_table(f)[0]["peak_rss_mb"] == ""
    with pytest.raises(ValueError):
        ur.append_record(f, {"nonsense": "1"})


def test_a_line_stays_below_pipe_buf(tmp_path):
    assert len(ur.format_line(rec()).encode()) < 4096


def _worker(args):
    path, i = args
    for k in range(25):
        ur.append_record(path, rec(task_id=str(i), job_id=str(k)))


def test_concurrent_appends_do_not_interleave(tmp_path):
    f = tmp_path / "t.tsv"
    with multiprocessing.Pool(6) as pool:
        pool.map(_worker, [(f, i) for i in range(6)])
    rows = ur.read_table(f)
    assert len(rows) == 150
    assert len({(r["task_id"], r["job_id"]) for r in rows}) == 150
    assert f.read_text().count("#schema_version") == 1


def test_cli_append_builds_a_record_from_key_values(tmp_path):
    f = tmp_path / "t.tsv"
    res = subprocess.run([sys.executable, "-I", str(SCRIPT), "append", "--file", str(f), "stage=path", "peak_rss_mb=10893",
                          "dataset=ADNI3", "exit_status=0"], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    r = ur.read_table(f)[0]
    assert r["stage"] == "path" and r["peak_rss_mb"] == "10893" and r["timestamp_utc"].endswith("Z")


def test_cli_never_fails_the_job_on_a_bad_target(tmp_path):
    res = subprocess.run([sys.executable, "-I", str(SCRIPT), "append", "--file", "/nonexistent/dir/t.tsv", "stage=prep"],
                         capture_output=True, text=True)
    assert res.returncode == 0
    assert "could not write" in res.stdout + res.stderr


def test_collect_merges_per_job_files_into_the_central_table_without_duplicates(tmp_path):
    a, b, central = tmp_path / "a.tsv", tmp_path / "b.tsv", tmp_path / "central.tsv"
    ur.append_record(a, rec(job_id="1", stage="prep"))
    ur.append_record(a, rec(job_id="1", stage="bedp"))
    ur.append_record(b, rec(job_id="2", stage="prep"))
    assert ur.collect([a, b], central) == 3
    assert ur.collect([a, b], central) == 0            # idempotent
    ur.append_record(b, rec(job_id="3", stage="prep"))
    assert ur.collect([a, b], central) == 1
    assert len(ur.read_table(central)) == 4


def test_context_string_for_the_rule(tmp_path):
    p = tmp_path / "d.nii.gz"
    import gzip, struct
    hdr = bytearray(348)
    struct.pack_into("<8h", hdr, 40, 4, 256, 256, 80, 36, 1, 1, 1)
    with gzip.open(p, "wb") as fh:
        fh.write(bytes(hdr))
    ctx = ur.usage_context(output_dir="/SAN/x/ADNI3/leukoquant/outputs", subject="sub-01/ses-02", parcellation="gif",
                           rule="trac_all_gif", dwi_paths=[str(p)], req_mem_mb=10240, req_scratch_mb=9216)
    kv = dict(x.split("=", 1) for x in ctx.split(" "))
    assert kv["dataset"] == "ADNI3"
    assert kv["subject"] == "sub-01" and kv["session"] == "ses-02"
    assert kv["dwi_dims"] == "256x256x80x36"
    assert float(kv["dwi_mvv"]) == pytest.approx(256 * 256 * 80 * 36 / 1e6, rel=1e-3)
    assert kv["req_mem_mb"] == "10240" and kv["req_scratch_mb"] == "9216"
    assert " " not in kv["subject"]
