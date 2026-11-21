"""Per-job, per-stage usage records: the measured peaks of memory and scratch, one line per stage, in one table.

Why
---
Memory and scratch requests (tracula_resource_model.py) were fitted to what the logs happened to show: resident memory
and scratch at the END of a run. The limits act on something else (address space, scratch at its PEAK). Each heavy stage
now writes one record with the true peaks, so requests are driven by the right data.

What is recorded (schema COLUMNS, version SCHEMA_VERSION) and how it is measured
--------------------------------------------------------------------------------
  peak_rss_mb          exact: ru_maxrss of the stage (GNU time %M), the largest resident set of any ONE process of the stage.
  peak_vmpeak_mb       sampled: largest VmPeak (peak address space) of any ONE process. THIS is the unit ``ulimit -v`` /
                       SGE h_vmem acts on (a per-process address-space limit).
  peak_vmhwm_mb        sampled: largest VmHWM (peak resident) of any one process (cross-check of peak_rss_mb).
  peak_vm_sum_mb       sampled: largest SUM of VmSize over all processes of the stage at one moment, the quantity SGE adds
                       up for its job-wide vmem accounting (maxvmem). qstat/qacct cannot be used from inside the container.
  peak_scratch_mb      sampled: largest ``du`` of the job's scratch folder.
  min_free_scratch_mb  sampled: smallest free space (df) of the filesystem holding the scratch folder, at the same moments:
                       a near-full node shows here before it kills a job.
  elapsed_s, exit_status, killed (1 = exit 137 or an out-of-memory message), host, job and task id, requests.

Limits of sampling (read before trusting a margin)
--------------------------------------------------
The sampled columns are sampled every ``sample_interval_s`` seconds (recorded in the row with ``n_samples``). They are
never above the true peak and can be below it: a spike shorter than the interval, or a process that starts and ends
between two samples, is missed. Fast allocations (an image loaded and freed in seconds) are the likely misses, so a
sampled peak is a lower bound; peak_rss_mb is exact and a sampled VmHWM far below it shows how much was missed. A margin
fitted on sampled peaks should include the observed gap between sampled VmHWM and exact RSS.

Atomicity
---------
A line is written with ONE write() on an O_APPEND descriptor under an exclusive flock, well below PIPE_BUF (4096 bytes).
Inside the TRACULA container only the output folders are bound, so the record goes to the subject's logs folder
(``usage_records.tsv``); ``collect`` merges those files into the central table
(/SAN/medic/apollon/leukoquant/pipeline_status/usage/usage_table.tsv) without duplicates. Writing never fails the job.

The module needs only the standard library (it runs inside the container); ``usage_context`` imports numpy-based helpers
lazily and is called on the host when the workflow is built.
"""
import csv
import fcntl
import gzip
import os
import struct
import sys
import time

SCHEMA_VERSION = 1
CENTRAL_TABLE = "/SAN/medic/apollon/leukoquant/pipeline_status/usage/usage_table.tsv"

COLUMNS = ["dataset", "subject", "session", "parcellation", "rule", "stage", "job_id", "task_id", "host", "dwi_dims",
           "dwi_mvv", "req_mem_mb", "req_scratch_mb", "peak_rss_mb", "peak_vmpeak_mb", "peak_vmhwm_mb",
           "peak_vm_sum_mb", "peak_scratch_mb", "min_free_scratch_mb", "elapsed_s", "exit_status", "killed",
           "sample_interval_s", "n_samples", "timestamp_utc"]

# A record is identified by where and when it was taken; collect uses this to avoid duplicates.
_KEY = ("dataset", "subject", "session", "parcellation", "rule", "stage", "job_id", "task_id", "timestamp_utc")


def _clean(value):
    return " ".join(str(value).split())          # no tab, newline or repeated blank can break the line


def format_line(record):
    unknown = set(record) - set(COLUMNS)
    if unknown:
        raise ValueError(f"unknown field(s): {sorted(unknown)}")
    return "\t".join(_clean(record.get(c, "")) for c in COLUMNS) + "\n"


def _header():
    return f"#schema_version={SCHEMA_VERSION}\n" + "\t".join(COLUMNS) + "\n"


def append_record(path, record):
    """Append one record (header first when the file is new) with a single write under an exclusive lock."""
    line = format_line(record).encode()
    fd = os.open(str(path), os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        data = line if os.fstat(fd).st_size > 0 else _header().encode() + line
        os.write(fd, data)
    finally:
        try:
            fcntl.flock(fd, fcntl.LOCK_UN)
        finally:
            os.close(fd)


def read_table(path):
    """The records of a table as dicts of strings (comment line and header skipped)."""
    with open(path, newline="") as f:
        lines = [l for l in f if not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def collect(sources, central):
    """Add the records of the per-job files `sources` to the central table; returns how many were new."""
    seen = set()
    if os.path.exists(central):
        seen = {tuple(r.get(k, "") for k in _KEY) for r in read_table(central)}
    added = 0
    for src in sources:
        if not os.path.exists(src):
            continue
        for r in read_table(src):
            key = tuple(r.get(k, "") for k in _KEY)
            if key in seen:
                continue
            append_record(central, {c: r.get(c, "") for c in COLUMNS})
            seen.add(key)
            added += 1
    return added


def dwi_dims(paths):
    """'XxYxZxT' of the first NIfTI among `paths`, else '' (a zip of DICOMs has no header without a conversion)."""
    for p in paths or []:
        p = str(p)
        if p.endswith((".nii", ".nii.gz")):
            try:
                with (gzip.open if p.endswith(".gz") else open)(p, "rb") as f:
                    d = struct.unpack("<8h", f.read(348)[40:56])
                return "x".join(str(v) for v in d[1:1 + max(d[0], 1)])
            except (OSError, struct.error):
                return ""
    return ""


def usage_context(output_dir, subject, parcellation, rule, dwi_paths, req_mem_mb, req_scratch_mb):
    """Space-separated key=value string the rule exports so each record knows what it belongs to (no spaces inside)."""
    from leukoquant.utils.tracula_resource_model import dwi_megavoxel_volumes   # numpy: host side only
    parts = [p for p in str(output_dir).split("/") if p]
    dataset = parts[parts.index("leukoquant") - 1] if "leukoquant" in parts and parts.index("leukoquant") > 0 else ""
    subj, _, sess = str(subject).partition("/")
    mvv = dwi_megavoxel_volumes(dwi_paths)
    kv = {"dataset": dataset, "subject": subj, "session": sess, "parcellation": parcellation, "rule": rule,
          "dwi_dims": dwi_dims(dwi_paths), "dwi_mvv": "" if mvv is None else f"{mvv:.3f}",
          "req_mem_mb": int(req_mem_mb), "req_scratch_mb": int(req_scratch_mb)}
    return " ".join(f"{k}={_clean(v).replace(' ', '_')}" for k, v in kv.items())


def _cli(argv):
    """python3 usage_record.py append --file F [--file F2 ...] key=value ...   (never fails: a job must not die on this)"""
    if not argv or argv[0] != "append":
        print("usage: usage_record.py append --file FILE [--file FILE ...] key=value ...")
        return 0
    files, record, i = [], {}, 1
    while i < len(argv):
        if argv[i] == "--file" and i + 1 < len(argv):
            files.append(argv[i + 1])
            i += 2
            continue
        if "=" in argv[i]:
            k, v = argv[i].split("=", 1)
            record[k] = v
        i += 1
    record.setdefault("timestamp_utc", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    for f in files:
        try:
            append_record(f, record)
        except (OSError, ValueError) as e:
            print(f"[usage] could not write the usage record to {f}: {e}")
    return 0


if __name__ == "__main__":
    sys.exit(_cli(sys.argv[1:]))
