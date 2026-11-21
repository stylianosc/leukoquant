"""Report on the per-job usage table: merge new records, flag runs close to their request, refit from true peaks.

Run on the cluster (inside the project environment):
    python3 -m leukoquant.utils.usage_report [--central TABLE] [--source FILE ...] [--glob 'PATTERN']

1. collect: merges the per-job files (logs/usage_records.tsv, default glob over the ADNI3 outputs) into the central
   table, without duplicates (see usage_record.py for the schema).
2. near misses: every record above 80 % of its memory or scratch request, or killed, with the quantity that was close.
3. refit: per job (rule, subject, session, job id) the largest peak over its stages is regressed on the DWI size in
   megavoxel-volumes, with the same line + largest-residual margin as tracula_resource_model.py. Memory uses the exact
   peak RSS, scratch the sampled peak. n, R^2, margin and the request implied at the largest observed size are printed.

Honest limits. Scratch and the address-space peaks are SAMPLED (every sample_interval_s seconds), so they are lower
bounds: a short spike is missed, and a margin fitted on them is only as good as the sampling is dense. The margin is the
largest residual of the fitted runs, so it covers every observed run but nothing larger; with few rows, or a narrow
range of sizes, it extrapolates poorly (the report prints the size range). Peak RSS is exact but is the largest SINGLE
process, while the limit acts on address space, which is larger: compare peak_vmpeak_mb with the request too.
"""
import argparse
import glob
import sys

from leukoquant.utils import usage_record as ur

DEFAULT_GLOB = "/SAN/medic/apollon/DATA/ADNI3/leukoquant/outputs/subj-*/tracula-*/logs/usage_records.tsv"
NEAR_FRACTION = 0.8


def _f(row, key):
    try:
        return float(row.get(key, ""))
    except (TypeError, ValueError):
        return None


def near_misses(rows, fraction=NEAR_FRACTION):
    """[(row, reason)] for records above `fraction` of a request, or killed."""
    out = []
    for r in rows:
        why = []
        if r.get("killed") == "1":
            why.append("killed")
        for peak, req, name in (("peak_rss_mb", "req_mem_mb", "memory (rss)"), ("peak_vmpeak_mb", "req_mem_mb", "memory (address space)"),
                                ("peak_scratch_mb", "req_scratch_mb", "scratch")):
            p, q = _f(r, peak), _f(r, req)
            if p and q and p > fraction * q:
                why.append(f"{name} {p:.0f} of {q:.0f} MB")
        if why:
            out.append((r, "; ".join(why)))
    return out


def per_job_peaks(rows):
    """One row per job: the largest peak over its stages (floats), keeping the identifying fields."""
    jobs = {}
    for r in rows:
        key = (r.get("rule"), r.get("subject"), r.get("session"), r.get("job_id"), r.get("task_id"))
        j = jobs.setdefault(key, {k: r.get(k) for k in ("dataset", "subject", "session", "rule", "dwi_mvv")})
        for k in ("peak_rss_mb", "peak_vmpeak_mb", "peak_scratch_mb"):
            v = _f(r, k)
            if v is not None and v > j.get(k, 0):
                j[k] = v
    return list(jobs.values())


def refit(rows):
    """Fit memory (GB, from peak RSS) and scratch (GB, from peak scratch) on the DWI size, from true per-job peaks."""
    from leukoquant.utils.tracula_resource_model import MARGIN_QUANTILE, fit_line   # numpy: host side only
    jobs = [j for j in per_job_peaks(rows) if _f(j, "dwi_mvv")]
    out = {}
    for name, key in (("memory_gb", "peak_rss_mb"), ("scratch_gb", "peak_scratch_mb")):
        pts = [(_f(j, "dwi_mvv"), j[key] / 1024) for j in jobs if j.get(key)]
        if len(pts) >= 3:
            out[name] = fit_line([p[0] for p in pts], [p[1] for p in pts], MARGIN_QUANTILE)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--central", default=ur.CENTRAL_TABLE)
    ap.add_argument("--source", action="append", default=[])
    ap.add_argument("--glob", default=DEFAULT_GLOB)
    args = ap.parse_args(argv)
    sources = args.source or sorted(glob.glob(args.glob))
    import os
    os.makedirs(os.path.dirname(os.path.abspath(args.central)), exist_ok=True)
    added = ur.collect(sources, args.central)
    rows = ur.read_table(args.central) if os.path.exists(args.central) else []
    print(f"collected {added} new record(s) from {len(sources)} file(s); table has {len(rows)} record(s): {args.central}")
    flagged = near_misses(rows)
    print(f"\n{len(flagged)} record(s) above {NEAR_FRACTION:.0%} of a request, or killed:")
    for r, why in flagged:
        print(f"  {r.get('dataset')} {r.get('subject')} {r.get('session')} {r.get('stage')} job {r.get('job_id')}: {why}")
    fits = refit(rows)
    if not fits:
        print("\nrefit: fewer than 3 jobs with a known DWI size, nothing fitted yet")
    for name, f in fits.items():
        print(f"\nrefit {name}: y = {f['intercept']:.2f} + {f['slope']:.4f} x Mvv + margin {f['margin']:.2f}   "
              f"n={f['n']} R2={f['r2']:.3f} largest-size {f['x_max']:.0f} Mvv "
              f"(request there: {f['intercept'] + f['slope'] * f['x_max'] + f['margin']:.1f} GB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
