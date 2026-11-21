"""Fitted TRACULA resource requests: memory and node-local scratch as continuous functions of the DWI size.

Why a fitted line and not size classes
--------------------------------------
How much memory and scratch a TRACULA run needs depends on the size of the diffusion data, not on the cohort. The size
is measured as megavoxel-volumes (Mvv) = voxels x volumes / 1e6 and read from the zip (uncompressed bytes / 2, 16-bit
DICOM) or the NIfTI header, see ``dwi_megavoxel_volumes``. For each quantity y the model is

    y(Mvv) = intercept + slope * Mvv + margin

with intercept and slope from an ordinary least-squares fit and ``margin`` an upper quantile (MARGIN_QUANTILE, 1.0 = the
largest residual) of the fit residuals, so the request covers that share of the measured runs. There is no safety factor and no size class.
The request is the prediction rounded up to whole GB and clamped only at the two ends (see mem_override.py).

Three quantities are fitted (see ``select_rows`` for which runs feed each):
  scratch_gb    node-local scratch in use at the END of a completed run (``du`` of the job scratch folder, at cleanup).
                The peak during the run is higher and unknown: mem_override.py multiplies this fit by a provisional factor
                until the per-stage [node] df lines give real peaks, then this fit should be redone from those.
  prep_rss_gb   peak resident memory of ``trac-all -prep`` (any run that went on to bedpostx)
  path_rss_gb   peak resident memory of ``trac-all -path`` (completed runs only: a failed stage is cut short)
The memory request is the larger of the prep and path predictions (bedpostx stays near 2.1 GB at every size).

Data and limits (be honest when quoting the fit)
------------------------------------------------
Input is ``tracula_model/tracula_job_fit_input.tsv`` (one row per session x parcellation, pseudonymous ids removed), produced from
the run logs by ``export_fit_input``. ``python -m leukoquant.utils.tracula_resource_model`` refits and rewrites
``tracula_model/tracula_resource_model.json``. What the data cannot give:
  * SGE keeps no accounting here (no qacct), so there is no peak virtual memory (maxvmem) of finished jobs. The request
    is a limit on address space (h_vmem); the fit is on resident memory (RSS). One running job showed maxvmem 3.4 GB against
    3.37 GB RSS, but the two are not proven equal.
  * Peak-RSS lines exist only for runs since 2026-09-30, about 150 rows, almost all of them small images. Few runs of
    256x256 data have a memory line, one of them a 10.9 GB outlier (``trac-all -path``, 54 volumes), which sets the margin
    of the path fit.
  * Runs killed at their cap show the cap, not the need; they are excluded where the stage did not complete.

Why the margin is the largest residual and the 10.9 GB run stays in (checked 2026-11-17)
----------------------------------------------------------------------------------------
Path-stage RSS (GB) predicted at 189 Mvv (256x256x80x36) and 283 Mvv (256x256x80x54) by each variant of the fit:
    all runs, margin = largest residual (used)    9.2 and 10.6
    all runs, margin = 0.99 residual quantile      4.5 and  6.0
    without the 10.9 GB run (any margin)           2.9 and  3.4   (quantile regression at 0.95 without it: 2.7 and 3.1)
Only one path measurement exists above 138 Mvv, the 10.9 GB run; the 149 others are small images. Without it the model
asks for the 8 GB floor for 256x256 data, and the direct evidence says that is too small: trac-paths failed with
bad_alloc at the 8 GB cap for four 256x256x80x36 sessions (021-s-7045/7055/7062/7092, peak RSS 7996 MB = the cap), and
that 10.9 GB run completed. So dropping it, or taking a lower residual quantile, repeats a known failure; the variant
used is the only one consistent with both observations. It is evidence-limited, not well determined: it rests on one
point and should be refitted when the sampled VmPeak/VmHWM lines ([vm]) and more large-image runs exist.
"""
import csv
import gzip
import json
import math
import os
import struct
import sys
import zipfile
from functools import lru_cache
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parent / "tracula_model"   # not named "data": .gitignore excludes that
FIT_INPUT_PATH = DATA_DIR / "tracula_job_fit_input.tsv"
MODEL_PATH = DATA_DIR / "tracula_resource_model.json"

# Upper quantile of the fit residuals added as margin. 1.0 = the largest residual, so the request covers every measured
# run (in-sample): an under-sized request kills a run that may already have computed for days, an over-sized one only
# reserves resources. With 150 memory rows a lower quantile would sit between the two largest residuals anyway.
MARGIN_QUANTILE = 1.0

# Header volume counts of 72 belong to an older run than the log that was parsed (two scans stacked in June, the rerun
# used one): their size does not match their memory lines.
STALE_HEADER_VOLUMES = 72

FIT_INPUT_COLUMNS = ["dataset", "parc", "dim_x", "dim_y", "dim_z", "n_vol", "prep_rss_mb", "bedp_rss_mb", "path_rss_mb",
                     "scratch_gb_at_cleanup", "completed"]


# ---------------------------------------------------------------------------------------------------------------- size
def dwi_megavoxel_volumes(paths):
    """Voxels x volumes (in millions) of the DWI inputs, or None if no input can be sized.

    A zip of DICOMs is sized by its uncompressed bytes / 2 (16-bit data; a header read would need a conversion), a NIfTI
    by its header. Several inputs add up, which over-estimates a session whose duplicate scan is dropped later."""
    return _megavoxel_volumes(tuple(str(p) for p in (paths or [])))


@lru_cache(maxsize=4096)
def _megavoxel_volumes(paths):
    total, found = 0.0, False
    for path in paths:
        try:
            if path.endswith(".zip"):
                with zipfile.ZipFile(path) as z:
                    total += sum(i.file_size for i in z.infolist()) / 2 / 1e6
                found = True
            elif path.endswith((".nii", ".nii.gz")):
                opener = gzip.open if path.endswith(".gz") else open
                with opener(path, "rb") as f:
                    dim = struct.unpack("<8h", f.read(348)[40:56])
                n = 1
                for d in dim[1:1 + max(dim[0], 1)]:
                    n *= max(d, 1)
                total += n / 1e6
                found = True
        except (OSError, zipfile.BadZipFile, struct.error):
            continue
    return total if found else None


# ---------------------------------------------------------------------------------------------------------------- fit
def fit_line(x, y, quantile=MARGIN_QUANTILE):
    """Least-squares line y = intercept + slope * x and the upper `quantile` of its residuals as margin."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 3:
        raise ValueError("need at least 3 points to fit a line")
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (intercept + slope * x)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float((resid ** 2).sum()) / ss_tot if ss_tot > 0 else 0.0
    return {"intercept": float(intercept), "slope": float(slope), "margin": float(np.quantile(resid, quantile)),
            "n": int(x.size), "r2": float(r2), "quantile": float(quantile),
            "residual_max": float(resid.max()), "x_max": float(x.max())}


def _num(value):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def select_rows(rows):
    """Split table rows into the data of each fitted quantity: {name: [{"mvv", "y", "completed"}, ...]}."""
    out = {"scratch_gb": [], "prep_rss_gb": [], "path_rss_gb": []}
    for r in rows:
        dims = [_num(r.get(k)) for k in ("dim_x", "dim_y", "dim_z", "n_vol")]
        if any(d is None or d <= 0 for d in dims) or dims[3] == STALE_HEADER_VOLUMES:
            continue
        mvv = dims[0] * dims[1] * dims[2] * dims[3] / 1e6
        done = int(_num(r.get("completed")) or 0) == 1
        scratch, prep, bedp, path = (_num(r.get(k)) for k in
                                     ("scratch_gb_at_cleanup", "prep_rss_mb", "bedp_rss_mb", "path_rss_mb"))
        if done and scratch is not None and scratch > 0:
            out["scratch_gb"].append({"mvv": mvv, "y": scratch, "completed": 1})
        if prep is not None and bedp is not None:           # prep finished: the run went on to bedpostx
            out["prep_rss_gb"].append({"mvv": mvv, "y": prep / 1024, "completed": int(done)})
        if done and path is not None:
            out["path_rss_gb"].append({"mvv": mvv, "y": path / 1024, "completed": 1})
    return out


def build_model(rows, quantile=MARGIN_QUANTILE):
    """The three fits plus provenance, from table rows (dicts with the FIT_INPUT_COLUMNS)."""
    sel = select_rows(rows)
    model = {name: fit_line([r["mvv"] for r in data], [r["y"] for r in data], quantile) for name, data in sel.items()}
    model["limits"] = {"memory_gb_min": 8, "memory_gb_max": 16, "memory_override_gb_max": 24,
                       "scratch_gb_min": 3, "scratch_gb_max": 40}
    model["rows"] = {name: len(data) for name, data in sel.items()}
    model["note"] = ("y = intercept + slope * Mvv + margin; margin = residual quantile; fitted on resident memory (RSS) "
                     "and scratch use at the end of the run, no maxvmem available; see module docstring")
    return model


def predict_gb(fit, mvv):
    """Prediction in GB of one fitted quantity at `mvv` megavoxel-volumes, margin included."""
    return fit["intercept"] + fit["slope"] * mvv + fit["margin"]


def request_gb(fit, mvv, low, high):
    """Whole-GB request: the prediction rounded up, clamped at the two ends only."""
    return int(min(max(math.ceil(predict_gb(fit, mvv)), low), high))


@lru_cache(maxsize=1)
def load_model():
    with open(MODEL_PATH) as f:
        return json.load(f)


# -------------------------------------------------------------------------------------------------------- data files
def export_fit_input(table_path, out_path=FIT_INPUT_PATH):
    """Reduce the run table (tracula_job_table.tsv from the cluster logs) to the columns the fit needs, without ids."""
    with open(table_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(FIT_INPUT_COLUMNS)
        for r in rows:
            w.writerow([r.get(c, "") for c in FIT_INPUT_COLUMNS])
    return len(rows)


def main(argv=None):
    """Refit from the committed fit input and rewrite the model JSON. With a table path: re-export the fit input first."""
    argv = sys.argv[1:] if argv is None else argv
    if argv:
        print(f"exported {export_fit_input(argv[0])} rows to {FIT_INPUT_PATH}")
    with open(FIT_INPUT_PATH, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    model = build_model(rows, MARGIN_QUANTILE)
    with open(MODEL_PATH, "w") as f:
        json.dump(model, f, indent=2, sort_keys=True)
        f.write("\n")
    for name in ("scratch_gb", "prep_rss_gb", "path_rss_gb"):
        m = model[name]
        print(f"{name:12s} n={m['n']:5d}  y = {m['intercept']:.3f} + {m['slope']:.5f} * Mvv + margin {m['margin']:.3f}"
              f"  (R2 {m['r2']:.2f}, quantile {m['quantile']}, max residual {m['residual_max']:.2f})")


if __name__ == "__main__":
    main()
