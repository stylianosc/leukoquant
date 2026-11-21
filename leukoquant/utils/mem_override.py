"""Resource requests of the TRACULA job: memory and node-local scratch from the DWI size, with optional overrides.

The metrics step computes its memory from the image size and has its own override (metrics_utils.py). TRACULA used a flat
8 GB; peaks of 4.7-10.9 GB were logged in ADNI-3 (large 256x256 DWI). Both requests are now fitted lines of the DWI size
(tracula_resource_model.py: y = intercept + slope x megavoxel-volumes + margin, coefficients in tracula_model/*.json, refit with
``python -m leukoquant.utils.tracula_resource_model``), rounded up to whole GB and clamped only at the two ends:

  memory   8 to 16 GB (8 GB is also the request when the DWI size is unknown)
  scratch  3 to 40 GB (3 GB when the DWI size is unknown)

Scratch is PROVISIONAL: the fit is on the scratch in use when a run ended (``du`` at cleanup), but the peak during the run
is higher and unknown, and exceeding the request can kill a job on a node without free space. The fitted value is
therefore multiplied by TRACULA_SCRATCH_PEAK_FACTOR (2.0, a stated guess, not a measurement; about 10 GB for 256x256
data). Refit it from the per-stage ``[node]`` df lines of the TRACULA log (utils/bash_utils.sh log_node_disk_state) as soon
as runs with them exist, then replace the fit input and drop the factor.

LEUKOQUANT_TRACULA_MEM_GB (4-24 GB, the user's cap) and LEUKOQUANT_TRACULA_SCRATCH_GB (3-40 GB) override a request for one
run; unset or empty uses the fitted value.
"""
import math
import os

from leukoquant.utils.tracula_resource_model import dwi_megavoxel_volumes, load_model, predict_gb

TRACULA_MEM_MIN_GB, TRACULA_MEM_MAX_GB = 8, 16            # clamp of the fitted memory request
TRACULA_MEM_DEFAULT_MB = TRACULA_MEM_MIN_GB * 1024        # request when the DWI size is unknown
TRACULA_MEM_OVERRIDE_ENV = "LEUKOQUANT_TRACULA_MEM_GB"
TRACULA_MEM_OVERRIDE_RANGE_GB = (4, 24)                   # 24 GB is the user's cap for any TRACULA run

# Node-local scratch requested by the TRACULA job (SGE tscratch; a scheduling number, not a quota).
TRACULA_SCRATCH_DEFAULT_MB = 3 * 1024
TRACULA_SCRATCH_OVERRIDE_ENV = "LEUKOQUANT_TRACULA_SCRATCH_GB"
TRACULA_SCRATCH_OVERRIDE_RANGE_GB = (3, 40)
# Provisional headroom of the peak scratch over the scratch measured at cleanup (see the module docstring).
TRACULA_SCRATCH_PEAK_FACTOR = 2.0


def apply_mem_override(default_mb: int, env_name: str, range_gb: tuple, label: str, environ=None) -> int:
    """Memory request (MB) after the optional override named env_name; it must be a whole number of GB inside range_gb.

    Anything else is ignored with a printed warning and the default is kept."""
    environ = os.environ if environ is None else environ
    raw = environ.get(env_name)
    if raw is None or str(raw).strip() == "":
        return default_mb
    low, high = range_gb
    try:
        gb = int(str(raw).strip())
    except ValueError:
        gb = None
    if gb is None or not low <= gb <= high:
        print(f"[{label} memory] ignoring {env_name}={raw!r}: need a whole number of GB between {low} and {high}; "
              f"using the default {default_mb // 1024} GB", flush=True)
        return default_mb
    print(f"[{label} memory] override {gb} GB (default {default_mb // 1024} GB)", flush=True)
    return gb * 1024


def tracula_mem_mb(environ=None, dwi_paths=None) -> int:
    """Memory request (MB) of one TRACULA job.

    LEUKOQUANT_TRACULA_MEM_GB (4-24) wins. Otherwise the larger of the fitted prep and path peaks at this DWI size, rounded
    up to whole GB and clamped to 8-16 GB; 8 GB when the size is unknown."""
    environ = os.environ if environ is None else environ
    if str(environ.get(TRACULA_MEM_OVERRIDE_ENV, "")).strip() != "":
        return apply_mem_override(TRACULA_MEM_DEFAULT_MB, TRACULA_MEM_OVERRIDE_ENV, TRACULA_MEM_OVERRIDE_RANGE_GB,
                                  "TRACULA", environ)
    mvv = dwi_megavoxel_volumes(dwi_paths)
    if mvv is None:
        return TRACULA_MEM_DEFAULT_MB
    model = load_model()
    peak = max(predict_gb(model["prep_rss_gb"], mvv), predict_gb(model["path_rss_gb"], mvv))
    return int(min(max(math.ceil(peak), TRACULA_MEM_MIN_GB), TRACULA_MEM_MAX_GB)) * 1024


def tracula_scratch_mb(dwi_paths=None, environ=None) -> int:
    """Scratch request (MB) of one TRACULA job.

    LEUKOQUANT_TRACULA_SCRATCH_GB (3-40) wins. Otherwise the fitted scratch use at cleanup at this DWI size (line plus
    residual margin) times the provisional peak factor, rounded up to whole GB and clamped to 3-40 GB; 3 GB when the size
    is unknown."""
    environ = os.environ if environ is None else environ
    low, high = TRACULA_SCRATCH_OVERRIDE_RANGE_GB
    if str(environ.get(TRACULA_SCRATCH_OVERRIDE_ENV, "")).strip() != "":
        return apply_mem_override(TRACULA_SCRATCH_DEFAULT_MB, TRACULA_SCRATCH_OVERRIDE_ENV,
                                  TRACULA_SCRATCH_OVERRIDE_RANGE_GB, "TRACULA scratch", environ)
    mvv = dwi_megavoxel_volumes(dwi_paths)
    if mvv is None:
        return TRACULA_SCRATCH_DEFAULT_MB
    model = load_model()
    fit = model["scratch_gb"]
    # The factor only compensates a fit made on scratch measured at cleanup; a fit on true sampled peaks
    # (usage_report.py refit, "scratch_basis": "peak" in the model JSON) is used as it is.
    factor = 1.0 if model.get("scratch_basis") == "peak" else TRACULA_SCRATCH_PEAK_FACTOR
    gb = math.ceil(factor * predict_gb(fit, mvv))
    return int(min(max(gb, low), high)) * 1024
