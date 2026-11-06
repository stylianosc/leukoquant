"""Optional memory override for workflow steps whose request is a fixed number.

The metrics step computes its memory from the image size and has its own override (metrics_utils.py). TRACULA asks for a
fixed 8 GB, measured on EPAD where the whole job peaks near 2 GB; on other cohorts a session can need more (peaks of
4.7-8.1 GB were logged in ADNI-3 and some of those sessions were killed). LEUKOQUANT_TRACULA_MEM_GB lets one run ask for
more (or less) without a code change; unset or empty keeps the default, so nothing changes for runs that do not set it.
"""
import os

TRACULA_MEM_DEFAULT_MB = 8 * 1024
TRACULA_MEM_OVERRIDE_ENV = "LEUKOQUANT_TRACULA_MEM_GB"
TRACULA_MEM_OVERRIDE_RANGE_GB = (4, 64)


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


def tracula_mem_mb(environ=None) -> int:
    """Memory request (MB) of one TRACULA job: 8 GB unless LEUKOQUANT_TRACULA_MEM_GB sets 4-64 GB."""
    return apply_mem_override(TRACULA_MEM_DEFAULT_MB, TRACULA_MEM_OVERRIDE_ENV, TRACULA_MEM_OVERRIDE_RANGE_GB,
                              "TRACULA", environ)
