"""Memory sizing for the metrics rule (used at Snakemake workflow parse time).

The metrics step runs ``mri_synthstrip`` on the T1 (twice: the registration T1 and the TIV mask). The tool resamples
the image to 1 mm isotropic voxels and pads every axis up to a multiple of 64, so its memory scales with the physical
size of the T1, not with its voxel count. Measured on the cluster (12 T1s, 6 to 24 GB limits, 1 GB steps):

    peak RSS            ~= 0.33 + 0.43 * padded_voxels_M   [GB]
    needed address space ~= 4.0 + 0.40 * padded_voxels_M   [GB]

SGE's h_vmem limits the address space, not the resident set, and the address space is 2.3 to 3.2 GB above the
peak RSS, so the address-space line is the one the request must cover. Sessions with large T1s (21 M and 26 M padded
voxels: ADNI-3 and EPAD at 1.2 x 1.05 x 1.05 mm, OASIS-3 at 1 mm) failed at 12 GB and needed up to 14 GB.

The line above is a conservative upper bound. The 2026-10 rerun then ran the sessions that had failed at 12 GB
at 14 GB: every 26.2 M session passed (ADNI-3 peak RSS 11.4 GB, EPAD 9.3 GB; the line predicted 14.5 GB of address
space), and the 21 M sessions that failed at 12 GB pass at 14 GB. So the request follows what was verified:

    padded size below 20 M      -> 12 GB (the floor; all such sessions passed at 12 GB)
    20 M up to 26.3 M           -> 14 GB (the largest size seen in ADNI-3, EPAD, OASIS-3)
    above 26.3 M (unverified)   -> the line above + the margin, so an unseen larger image is not under-requested

Like ``compute_zscore_mem_mb`` all metrics jobs of one run share one request, sized for the LARGEST T1 of the
dataset (each dataset has its own submission).
"""
import glob
import math
import os

# Fit of the address space the skull strip needs: base + slope * padded voxels (millions of 1 mm voxels), in GB.
METRICS_MEM_BASE_GB = 4.0
METRICS_MEM_GB_PER_MVOXEL = 0.40
# A fixed margin on top of the measured need (node-to-node variation of the limit, which the 12 samples cannot
# bound), instead of a multiple: 2 GB over the need of the largest T1. The need is first rounded to a whole GB.
METRICS_MEM_MARGIN_GB = 2
METRICS_MEM_VERIFIED_FROM_MVOXEL = 20.0   # from here, the verified request applies (12 GB failed at 21 M)
METRICS_MEM_VERIFIED_MB = 14 * 1024       # passed for every session up to METRICS_MEM_VERIFIED_UP_TO_MVOXEL
METRICS_MEM_VERIFIED_UP_TO_MVOXEL = 26.3  # the largest padded size in the data (196x256x256 at 1.2/1.05/1.05)
METRICS_MEM_FALLBACK_MB = 16 * 1024  # used when no T1 header can be read: the request the largest known T1 gets, not the floor
METRICS_MEM_FLOOR_MB = 12 * 1024    # the rest of the step (registration, metrics tables) peaked at 6.7-9.7 GB RSS
METRICS_MEM_CAP_MB = 32 * 1024
_PAD = 64                           # mri_synthstrip pads each axis to a multiple of 64 voxels


def padded_mvoxels(shape, zooms) -> float:
    """Millions of 1 mm voxels of the T1 after the skull strip's resampling and padding."""
    n = 1
    for size, zoom in zip(shape[:3], zooms[:3]):
        n *= math.ceil(size * zoom / _PAD) * _PAD
    return n / 1e6


def resolve_t1_candidates(spec, subject):
    """Concrete T1 files for one subject from a path spec, without failing on anything unexpected.

    *spec* is either a file path, or the "base:glob" form of the metrics config, where the image is looked up as
    ``base/<subject>/<glob>`` (the rule's own ``_sif_t1_path`` builds the same path inside the container)."""
    if not spec:
        return []
    spec = str(spec).replace("{subject}", subject)
    base, _, pattern = spec.partition(":")
    if not pattern.strip(":"):
        return [base] if os.path.isfile(base) else []
    return sorted(glob.glob(os.path.join(base, subject, pattern.split(":")[0])))


def compute_metrics_mem_mb(t1_paths=None) -> int:
    """Memory request (MB) for the metrics jobs, sized for the largest T1 among *t1_paths* (a path or a list).

    Reads NIfTI headers only. A missing or unreadable path is skipped with a printed warning; if no header can be
    read the fallback (sized for the largest known T1) is returned, so a path problem cannot silently under-request. The need is rounded to a whole GB and the fixed margin is added; never below the floor, never above the cap."""
    if isinstance(t1_paths, str):
        t1_paths = [t1_paths]
    largest, unreadable = 0.0, 0
    for path in t1_paths or []:
        if not (path and os.path.isfile(path)):
            unreadable += 1
            continue
        try:
            import nibabel as nib
            header = nib.load(path).header
            largest = max(largest, padded_mvoxels(header.get_data_shape(), header.get_zooms()))
        except Exception:
            unreadable += 1
    if unreadable:
        print(f"[metrics memory] {unreadable} T1 path(s) missing or unreadable; skipped", flush=True)
    if not largest:
        print(f"[metrics memory] no T1 header readable; using the fallback of {METRICS_MEM_FALLBACK_MB // 1024} GB", flush=True)
        return METRICS_MEM_FALLBACK_MB
    if largest <= METRICS_MEM_VERIFIED_UP_TO_MVOXEL:
        # Sizes with a verified request: below the threshold the floor applies, from it the verified 14 GB.
        mem_mb = METRICS_MEM_VERIFIED_MB if largest >= METRICS_MEM_VERIFIED_FROM_MVOXEL else METRICS_MEM_FLOOR_MB
    else:
        needed_gb = METRICS_MEM_BASE_GB + METRICS_MEM_GB_PER_MVOXEL * largest
        mem_mb = (math.floor(needed_gb + 0.5) + METRICS_MEM_MARGIN_GB) * 1024  # whole GB: need rounded, plus the margin
    return min(max(mem_mb, METRICS_MEM_FLOOR_MB), METRICS_MEM_CAP_MB)


METRICS_MEM_OVERRIDE_ENV = "LEUKOQUANT_METRICS_MEM_GB"
METRICS_MEM_OVERRIDE_RANGE_GB = (8, 32)


def apply_metrics_mem_override(computed_mb: int, environ=None) -> int:
    """Memory request (MB) for the metrics jobs after the optional LEUKOQUANT_METRICS_MEM_GB override.

    The override lets a run use less (or more) than the computed value, e.g. to try 12 GB first for sessions whose
    failure had nothing to do with memory. It must be a whole number of GB between 8 and 32; anything else is
    ignored with a printed warning and the computed value is kept."""
    environ = os.environ if environ is None else environ
    raw = environ.get(METRICS_MEM_OVERRIDE_ENV)
    if raw is None or str(raw).strip() == "":
        return computed_mb
    low, high = METRICS_MEM_OVERRIDE_RANGE_GB
    try:
        gb = int(str(raw).strip())
    except ValueError:
        gb = None
    if gb is None or not low <= gb <= high:
        print(f"[metrics memory] ignoring {METRICS_MEM_OVERRIDE_ENV}={raw!r}: need a whole number of GB between "
              f"{low} and {high}; using the computed {computed_mb // 1024} GB", flush=True)
        return computed_mb
    print(f"[metrics memory] override {gb} GB (computed {computed_mb // 1024} GB)", flush=True)
    return gb * 1024

