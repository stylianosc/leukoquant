#!/usr/bin/env python3
"""Metrics extraction for brain imaging (CSVD white matter damage).

All inputs are expected to be in T1 space (resampled upstream by metrics_calc.sh).
Metrics are organised into families (one function per family).  Results are written
as a single CSV row per subject so that many subjects can be concatenated later.

CLI
---
--subject         Subject ID (written into every output row).
--t1-path         T1 file / directory / glob (already in T1 space).
--lesion-path     Lesion mask file / directory / glob (optional).
--t1-brain-path   Skull-stripped T1 brain image used for TIV estimation (optional).
--skip-skullstrip If set, TIV is not computed even if a brain file is supplied.
--maps            Diffusion / other map spec: name;path,name2;path2 (optional).
--tract-mode      Tract analysis mode: tractography, tractography-atlas (default), or atlas.
--output          Output CSV path (default: metrics.csv).
--verbose         Enable verbose logging.
"""

from __future__ import annotations

import warnings
import logging
# Suppress verbose/repetitive nibabel warnings (e.g. qfac/pixdim[0] warnings)
logging.getLogger('nibabel').setLevel(logging.ERROR)
warnings.filterwarnings("ignore", category=UserWarning, module="nibabel")
# Set root logger to WARNING by default (only show errors/warnings unless --verbose)
logging.getLogger().setLevel(logging.WARNING)

import argparse
import glob
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, cast

# Allow imports of leukoquant when run from containers where the module is
# bind-mounted at /leukoquant but not installed in site-packages.
_leukoquant_paths = ["/leukoquant", "/leukoquant"]
for path in _leukoquant_paths:
    if path not in sys.path and os.path.isdir(path):
        sys.path.insert(0, path)

from scipy import ndimage
from scipy.stats import rankdata
from skimage.morphology import skeletonize

from leukoquant.utils.tract_qc_util import qc_tract_file as _run_qc_tract_file

import importlib.util as _importlib_util


def _load_optional_extension():
    """Load leukoquant/utils/extra/extended_analysis.py if present.

    That path is not part of the public distribution -- this is a generic
    extension point for additional per-subject analysis modules dropped in
    locally. Returns None (silent no-op everywhere it's checked) when the
    file doesn't exist, which is always true for a public installation.
    """
    ext_path = Path(__file__).parent / "extra" / "extended_analysis.py"
    if not ext_path.is_file():
        return None
    spec = _importlib_util.spec_from_file_location("_leukoquant_extended_analysis", ext_path)
    if spec is None or spec.loader is None:
        return None
    module = _importlib_util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_extension = _load_optional_extension()

import nibabel as nib
from nibabel.nifti1 import Nifti1Image
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ============================================================================
# Path utilities
# ============================================================================

def _load_qc_report(qc_report_path: str) -> Dict[str, bool]:
    """Load a tract_qc qc_report.csv and return a dict of tract_dir_name → qc_pass.

    The CSV 'file' column contains the full path to path.pd.trk; the tract name
    is the parent directory (e.g. acomm_avg16_syn_bbr). The caller matches on
    this directory name, which is also how metrics_calc.py names tracts from the
    .nii.gz side. Maps both the full tract_dir name and base name (stripped of
    _avg16 suffix) so lookups by either key succeed. Also checks 'tract_name'
    if present in the CSV. Returns an empty dict if the file cannot be read.
    """
    try:
        df = pd.read_csv(qc_report_path)
        result: Dict[str, bool] = {}
        for _, row in df.iterrows():
            file_str = str(row.get("file", ""))
            tract_name_str = str(row.get("tract_name", ""))
            qc_pass = bool(row.get("qc_pass", False))
            if file_str and file_str != "nan":
                tract_dir = str(Path(file_str).parent.name)
                result[tract_dir] = qc_pass
                result[tract_dir.split("_avg16", 1)[0]] = qc_pass
            if tract_name_str and tract_name_str != "nan":
                result[tract_name_str] = qc_pass
                result[tract_name_str.split("_avg16", 1)[0]] = qc_pass
        logger.info("Loaded QC report from %s: %d tracts", qc_report_path, len(result))
        return result
    except Exception as exc:
        logger.warning("Could not load QC report %s: %s - falling back to inline QC", qc_report_path, exc)
        return {}


def _qc_pass_for_tract(
    tract_dir_name: str,
    tract_file: str,
    qc_map: Dict[str, bool],
) -> Optional[bool]:
    """Return the QC pass/fail for a tract.

    Looks up ``tract_dir_name`` or its base tract name in the pre-loaded ``qc_map``.
    If not found (report absent or tract missing from it) falls back to running QC
    inline on the sibling tracts.trk file. Returns None when QC cannot be determined.
    """
    if tract_dir_name in qc_map:
        return qc_map[tract_dir_name]
    base_name = tract_dir_name.split("_avg16", 1)[0]
    if base_name in qc_map:
        return qc_map[base_name]
    # Fallback: run QC inline (covers the case where no report was passed)
    trk_file = str(Path(tract_file).parent / "tracts.trk")
    if Path(trk_file).exists():
        result = _run_qc_tract_file(trk_file)
        return result.get("qc_pass")
    return None


def _expand_path_spec(spec: str) -> List[str]:
    """Expand a file path, directory, or glob pattern into a sorted file list."""
    p = Path(spec)
    if p.exists() and p.is_file():
        return [str(p.resolve())]
    if p.exists() and p.is_dir():
        files = sorted(
            [str(f.resolve()) for f in p.rglob("*.nii")]
            + [str(f.resolve()) for f in p.rglob("*.nii.gz")]
        )
        return files
    matches = sorted(str(Path(x).resolve()) for x in glob.glob(spec, recursive=True))
    return [f for f in matches if Path(f).is_file()]

def _parse_lesion_spec(spec: str, default_name: str = "lesion") -> List[Tuple[str, str]]:
    """Parse lesion spec into list of (name, path) tuples.

    Format (comma-separated entries)::

        [name=]path or [name:]path

    ``name`` is separated from the path spec by ``=`` or ``:`` (optional, default ``lesion``).
    """
    entries: List[Tuple[str, str]] = []
    if not spec:
        return entries
    for token in [t.strip() for t in spec.split(",") if t.strip()]:
        if "=" in token:
            name, path = token.split("=", 1)
            name, path = name.strip(), path.strip()
        elif ":" in token:
            name, path = token.split(":", 1)
            name, path = name.strip(), path.strip()
        else:
            name, path = default_name, token.strip()
        if not path:
            raise ValueError(f"Lesion spec '{token}': path must be non-empty.")
        entries.append((name, path))
    return entries

def _parse_named_spec(spec: str, default_name: str = "") -> List[Tuple[str, str]]:
    """Parse named specs into ``(name, path)`` tuples.

    Accepted forms (comma-separated):

    - ``name=path`` (preferred, shell-safe)
    - ``path`` (uses ``default_name`` if provided)
    
    If name not provided and default_name given, use default_name.
    Otherwise raise error if name required.
    """
    entries: List[Tuple[str, str]] = []
    if not spec:
        return entries
    for token in [t.strip() for t in spec.split(",") if t.strip()]:
        if "=" in token:
            name, path = token.split("=", 1)
            name, path = name.strip(), path.strip()
        else:
            name, path = default_name, token.strip()
        if not path:
            raise ValueError(f"Spec '{token}': path must be non-empty.")
        if not name:
            raise ValueError(f"Spec '{token}': name must be non-empty (or provide default).")
        entries.append((name, path))
    return entries

def _load(path: str) -> Any:
    return cast(Any, nib).load(path)

def _load_array(path: str) -> np.ndarray:
    """Load a NIfTI data array as a C-contiguous ndarray.

    nibabel returns Fortran-ordered data; every downstream step works on C-order
    flat voxel indices (so that gathered values reach np.sum / np.percentile in
    the same order as the original ``values[np.nonzero(values)]``), and a
    C-contiguous array lets ``arr.reshape(-1)`` be a free view.
    """
    return _treat_non_finite_as_missing(np.ascontiguousarray(np.asarray(_load(path).dataobj)), path)


def _treat_non_finite_as_missing(arr: np.ndarray, label: str = "") -> np.ndarray:
    """Set NaN and +-Inf voxels to 0, the value that every summary excludes (a voxel with no data).

    A NaN left in a map would turn every statistic of every region it touches into NaN. Treating such voxels as missing
    keeps all results finite; the number of voxels affected is logged.
    """
    if arr.dtype.kind not in "fc":
        return arr
    bad = ~np.isfinite(arr)
    n_bad = int(bad.sum())
    if n_bad:
        logger.warning("%d non-finite voxel(s) in %s treated as missing (0).", n_bad, label or "an input array")
        arr = np.where(bad, 0, arr).astype(arr.dtype, copy=False)
    return arr

def _voxel_vol_mm3(img: Any) -> float:
    """Return voxel volume in mm³ from a nibabel image."""
    return float(np.prod([float(z) for z in img.header.get_zooms()[:3]]))

def _assert_shape_compatibility(
    reference_shape: Tuple[int, ...],
    arrays: Dict[str, np.ndarray],
    group_name: str,
) -> None:
    """Raise ValueError if any array in a group does not match reference shape."""
    mismatches: List[str] = []
    for name, arr in arrays.items():
        if tuple(arr.shape) != tuple(reference_shape):
            mismatches.append(
                f"{group_name} '{name}': expected {reference_shape}, got {tuple(arr.shape)}"
            )
    if mismatches:
        details = "\n".join(mismatches)
        raise ValueError(f"Shape mismatch detected:\n{details}")

# ============================================================================
# Map-based metric utilities
# ============================================================================

def _compute_load(mask_array: np.ndarray, region_array: np.ndarray, voxel_vol_mm3: float) -> float:
    """Compute load (intersection volume) between two binary masks.

    Parameters
    ----------
    mask_array : np.ndarray
        Binary mask (lesion, etc.) with dtype allowing numeric comparison.
    region_array : np.ndarray
        Binary region mask (tract, whole brain, etc.).
    voxel_vol_mm3 : float
        Single voxel volume in mm³.

    Returns
    -------
    float
        Intersection volume in mm³.
    """
    intersection = np.logical_and(mask_array > 0, region_array > 0)
    n_voxels = int(np.sum(intersection))
    return float(n_voxels) * voxel_vol_mm3

def write_csv(df: pd.DataFrame, path: str) -> None:
    """
    Write DataFrame to CSV with logging.

    Written to a temp path first, then atomically renamed into place
    (os.replace, same guarantee as POSIX rename) -- if the process gets
    killed mid-write (walltime, OOM, an unrelated `qdel`), the final path
    either still holds the last fully-valid CSV or never existed, never a
    truncated one. Paired with metrics_workflow.smk's update() flag on
    these outputs: together, a subject that already has valid metrics is
    never left with nothing just because Snakemake schedules the rule to
    run again for an unrelated reason.
    """
    if df.empty:
        logger.warning("DataFrame is empty; no data to write.")
        return

    # "_tmp" spliced in before the extension (not appended after, e.g. not
    # "name.csv.tmp") -- pandas.to_csv() writes to the literal path given
    # with no extension-sniffing, so this isn't strictly required the way
    # it was for fslmaths, but keeping the same convention everywhere
    # avoids ever having to re-derive that reasoning per call site.
    root, ext = os.path.splitext(path)
    tmp_path = f"{root}_tmp{ext}"
    df.to_csv(tmp_path, index=False)
    os.replace(tmp_path, path)
    logger.info("DataFrame written to: %s", path)

def _compute_tiv_from_array(brain_array: np.ndarray, voxel_vol_mm3: float) -> float:
    """Estimate total intracranial volume from a skull-stripped brain array.

    Parameters
    ----------
    brain_array : np.ndarray
        Skull-stripped T1 brain image data.
    voxel_vol_mm3 : float
        Single voxel volume in mm³.

    Returns
    -------
    float
        TIV in mm³.
    """
    n_voxels = int(np.sum(brain_array > 0))
    tiv = float(n_voxels) * voxel_vol_mm3
    logger.info("TIV: %d voxels × %.4f mm³ = %.1f mm³", n_voxels, voxel_vol_mm3, tiv)
    return tiv

# Distribution-shape statistics. They are summary statistics of every region like the others, but a
# lesion-versus-outside ratio of two signed shape values (_summary_change_nonzero) is not meaningful.
SHAPE_STATS = ("skew", "kurtosis")


def _shape_stats(values: np.ndarray, mean: float) -> Tuple[float, float]:
    """Skewness and excess kurtosis of a 1-D sample (moment estimators).

    Same definitions as scipy.stats.skew / scipy.stats.kurtosis with their defaults (bias=True,
    fisher=True): skew = m3 / m2^1.5 and kurtosis = m4 / m2^2 - 3, with m_k the k-th central moment.
    A sample with no spread has no shape: when the variance is zero, or lies below the floating-point
    resolution of the mean (scipy's own guard, where it returns NaN with a precision-loss warning),
    both are reported as 0.0, the value that a symmetric, normal-tailed sample would have.
    """
    d = values - mean
    d2 = d * d
    m2 = float(np.mean(d2))
    if m2 <= (np.finfo(np.float64).resolution * mean) ** 2 or m2 == 0.0:
        return 0.0, 0.0
    m3 = float(np.mean(d2 * d))
    m4 = float(np.mean(d2 * d2))
    return m3 / m2 ** 1.5, m4 / (m2 * m2) - 3.0


# Statistics of every region summary, in column order. `n_voxels` (the number of non-zero values) is the first key
# so that a statistic of an empty region (all 0.0) can be told apart from a real zero, and so that the reliability of a
# statistic (a standard deviation of one voxel is 0) is visible to a model.
SUMMARY_STATS = ("p5", "25th", "median", "75th", "p95", "mean", "sum", "std", "IQR", "peak_width", "skew", "kurtosis")


def _summary_nonzero(values: np.ndarray, count_key: Optional[str] = "n_voxels") -> Dict[str, float]:
    """Return summary stats over non-zero values (zeros if empty).

    Keys: n_voxels (number of non-zero values; omitted when `count_key` is None), p5, 25th, median, 75th, p95
    (percentiles), mean, sum, std, IQR (75th - 25th), peak_width (95th - 5th percentile), skew and kurtosis
    (excess kurtosis; see _shape_stats). An empty region gives 0.0 everywhere (never NaN).
    """
    non_zero = values[np.nonzero(values)]
    out: Dict[str, float] = {}
    if count_key:
        out[count_key] = float(non_zero.size)
    if non_zero.size == 0:
        out.update({stat: 0.0 for stat in SUMMARY_STATS})
        return out

    # One np.percentile call (a single partition) instead of five; the values
    # are bit-identical to separate calls (same linear interpolation).
    p5, p25, p50, p75, p95 = (float(x) for x in np.percentile(non_zero, [5, 25, 50, 75, 95]))
    # Float64 for the moments, so float32 maps do not lose precision in the 3rd and 4th powers.
    non_zero64 = non_zero.astype(np.float64, copy=False)
    mean = float(np.mean(non_zero64))
    skew, kurtosis = _shape_stats(non_zero64, mean)
    out.update({
        "p5": p5,
        "25th": p25,
        "median": p50,
        "75th": p75,
        "p95": p95,
        "mean": float(np.mean(non_zero)),
        "sum": float(np.sum(non_zero)),
        "std": float(np.std(non_zero)),
        "IQR": p75 - p25,
        "peak_width": p95 - p5,
        "skew": skew,
        "kurtosis": kurtosis,
    })
    return out

def _weighted_percentiles(values: np.ndarray, weights: np.ndarray, quantiles) -> List[float]:
    """Weighted percentiles (quantiles in [0, 1]) by linear interpolation of the weighted cumulative distribution.

    Each value sits at the middle of its own weight interval, so with equal weights the result is close to (not identical
    with) np.percentile; _weighted_summary uses np.percentile itself when the weights are all equal.
    """
    order = np.argsort(values, kind="stable")
    v = values[order]
    w = weights[order]
    positions = (np.cumsum(w) - 0.5 * w) / w.sum()
    return [float(x) for x in np.interp(quantiles, positions, v)]


def _weighted_summary(values: np.ndarray, weights: np.ndarray, count_key: Optional[str] = "n_voxels") -> Dict[str, float]:
    """The summary statistics of `values` where each voxel counts in proportion to its `weights`.

    Use this instead of summarising `values * weights`: that product shrinks a value towards 0 wherever its weight is small
    (for a tract probability map, at the tract's edge), so the statistics describe the weights as much as the map. Here the
    values are untouched and the weights decide how much each voxel counts: weighted mean sum(w v) / sum(w), weighted
    standard deviation, weighted percentiles, weighted skewness and excess kurtosis (central moments with the same
    weighted mean). `sum` stays sum(w v), the weighted total (an expected volume for probability weights).

    Voxels where the value or the weight is 0 are excluded, exactly as the zeros of a product are excluded by
    _summary_nonzero. With equal weights the percentiles are np.percentile's and the mean, std, skewness and kurtosis
    equal the unweighted ones. An empty region gives 0.0 everywhere.
    """
    keep = (values != 0) & (weights != 0)
    v = values[keep].astype(np.float64, copy=False)
    w = weights[keep].astype(np.float64, copy=False)
    out: Dict[str, float] = {}
    if count_key:
        out[count_key] = float(v.size)
    if v.size == 0:
        out.update({stat: 0.0 for stat in SUMMARY_STATS})
        return out

    total = float(w.sum())
    mean = float((w * v).sum() / total)
    d = v - mean
    d2 = d * d
    m2 = float((w * d2).sum() / total)
    if m2 <= (np.finfo(np.float64).resolution * mean) ** 2 or m2 == 0.0:
        skew, kurtosis = 0.0, 0.0
    else:
        m3 = float((w * d2 * d).sum() / total)
        m4 = float((w * d2 * d2).sum() / total)
        skew, kurtosis = m3 / m2 ** 1.5, m4 / (m2 * m2) - 3.0
    if np.all(w == w[0]):
        p5, p25, p50, p75, p95 = (float(x) for x in np.percentile(v, [5, 25, 50, 75, 95]))
    else:
        p5, p25, p50, p75, p95 = _weighted_percentiles(v, w, [0.05, 0.25, 0.5, 0.75, 0.95])
    out.update({
        "p5": p5, "25th": p25, "median": p50, "75th": p75, "p95": p95,
        "mean": mean, "sum": float((w * v).sum()), "std": float(np.sqrt(m2)),
        "IQR": p75 - p25, "peak_width": p95 - p5, "skew": skew, "kurtosis": kurtosis,
    })
    return out


# ============================================================================
# Metric family: whole-brain
# ============================================================================

def _compute_whole_brain_metrics(
    subject: str,
    brain_array: Optional[np.ndarray],
    lesion_dict: Dict[str, np.ndarray],
    voxel_vol_mm3: float,
) -> pd.DataFrame:
    """Compute whole-brain family metrics from numpy arrays, per lesion type.

    Columns
    -------
    subject            Subject ID.
    tiv_mm3            Total intracranial volume (mm³).
    {lesion_name}_volume_mm3
                       Lesion volume in whole brain (per lesion type).
    {lesion_name}_load_pct
                       Lesion load % in whole brain (per lesion type): volume / TIV x 100.
    {lesion_name}_n_lesions
                       Number of lesions (26-connected components) per lesion type.
    {lesion_name}_mean_lesion_volume_mm3
                       Mean volume of a lesion (mm³); 0.0 when there is no lesion.
    """
    row: Dict[str, Any] = {"subject": subject}

    tiv_mm3: Optional[float] = None
    if brain_array is not None:
        tiv_mm3 = _compute_tiv_from_array(brain_array, voxel_vol_mm3)
    else:
        logger.warning("No brain array available; tiv_mm3 will be None.")
    row["tiv_mm3"] = tiv_mm3

    # Compute metrics per lesion type.
    for lesion_name, lesion_array in lesion_dict.items():
        if brain_array is not None:
            lesion_vol = _compute_load(lesion_array, brain_array, voxel_vol_mm3)
        else:
            lesion_vol = float(np.sum(lesion_array > 0)) * voxel_vol_mm3
        row[f"{lesion_name}_volume_mm3"] = lesion_vol
        if tiv_mm3 and tiv_mm3 > 0:
            row[f"{lesion_name}_load_pct"] = lesion_vol / tiv_mm3 * 100.0
        else:
            row[f"{lesion_name}_load_pct"] = None
        _, n_lesions = ndimage.label(lesion_array > 0, structure=STRUCTURE_26)
        row[f"{lesion_name}_n_lesions"] = int(n_lesions)
        row[f"{lesion_name}_mean_lesion_volume_mm3"] = lesion_vol / n_lesions if n_lesions else 0.0

    return pd.DataFrame([row])

# ============================================================================
# Metric family: periventricular versus deep white matter lesions
# ============================================================================

# Lesions within this distance of the ventricles are periventricular, the rest deep. The two have different causes and
# are graded separately on the Fazekas scale.
PERIVENTRICULAR_MM = 10.0


def _ventricle_mask(parcellation_path: str, labels: str, ref_img: Any) -> np.ndarray:
    """Binary ventricle mask on the reference (T1) voxel grid, from a label image.

    The label image is resampled to the reference grid with nearest-neighbour interpolation when its shape or affine
    differ (a FreeSurfer aparc+aseg is conformed, a GIF parcellation is usually already on the T1 grid), then the voxels
    whose label is in `labels` ("50,51,52,53") are kept. An empty mask is an error: it would silently make every lesion
    "deep".
    """
    import nibabel.processing as nibproc
    wanted = [int(x) for x in str(labels).split(",") if x.strip()]
    img = _load(parcellation_path)
    if tuple(img.shape[:3]) != tuple(ref_img.shape[:3]) or not np.allclose(img.affine, ref_img.affine, atol=1e-3):
        img = nibproc.resample_from_to(img, (ref_img.shape[:3], ref_img.affine), order=0)
    mask = np.isin(np.asarray(img.dataobj), wanted)
    if not mask.any():
        raise ValueError(f"None of the ventricle labels {wanted} is present in {parcellation_path}")
    return mask.astype(np.uint8)


def _periventricular_deep_metrics(
    lesion_name: str,
    lesion_array: np.ndarray,
    ventricle_array: np.ndarray,
    voxel_sizes,
    voxel_vol_mm3: float,
) -> Dict[str, float]:
    """Lesion volume within PERIVENTRICULAR_MM of the ventricles, beyond it, and the periventricular share.

    {lesion}_pvwmh_volume_mm3   lesion volume (mm3) within PERIVENTRICULAR_MM of the ventricular system
    {lesion}_dwmh_volume_mm3    lesion volume (mm3) farther away
    {lesion}_pvwmh_fraction     pvwmh / (pvwmh + dwmh); 0.0 for a subject without lesion

    `ventricle_array` is a binary mask (non-zero = ventricle) on the lesion's voxel grid.
    """
    lesion = np.asarray(lesion_array) > 0
    ventricles = np.asarray(ventricle_array) > 0
    if lesion.shape != ventricles.shape:
        raise ValueError(f"Ventricle mask shape {ventricles.shape} does not match the lesion shape {lesion.shape}")
    if not lesion.any():
        pv = dw = 0.0
    else:
        near = _distance_to_lesion(ventricles, voxel_sizes) <= PERIVENTRICULAR_MM if ventricles.any() else np.zeros_like(lesion)
        pv = float(np.count_nonzero(lesion & near)) * voxel_vol_mm3
        dw = float(np.count_nonzero(lesion & ~near)) * voxel_vol_mm3
    total = pv + dw
    return {
        f"{lesion_name}_pvwmh_volume_mm3": pv,
        f"{lesion_name}_dwmh_volume_mm3": dw,
        f"{lesion_name}_pvwmh_fraction": pv / total if total > 0 else 0.0,
    }


# ============================================================================
# Metric family: whole-brain lesion microstructure
# ============================================================================

# ============================================================================
# Metric family: tract-lesion load
# ============================================================================

# ============================================================================
# SPARSE IMPLEMENTATIONS of the per-tract / per-lesion statistics.
#
# Every statistic below is _summary_nonzero(product of a few arrays), i.e. it only
# ever looks at the NON-ZERO entries of an elementwise product. A product can only
# be non-zero where each factor is non-zero, so it is enough to evaluate the same
# elementwise expression on the voxels where the sparsest factor is non-zero (a
# tract's support, a lesion's support ...), taken in ascending C-order flat index.
# Non-finite voxels are set to 0 when an array is read (_treat_non_finite_as_missing),
# so no NaN/Inf can reach these functions.
# ============================================================================

# Distance-based tissue around lesions. The penumbra is the ring of tissue within PENUMBRA_RINGS_MM[0] of a lesion; the
# dilated lesion is the lesion plus that ring. Two outer rings show how far the change in a map extends beyond the lesion.
# The normal-appearing white matter (NAWM) of a tract is the part of the tract outside the dilated lesion: the reference
# tissue that the lesion is compared with (the tissue just beside a lesion is already partly damaged, so "everything outside
# the lesion" is not a clean reference).
# Lesion components are 26-connected (faces, edges and corners), here and in the per-lesion table.
STRUCTURE_26 = np.ones((3, 3, 3), dtype=bool)

PENUMBRA_RINGS_MM = ((0.0, 3.0), (3.0, 6.0), (6.0, 9.0))
RING_SUFFIXES = ("penumbra", "penumbra_3_6", "penumbra_6_9")
DILATION_MM = PENUMBRA_RINGS_MM[0][1]
# Distance to the nearest lesion is capped (and is this value when a subject has no lesion at all).
DISTANCE_CAP_MM = 100.0


def _distance_to_lesion(lesion_mask: np.ndarray, voxel_sizes) -> np.ndarray:
    """Euclidean distance (mm, float32) from every voxel to the nearest lesion voxel, capped at DISTANCE_CAP_MM.

    scipy's distance transform of an array without any lesion voxel returns distances to the array border, which are
    meaningless here; a subject without a lesion therefore gets DISTANCE_CAP_MM everywhere.
    """
    mask = np.asarray(lesion_mask) != 0
    if not mask.any():
        return np.full(mask.shape, DISTANCE_CAP_MM, dtype=np.float32)
    dist = ndimage.distance_transform_edt(~mask, sampling=voxel_sizes)
    return np.minimum(dist, DISTANCE_CAP_MM).astype(np.float32)


def _summary_reduced(values: np.ndarray) -> Dict[str, float]:
    """The short summary used for the outer rings: n_voxels, p5, median, p95, mean, std (zeros if empty)."""
    full = _summary_nonzero(values)
    return {k: full[k] for k in ("n_voxels", "p5", "median", "p95", "mean", "std")}


def _flat(arr: np.ndarray) -> np.ndarray:
    """C-order 1-D view of an array (a copy only if it is not C-contiguous)."""
    return arr.reshape(-1)


class _LesionContext:
    """Per-lesion-set data computed once and reused for every tract and map."""

    def __init__(self, lesion_array, binary, dilated, penumbra, distance=None):
        self.lesion = _flat(lesion_array)
        self.binary = _flat(binary)
        self.dilated = _flat(dilated)
        self.penumbra = _flat(penumbra)
        self.shape = lesion_array.shape
        # supports (ascending C-order flat indices) of every region used
        self.idx_pd = np.flatnonzero(self.lesion)
        self.idx_bin = np.flatnonzero(self.binary)
        self.idx_dil = np.flatnonzero(self.dilated)
        self.idx_pen = np.flatnonzero(self.penumbra)
        # lesion_map of _summary_change_nonzero / distance code: np.where(x < 0.5, 0, 1)
        self.lm = ~(self.lesion < 0.5)
        self.has_negative = bool((self.lesion < 0).any())
        # 26-connected components of the binary lesion
        labels, self.n_comp = ndimage.label(binary, structure=STRUCTURE_26)
        self.labels = _flat(labels)
        self.labels_at_bin = self.labels[self.idx_bin]
        # components used by _compute_lesion_tract_intersection_details: label
        # of lesion_binary derived from lesion_array. For the arrays built by
        # _load_inputs this is the same mask, so reuse the labelling.
        if lesion_array.dtype in [np.float32, np.float64]:
            det_binary = ~(self.lesion < 0.5)
        else:
            det_binary = self.lesion > 0
        if np.array_equal(det_binary, self.binary != 0):
            self.det_labels, self.det_n = self.labels, self.n_comp
        else:
            dl, self.det_n = ndimage.label(det_binary.reshape(self.shape), structure=STRUCTURE_26)
            self.det_labels = _flat(dl)
        self.det_comp_voxels = np.bincount(self.det_labels, minlength=self.det_n + 1)[1:]
        self._dist = None if distance is None else _flat(np.asarray(distance, dtype=np.float32))

    def distance_flat(self, voxel_sizes):
        """Distance (mm, capped; see _distance_to_lesion) to the nearest lesion voxel, once per lesion set."""
        if self._dist is None:
            self._dist = _flat(_distance_to_lesion(self.lm.reshape(self.shape), voxel_sizes))
        return self._dist


class _TractContext:
    """Per-tract supports, computed once and reused for every lesion and map."""

    def __init__(self, tract_array, tract_binary, tract_skeleton):
        self.shape = tract_array.shape
        flat = _flat(tract_array)
        self.flat = flat
        self.idx_nz = np.flatnonzero(flat)            # where the tract is non-zero
        self.vals = flat[self.idx_nz]                  # tract values there
        self.binary = _flat(tract_binary)
        self.idx_bin = np.flatnonzero(self.binary)     # tract > 0
        self.idx_skel = np.flatnonzero(_flat(tract_skeleton) > 0.5)
        self.has_negative = bool((self.vals < 0).any())
        # lesion-independent products with each map, filled lazily
        self.map_pd = {}


# Statistics compared between the lesion and the reference tissue (no n_voxels, no shape statistics: the difference of two
# skewness values is not a meaningful quantity).
CHANGE_STATS = ("p5", "25th", "median", "75th", "p95", "mean", "sum", "std", "IQR", "peak_width")


def _auc(inside: np.ndarray, outside: np.ndarray) -> float:
    """Probability that a random `inside` value exceeds a random `outside` value (ties count one half).

    This is the Mann-Whitney U statistic divided by the number of pairs. It lies in [0, 1] whatever the values are (0.5 = no
    difference, 1 = every inside value is larger), needs no distribution, and a single extreme value hardly moves it. If
    either group is empty there is nothing to compare and the result is the neutral 0.5.
    """
    n1, n2 = len(inside), len(outside)
    if n1 == 0 or n2 == 0:
        return 0.5
    ranks = rankdata(np.concatenate([inside, outside]))          # average ranks: ties count one half
    return float((ranks[:n1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n2))


def _lesion_vs_reference(inside: np.ndarray, reference: np.ndarray) -> Dict[str, float]:
    """How the values inside the lesion differ from those of the reference tissue (the tract's NAWM).

    For each statistic S in CHANGE_STATS: `{S}_diff` = S(inside) - S(reference), in the units of the map (for a Z-score map,
    standard deviations). Plus `auc` (see _auc). Zero values are excluded, as in every region summary. Where a group is
    empty the differences are 0.0 and the AUC 0.5: no values, no detectable change, never NaN.
    """
    inside = inside[inside != 0].astype(np.float64, copy=False)
    reference = reference[reference != 0].astype(np.float64, copy=False)
    a = _summary_nonzero(inside, count_key=None)
    b = _summary_nonzero(reference, count_key=None)
    both = len(inside) > 0 and len(reference) > 0
    out = {f"{stat}_diff": (a[stat] - b[stat]) if both else 0.0 for stat in CHANGE_STATS}
    out["auc"] = _auc(inside, reference)
    return out


def _fast_lesion_tract_intersection_details(lc: _LesionContext, tc: _TractContext, voxel_vol_mm3: float):
    """Fast equivalent of _compute_lesion_tract_intersection_details.

    The connected components of the lesion mask and their volumes are
    tract-independent (computed once in lc); per tract only the components'
    overlap counts are needed: one bincount of the component labels over the
    tract's voxels.
    """
    num_components = lc.det_n
    component_volumes_mm3 = np.asarray(lc.det_comp_voxels, dtype=float) * voxel_vol_mm3

    result: Dict[str, Any] = {
        "in_tract_mm3": 0.0,
        "load_in_tract_pct": 0.0,
        "overall_percentage_overlap": 0.0,
        "intersecting_lesions_total_volume": 0.0,
        "number_lesions_tract_overlap": 0,
    }

    tract_volume_mm3 = float(len(tc.idx_bin)) * voxel_vol_mm3
    # _compute_load(lesion_array, tract_array): voxels with lesion > 0 and tract > 0
    n_load = int(np.count_nonzero(lc.lesion[tc.idx_bin] > 0))
    tract_lesion_load_mm3 = float(n_load) * voxel_vol_mm3
    result["in_tract_mm3"] = tract_lesion_load_mm3
    if tract_volume_mm3 > 0:
        result["load_in_tract_pct"] = tract_lesion_load_mm3 / tract_volume_mm3 * 100.0

    lesion_sizes: List[float] = []
    overlap_percentages: List[float] = []
    overlap_volumes: List[float] = []
    intersecting_total_volume = 0.0
    intersecting_components = 0

    # voxels of every lesion component inside this tract, one bincount
    overlap_counts = np.bincount(lc.det_labels[tc.idx_bin], minlength=num_components + 1)
    for component_id in np.flatnonzero(overlap_counts[1:]) + 1:
        component_volume = float(component_volumes_mm3[component_id - 1])
        overlap_volume_count = int(overlap_counts[component_id]) * voxel_vol_mm3
        if overlap_volume_count > 0:
            lesion_sizes.append(component_volume)
            overlap_percentages.append(float(overlap_volume_count) / component_volume * 100.0)
            overlap_volumes.append(overlap_volume_count)
            intersecting_total_volume += component_volume
            intersecting_components += 1

    lesion_sizes_array = np.array(lesion_sizes)
    overlap_percentages_array = np.array(overlap_percentages)
    overlap_volumes_array = np.array(overlap_volumes)

    total_overlap_voxels = int(overlap_counts[1:].sum())
    total_overlap_volume = total_overlap_voxels * voxel_vol_mm3
    if total_overlap_volume != overlap_volumes_array.sum():
        logger.warning(
            "Total overlap volume (%.2f mm³) does not match sum of component overlaps (%.2f mm³).",
            total_overlap_volume, overlap_volumes_array.sum()
        )

    result["overall_percentage_overlap"] = (
        float(total_overlap_volume) / intersecting_total_volume * 100.0
        if intersecting_total_volume > 0 else 0.0
    )
    result["intersecting_lesions_total_volume"] = intersecting_total_volume
    result["number_lesions_tract_overlap"] = intersecting_components

    distribution_inputs = (
        (lesion_sizes_array, "lesion_sizes_mm3"),
        (overlap_percentages_array, "percentage_tract_lesion_overlap"),
        (overlap_volumes_array, "lesions_overlap_volumes_mm3"),
    )
    for values, suffix in distribution_inputs:
        stats = _summary_nonzero(values, count_key=None)
        for stat_name, stat_value in stats.items():
            result[f"{stat_name}_{suffix}"] = stat_value

    return result


FRACTION_WITHIN_MM = 5.0


def _fast_tract_distance_metrics(lc: _LesionContext, tc: _TractContext, voxel_sizes):
    """How close the tract is to the lesions (mm; capped at DISTANCE_CAP_MM, which is also the value without any lesion).

    tract_nearest_lesion_distance            smallest distance of any tract voxel to a lesion (0 if the tract contains lesion)
    tract_mean_lesion_distance               mean over the tract's voxels of the distance to the nearest lesion
    tract_fraction_within_5mm_of_lesion      share of the tract's voxels that lie within 5 mm of a lesion (lesion included)

    The distance transform depends only on the lesion mask, so it is computed once per lesion set (lc.distance_flat).
    These replace the distance of the tract's centre of mass, which for a curved tract can lie outside the tract.
    A tract without voxels has no distance: the cap, and a fraction of 0.
    """
    if len(tc.idx_bin) == 0:
        return {
            "tract_nearest_lesion_distance": DISTANCE_CAP_MM,
            "tract_mean_lesion_distance": DISTANCE_CAP_MM,
            f"tract_fraction_within_{int(FRACTION_WITHIN_MM)}mm_of_lesion": 0.0,
        }
    d = lc.distance_flat(voxel_sizes)[tc.idx_bin].astype(np.float64)
    return {
        "tract_nearest_lesion_distance": float(d.min()),
        "tract_mean_lesion_distance": float(d.mean()),
        f"tract_fraction_within_{int(FRACTION_WITHIN_MM)}mm_of_lesion": float(np.mean(d <= FRACTION_WITHIN_MM)),
    }


def _compute_tract_level_metrics_fast(
    subject: str,
    tract_arrays: Dict[str, np.ndarray],
    tract_binary_arrays: Dict[str, np.ndarray],
    tract_skeletons_arrays: Dict[str, np.ndarray],
    tract_source: Dict[str, str],
    lesion_arrays: Dict[str, np.ndarray],
    binary_lesion_arrays: Dict[str, np.ndarray],
    dilated_lesion_arrays: Dict[str, np.ndarray],
    penumbras_arrays: Dict[str, np.ndarray],
    map_arrays: Dict[str, np.ndarray],
    voxel_vol: float,
    voxel_sizes,
    lesion_contexts: Optional[Dict[str, _LesionContext]] = None,
) -> pd.DataFrame:
    """Support-based implementation of the tract-level metrics loop.

    Same rows, same columns in the same order, same values as the dense loop
    (_compute_tract_level_metrics_legacy). See the FAST PATH note above.
    """
    if lesion_contexts is None:
        lesion_contexts = {
            name: _LesionContext(lesion_arrays[name], binary_lesion_arrays[name],
                                 dilated_lesion_arrays[name], penumbras_arrays[name])
            for name in lesion_arrays
        }
    map_flat = {name: _flat(arr) for name, arr in map_arrays.items()}
    map_negative = {name: bool((arr < 0).any()) for name, arr in map_arrays.items()}
    rows: List[pd.DataFrame] = []

    for tract_name, tract_array in tract_arrays.items():
        source = tract_source.get(tract_name, "tractography")
        tc = _TractContext(tract_array, tract_binary_arrays[tract_name], tract_skeletons_arrays[tract_name])
        tbin_flat = tc.binary
        tract_row = {
            "subject": subject,
            "tract_name": tract_name,
            "tract_source": source,
            "used_atlas": source == "atlas",
            "tract_volume_mm3": float(len(tc.idx_bin)) * voxel_vol,
        }

        if tc.has_negative:
            logger.warning(f"Warning: Tract '{tract_name}' contains negative values; check input data.")


        for lesion_name, lc in lesion_contexts.items():
            if lc.has_negative:
                logger.warning(f"Warning: Lesion '{lesion_name}' contains negative values; check input data.")

            tract_lesion_load_details = _fast_lesion_tract_intersection_details(lc, tc, voxel_vol)
            tract_lesion_distance_metrics = _fast_tract_distance_metrics(lc, tc, voxel_sizes)

            # Region products, evaluated on each region's own support. The factor
            # order and dtypes are those of the dense expressions.
            t_pd = tc.flat
            tpl_pd = t_pd[lc.idx_pd] * lc.lesion[lc.idx_pd]                       # tract PD x lesion PD
            tpl_bin = t_pd[lc.idx_bin] * lc.binary[lc.idx_bin]                    # tract PD x lesion binary
            tpl_dil = t_pd[lc.idx_dil] * lc.dilated[lc.idx_dil]                   # tract PD x dilated lesion
            tpl_pen = t_pd[lc.idx_pen] * lc.penumbra[lc.idx_pen]                  # tract PD x penumbra
            tbl_bin = tbin_flat[lc.idx_bin] * lc.binary[lc.idx_bin]               # tract binary x lesion binary
            tbl_dil = tbin_flat[lc.idx_dil] * lc.dilated[lc.idx_dil]              # tract binary x dilated
            tbl_pen = tbin_flat[lc.idx_pen] * lc.penumbra[lc.idx_pen]             # tract binary x penumbra

            # Components of the lesion mask that touch this tract (tract binary x lesion binary)
            overlap_component_ids = np.unique(lc.labels_at_bin[tbl_bin != 0])
            overlap_component_ids = overlap_component_ids[overlap_component_ids != 0]
            idx_whole = lc.idx_bin[np.isin(lc.labels_at_bin, overlap_component_ids)]
            n_whole = len(idx_whole)

            stats_tract_pd_lesion_pd = _summary_nonzero(tpl_pd)
            stats_tract_pd_lesion_binary = _summary_nonzero(tpl_bin)
            stats_tract_pd_lesion_dilated = _summary_nonzero(tpl_dil)
            stats_tract_pd_lesion_penumbra = _summary_nonzero(tpl_pen)

            # Tissue by distance to the nearest lesion: the NAWM (the tract outside the dilated lesion) and the two outer rings.
            dist_all = lc.distance_flat(voxel_sizes)
            d_nz, d_bin = dist_all[tc.idx_nz], dist_all[tc.idx_bin]
            nawm_nz = d_nz > DILATION_MM                                # over the tract's probability support
            nawm_bin = tc.idx_bin[d_bin > DILATION_MM]                  # over the binary tract
            outer_rings = {
                suffix: tc.idx_bin[(d_bin > lo) & (d_bin <= hi)]
                for suffix, (lo, hi) in zip(RING_SUFFIXES[1:], PENUMBRA_RINGS_MM[1:])
            }
            # Share of the tract's streamline density that lies inside lesions: a voxel-based proxy for how much of the tract
            # a lesion cuts (the tract probability map is the normalised streamline density). Not the share of streamlines
            # touching a lesion, which needs the streamlines themselves.
            density_total = float(tc.vals.sum())
            tract_row[f"{lesion_name}_tract_density_fraction_in_lesion"] = (
                float(tc.vals[lc.lm[tc.idx_nz]].sum()) / density_total if density_total > 0 else 0.0)
            tract_row[f"{lesion_name}_tract_nawm_n_voxels"] = float(len(nawm_bin))
            for suffix, idx in outer_rings.items():
                tract_row[f"{lesion_name}_tract_{suffix}_n_voxels"] = float(len(idx))

            # Lesion versus the tract's NAWM, for the tract probability itself and (below) for every map.
            inside_nz = lc.lm[tc.idx_nz]
            change_tract_pd_lesion = _lesion_vs_reference(tc.vals[inside_nz], tc.vals[nawm_nz])
            inside_bin = tc.idx_bin[lc.lm[tc.idx_bin]]

            tract_lesion_metrics = {
                "tract_pd_lesion_pd": stats_tract_pd_lesion_pd,
                "tract_pd_lesion_binary": stats_tract_pd_lesion_binary,
                "tract_pd_lesion_dilated": stats_tract_pd_lesion_dilated,
                "tract_pd_lesion_penumbra": stats_tract_pd_lesion_penumbra,
                "tract_lesion_vs_nawm": change_tract_pd_lesion,
            }

            for detail_key, detail_val in tract_lesion_load_details.items():
                tract_row[f"{lesion_name}_{detail_key}"] = detail_val
            for distance_key, distance_val in tract_lesion_distance_metrics.items():
                tract_row[f"{lesion_name}_{distance_key}"] = distance_val
            for suffix, metrics in tract_lesion_metrics.items():
                for metric_name, metric_value in metrics.items():
                    tract_row[f"{lesion_name}_{suffix}_{metric_name}"] = metric_value

            ones_whole = np.ones(n_whole, dtype=np.uint8)
            for map_name, mflat in map_flat.items():
                stats_map_tract_pd_lesion_pd = _weighted_summary(mflat[lc.idx_pd], tpl_pd)
                stats_map_tract_binary_lesion_binary = _summary_nonzero(mflat[lc.idx_bin] * tbl_bin)
                stats_map_tract_pd_dilated_lesion = _weighted_summary(mflat[lc.idx_dil], tpl_dil)
                stats_map_tract_binary_dilated_lesion = _summary_nonzero(mflat[lc.idx_dil] * tbl_dil)
                stats_map_tract_pd_penumbra = _weighted_summary(mflat[lc.idx_pen], tpl_pen)
                stats_map_tract_binary_penumbra = _summary_nonzero(mflat[lc.idx_pen] * tbl_pen)
                stats_map_tract_binary_whole_overlapping_lesions = _summary_nonzero(mflat[idx_whole] * ones_whole)
                change_map_lesion = _lesion_vs_reference(mflat[inside_bin], mflat[nawm_bin])

                tract_lesion_metrics_map = {
                    "map_tract_pd_lesion_pd": stats_map_tract_pd_lesion_pd,
                    "map_tract_binary_lesion_binary": stats_map_tract_binary_lesion_binary,
                    "map_tract_pd_dilated_lesion": stats_map_tract_pd_dilated_lesion,
                    "map_tract_binary_dilated_lesion": stats_map_tract_binary_dilated_lesion,
                    "map_tract_pd_penumbra": stats_map_tract_pd_penumbra,
                    "map_tract_binary_penumbra": stats_map_tract_binary_penumbra,
                    "map_tract_lesion_vs_nawm": change_map_lesion,
                    "map_tract_binary_whole_overlapping_lesions": stats_map_tract_binary_whole_overlapping_lesions,
                    "map_tract_pd_nawm": _weighted_summary(mflat[tc.idx_nz][nawm_nz], tc.vals[nawm_nz]),
                    "map_tract_binary_nawm": _summary_nonzero(mflat[nawm_bin]),
                    **{f"map_tract_binary_{suffix}": _summary_reduced(mflat[idx]) for suffix, idx in outer_rings.items()},
                }
                for suffix, metrics in tract_lesion_metrics_map.items():
                    for metric_name, metric_value in metrics.items():
                        tract_row[f"{map_name}_{lesion_name}_{suffix}_{metric_name}"] = metric_value

        for map_name, mflat in map_flat.items():
            if map_negative[map_name]:
                logger.warning(f"Warning: Map '{map_name}' contains negative values; check input data.")

            # Whole-tract map stats (tract PD / tract binary / tract skeleton)
            stats_map_tract_pd = _weighted_summary(mflat[tc.idx_nz], tc.vals)
            stats_map_tract_binary = _summary_nonzero(mflat[tc.idx_bin] * tbin_flat[tc.idx_bin])
            stats_map_tract_skeletonised = _summary_nonzero(mflat[tc.idx_skel])     # the skeleton is binary
            tract_array_metrics_map = {
                "map_tract_pd": stats_map_tract_pd,
                "map_tract_binary": stats_map_tract_binary,
                "map_tract_skeletonised": stats_map_tract_skeletonised,
            }
            for suffix, metrics in tract_array_metrics_map.items():
                for metric_name, metric_value in metrics.items():
                    tract_row[f"{map_name}_{suffix}_{metric_name}"] = metric_value

        rows.append(pd.DataFrame([tract_row]))

    return pd.concat([pd.DataFrame()] + rows, ignore_index=True) if rows else pd.DataFrame()


def _compute_whole_brain_lesion_microstructure_metrics_fast(
    subject: str,
    map_arrays: Dict[str, np.ndarray],
    lesion_arrays: Dict[str, np.ndarray],
    binary_lesion_arrays: Dict[str, np.ndarray],
    penumbras_arrays: Dict[str, np.ndarray],
    dilated_lesion_arrays: Dict[str, np.ndarray],
    lesion_contexts: Optional[Dict[str, _LesionContext]] = None,
) -> pd.DataFrame:
    """Whole-brain lesion-restricted microstructure (support-based): each map summarised over the lesion regions."""
    row: Dict[str, Any] = {"subject": subject}
    map_flat = {name: _flat(arr) for name, arr in map_arrays.items()}

    for lesion_name, binary_array in binary_lesion_arrays.items():
        # Optional regions that were not supplied are empty (and are not written).
        empty = np.zeros_like(binary_array)
        lc = (lesion_contexts or {}).get(lesion_name) or _LesionContext(
            lesion_arrays.get(lesion_name, empty), binary_array,
            dilated_lesion_arrays.get(lesion_name, empty), penumbras_arrays.get(lesion_name, empty))

        for map_name, mflat in map_flat.items():
            # (suffix, summary): the lesion probability weights the lesion_pd region; the others are binary regions.
            summaries = [("map_wb_lesion_binary", _summary_nonzero(mflat[lc.idx_bin] * lc.binary[lc.idx_bin]))]
            if lesion_name in lesion_arrays:
                summaries.append(("map_wb_lesion_pd", _weighted_summary(mflat[lc.idx_pd], lc.lesion[lc.idx_pd])))
            if lesion_name in penumbras_arrays:
                summaries.append(("map_wb_penumbra", _summary_nonzero(mflat[lc.idx_pen] * lc.penumbra[lc.idx_pen])))
            if lesion_name in dilated_lesion_arrays:
                summaries.append(("map_wb_dilated_lesion", _summary_nonzero(mflat[lc.idx_dil] * lc.dilated[lc.idx_dil])))

            for suffix, stats in summaries:
                for stat_name, stat_value in stats.items():
                    row[f"{map_name}_{lesion_name}_{suffix}_{stat_name}"] = stat_value

    return pd.DataFrame([row])


def _skeleton_metrics_fast(subject: str, map_arrays: Dict[str, np.ndarray], skeleton_array: np.ndarray) -> pd.DataFrame:
    """map x skeleton statistics over the skeleton's support only."""
    # The skeleton is a binary mask: a voxel is on it or not. Weighting the map by the mask's values (255 for some scikit-image
    # versions, fractions after a linear resampling) would make every statistic, PSMD included, one of map x weight.
    idx = np.flatnonzero(_flat(skeleton_array) > 0.5)
    rows: List[Dict[str, Any]] = []
    for map_name, map_array in map_arrays.items():
        stats_skeleton_map = _summary_nonzero(_flat(map_array)[idx])
        row = {"subject": subject, "map_name": map_name}
        for stat_name, stat_value in stats_skeleton_map.items():
            row[f"skeleton_{stat_name}"] = stat_value
        rows.append(row)
    return pd.DataFrame(rows)


def _remove_small_or_flat_components(binary_array: np.ndarray) -> None:
    """In place: drop 26-connected components with < 4 voxels or flat in one axis.

    Equivalent to the legacy per-component loop (boolean mask + argwhere of the
    whole volume for every component) using one bincount and the components'
    bounding boxes: ptp(coords, axis) == 0 on some axis  <=>  bbox extent == 1.
    """
    labeled, n = ndimage.label(binary_array, structure=STRUCTURE_26)
    if n == 0:
        return
    sizes = np.bincount(labeled.ravel(), minlength=n + 1)
    boxes = ndimage.find_objects(labeled)
    drop = np.zeros(n + 1, dtype=bool)
    for comp_id, box in enumerate(boxes, start=1):
        if sizes[comp_id] < 4 or any(s.stop - s.start == 1 for s in box):
            drop[comp_id] = True
    if drop.any():
        binary_array[drop[labeled]] = 0


def _skeletonize_cropped(mask: np.ndarray) -> np.ndarray:
    """skeletonize() of a binary mask, computed inside its bounding box.

    The 3-D thinning is local and treats everything outside the array as
    background, so a bounding box with a 2-voxel zero margin gives the same
    skeleton as the full volume for a fraction of the cost.
    """
    nz = np.nonzero(mask)
    if len(nz[0]) == 0:
        return skeletonize(mask)
    box = tuple(slice(max(int(a.min()) - 2, 0), min(int(a.max()) + 3, s)) for a, s in zip(nz, mask.shape))
    sub = skeletonize(mask[box])
    # Keep skeletonize()'s own dtype and values: depending on the scikit-image
    # version it is bool (True) or uint8 (255), and downstream code multiplies
    # the skeleton with the maps, so the value 255 vs 1 matters.
    out = np.zeros(mask.shape, dtype=sub.dtype)
    out[box] = sub
    return out


# ============================================================================
# Input Bundle Dataclass
# ============================================================================

@dataclass
class InputBundle:
    """Container for all loaded and processed input arrays."""
    ref_img: Nifti1Image
    brain_array: Optional[np.ndarray]
    lesion_arrays: Dict[str, np.ndarray]
    binary_lesion_arrays: Dict[str, np.ndarray]
    dilated_lesion_arrays: Dict[str, np.ndarray]
    penumbras_arrays: Dict[str, np.ndarray]
    distance_maps: Dict[str, np.ndarray]    # mm to the nearest lesion voxel, per lesion set
    tract_arrays: Dict[str, np.ndarray]
    thresholded_tract_arrays: Dict[str, np.ndarray]
    tract_binary_arrays: Dict[str, np.ndarray]
    tract_skeletons_arrays: Dict[str, np.ndarray]
    map_arrays: Dict[str, np.ndarray]
    skeleton_array: Optional[np.ndarray]
    voxel_vol: float
    tract_atlas_files_dict: Dict[str, str]
    ventricle_array: Optional[np.ndarray]
    # Path to tracts.trk for each tract (None when atlas fallback was used)
    tract_trk_files: Dict[str, Optional[str]]
    # Origin of each tract: "tractography" or "atlas"
    tract_source: Dict[str, str] = field(default_factory=dict)

# ============================================================================
# Input Loading
# ============================================================================

def _load_inputs(args: argparse.Namespace) -> InputBundle:
    """Load, validate, and process all input files into arrays.

    Handles path expansion, array loading, shape validation, and per-array
    transformations (binarisation, dilation, skeletonisation, etc.).
    """
    # Resolve inputs ───────────────────────────────────────────────────────
    t1_files = _expand_path_spec(args.t1_path)
    if not t1_files:
        raise FileNotFoundError(f"No T1 files found: {args.t1_path}")
    logger.info("T1 files (%d): %s ...", len(t1_files), t1_files[0])

    lesion_dict: Dict[str, List[str]] = {}
    if args.lesion_path:
        for lesion_name, lesion_path in _parse_lesion_spec(args.lesion_path):
            lesion_files = _expand_path_spec(lesion_path)
            if lesion_files:
                lesion_dict[lesion_name] = lesion_files
                logger.info("Lesion '%s': %d file(s)", lesion_name, len(lesion_files))

    tract_files: List[str] = []
    if args.tracts_path:
        tract_files = _expand_path_spec(args.tracts_path)
        logger.info("Tract files (%d)", len(tract_files))

    skeleton_files: List[str] = []
    skeleton_array = None
    if args.skeleton_path:
        skeleton_files = _expand_path_spec(args.skeleton_path)
        if skeleton_files:
            skeleton_array = _load_array(skeleton_files[0])
        logger.info("Skeleton files (%d)", len(skeleton_files))

    tract_atlas_files_dict: Dict[str, str] = {}
    if args.tract_atlas_path:
        tract_atlas_files = _expand_path_spec(args.tract_atlas_path)
        logger.info("Tract atlas files (%d)", len(tract_atlas_files))
        for atlas_file in tract_atlas_files:
            atlas_name = str(Path(atlas_file).name).split(".nii.gz", 1)[0]
            tract_atlas_files_dict[atlas_name] = atlas_file

    # Load pre-computed QC report when available; inline QC is used as fallback.
    qc_map: Dict[str, bool] = {}
    if args.qc_report and Path(args.qc_report).is_file():
        qc_map = _load_qc_report(args.qc_report)
    elif args.qc_report:
        logger.warning("QC report not found at %s - falling back to inline QC per tract", args.qc_report)

    map_file_groups: Dict[str, List[str]] = {}
    if args.maps:
        for name, path in _parse_named_spec(args.maps):
            map_file_groups[name] = _expand_path_spec(path)
            logger.info("Map '%s': %d file(s)", name, len(map_file_groups[name]))

    # Validate reference T1 is readable.
    ref_img = _load(t1_files[0])
    logger.info("Reference T1 shape: %s  zooms: %s", ref_img.shape, ref_img.header.get_zooms()[:3])

    ventricle_array: Optional[np.ndarray] = None
    if getattr(args, "parcellation_path", None) and getattr(args, "ventricle_labels", None):
        parcellation_files = _expand_path_spec(args.parcellation_path)
        if parcellation_files:
            ventricle_array = _ventricle_mask(parcellation_files[0], args.ventricle_labels, ref_img)
        else:
            logger.warning("No parcellation found at %s; periventricular / deep lesion volumes are not computed.",
                           args.parcellation_path)

    # Load critical arrays once for reuse across metric families ───────────
    brain_array: Optional[np.ndarray] = None

    lesion_arrays: Dict[str, np.ndarray] = {}
    binary_lesion_arrays: Dict[str, np.ndarray] = {}
    dilated_lesion_arrays: Dict[str, np.ndarray] = {}
    penumbras_arrays: Dict[str, np.ndarray] = {}
    distance_maps: Dict[str, np.ndarray] = {}

    tract_arrays: Dict[str, np.ndarray] = {}
    thresholded_tract_arrays: Dict[str, np.ndarray] = {}
    tract_binary_arrays: Dict[str, np.ndarray] = {}
    tract_skeletons_arrays: Dict[str, np.ndarray] = {}
    tract_trk_files: Dict[str, Optional[str]] = {}

    map_arrays: Dict[str, np.ndarray] = {}
    voxel_vol = _voxel_vol_mm3(ref_img)

    # Brain (for TIV).
    if args.skip_skullstrip:
        brain_array = _load_array(t1_files[0])
    elif args.t1_brain_path:
        brain_array = _load_array(args.t1_brain_path)

    # Lesions (one array per named lesion type).
    for lesion_name, lesion_files in lesion_dict.items():
        lesion_arrs = [_load_array(p) for p in lesion_files]
        if len(lesion_arrs) > 1:
            logger.warning(
                "Lesion '%s': %d files found; only the first (%s) will be used. "
                "If you intended to combine them, provide a single pre-merged file.",
                lesion_name, len(lesion_arrs), lesion_files[0],
            )
        lesion_array = np.where(lesion_arrs[0] < 0.5, 0, lesion_arrs[0])
        binary_array = np.where(lesion_array < 0.5, 0, 1).astype(np.uint8)

        # Remove components that are too small (<4 voxels) or flat in one dimension (2D).
        # Matches the same thresholds applied in lesion_metrics.py so all downstream
        # arrays (lesion_array, dilated, penumbra, edges, tract metrics) see the same data.
        _remove_small_or_flat_components(binary_array)
        lesion_array = lesion_array * binary_array  # zero PD where components were removed

        # Tissue around the lesions, by distance in mm (not by a number of voxels): penumbra = within DILATION_MM of a
        # lesion but not lesion; dilated lesion = lesion + penumbra.
        distance_map = _distance_to_lesion(binary_array, ref_img.header.get_zooms()[:3])
        dilated_array = (distance_map <= DILATION_MM).astype(np.uint8)
        penumbra_array = ((distance_map > 0) & (distance_map <= DILATION_MM)).astype(np.uint8)
        distance_maps[lesion_name] = distance_map

        lesion_arrays[lesion_name] = lesion_array
        binary_lesion_arrays[lesion_name] = binary_array
        dilated_lesion_arrays[lesion_name] = dilated_array
        penumbras_arrays[lesion_name] = penumbra_array

        logger.info("Lesion '%s' array shape: %s", lesion_name, lesion_arrays[lesion_name].shape)

    # Maps (one array per named map).
    for map_name, map_paths in map_file_groups.items():
        map_arrs = [_load_array(p) for p in map_paths]
        if len(map_arrs) > 1:
            logger.warning(
                "Map '%s': %d files found; only the first (%s) will be used. "
                "If you intended to combine them, provide a single pre-merged file.",
                map_name, len(map_arrs), map_paths[0],
            )
        map_arrays[map_name] = map_arrs[0]
        logger.info("Map '%s' array shape: %s", map_name, map_arrays[map_name].shape)

    # Tracts loading logic depending on tract_mode
    tract_source: Dict[str, str] = {}

    if args.tract_mode == "atlas":
        if not tract_atlas_files_dict:
            raise FileNotFoundError(
                f"Tract mode is 'atlas', but no tract atlas files were found "
                f"(searched: '{args.tract_atlas_path}'). FreeSurfer atlas tracts are required for atlas mode."
            )
        for atlas_name, atlas_file in tract_atlas_files_dict.items():
            logger.info("Tract mode: atlas. Using atlas for '%s'.", atlas_name)
            tract_arr = _load_array(atlas_file)
            tract_arrays[atlas_name] = tract_arr
            tract_source[atlas_name] = "atlas"

            atlas_trk = str(Path(atlas_file).with_name(f"{atlas_name}_in_subject.trk"))
            if Path(atlas_trk).exists():
                tract_trk_files[atlas_name] = atlas_trk
            else:
                tract_trk_files[atlas_name] = None
                logger.warning("No atlas .trk file found for '%s'", atlas_name)

            # Binary tract mask. The original zeroed values below 0.0 * max(tract) and then took
            # '> 0'; that is exactly 'tract_arr > 0' (negatives and NaN are False either way), so the
            # thresholded float copy of every tract is no longer built (~2 GB saved over 42 tracts).
            tract_binary_arrays[atlas_name] = tract_arr > 0
            tract_skeletons_arrays[atlas_name] = _skeletonize_cropped(tract_binary_arrays[atlas_name])
    else:
        # Tracts from tractography outputs
        for tract_file in tract_files:
            tract_name = str(Path(tract_file).parent.name).split("_avg16", 1)[0]
            tract_dir_name = str(Path(tract_file).parent.name)
            used_atlas = False

            if args.tract_mode == "tractography":
                logger.info("Tract mode: tractography. Using tractography for '%s' (QC ignored).", tract_name)
                tract_arr = _load_array(tract_file)
                tract_source[tract_name] = "tractography"
            else:  # tractography-atlas (default)
                qc_pass = _qc_pass_for_tract(tract_dir_name, tract_file, qc_map)
                if qc_pass is True:
                    logger.info("Tract QC PASS for '%s' (tractography-atlas mode)", tract_name)
                    tract_arr = _load_array(tract_file)
                    tract_source[tract_name] = "tractography"
                else:
                    logger.warning("Tract QC FAIL for '%s' (tractography-atlas mode)", tract_name)
                    if tract_atlas_files_dict and tract_name in tract_atlas_files_dict:
                        logger.info("Loading tract atlas file for '%s' due to QC fail.", tract_name)
                        tract_arr = _load_array(tract_atlas_files_dict[tract_name])
                        used_atlas = True
                        tract_source[tract_name] = "atlas"
                    else:
                        logger.warning("No atlas file found for '%s'. Skipping tract '%s'", tract_name, tract_name)
                        continue

            tract_arrays[tract_name] = tract_arr

            if not used_atlas:
                trk_candidate = str(Path(tract_file).parent / "tracts.trk")
                tract_trk_files[tract_name] = trk_candidate if Path(trk_candidate).exists() else None
            else:
                atlas_nii = tract_atlas_files_dict.get(tract_name, "")
                atlas_trk = str(Path(atlas_nii).with_name(f"{tract_name}_in_subject.trk"))
                if Path(atlas_trk).exists():
                    tract_trk_files[tract_name] = atlas_trk
                else:
                    tract_trk_files[tract_name] = None
                    logger.warning("No atlas .trk file found for '%s'", tract_name)

            # Binary tract mask. The original zeroed values below 0.0 * max(tract) and then took
            # '> 0'; that is exactly 'tract_arr > 0' (negatives and NaN are False either way), so the
            # thresholded float copy of every tract is no longer built (~2 GB saved over 42 tracts).
            tract_binary_arrays[tract_name] = tract_arr > 0
            tract_skeletons_arrays[tract_name] = _skeletonize_cropped(tract_binary_arrays[tract_name])

        # For tractography-atlas mode, check for expected atlas tracts missing from tractography
        # (e.g. brainstem tracts like CST or MCP omitted by TRACULA due to low inferior DWI FOV coverage).
        # Fall back to atlas so downstream metrics cover the full tractography-atlas set.
        if args.tract_mode == "tractography-atlas" and tract_atlas_files_dict:
            for atlas_name, atlas_file in tract_atlas_files_dict.items():
                if atlas_name not in tract_arrays:
                    logger.warning(
                        "Tract '%s' not present in tractography outputs; "
                        "falling back to tract atlas (tractography-atlas mode).",
                        atlas_name,
                    )
                    tract_arr = _load_array(atlas_file)
                    tract_arrays[atlas_name] = tract_arr
                    tract_source[atlas_name] = "atlas"

                    atlas_trk = str(Path(atlas_file).with_name(f"{atlas_name}_in_subject.trk"))
                    if Path(atlas_trk).exists():
                        tract_trk_files[atlas_name] = atlas_trk
                    else:
                        tract_trk_files[atlas_name] = None
                        logger.warning("No atlas .trk file found for '%s'", atlas_name)

                    # Binary tract mask. The original zeroed values below 0.0 * max(tract) and then took
                    # '> 0'; that is exactly 'tract_arr > 0' (negatives and NaN are False either way), so the
                    # thresholded float copy of every tract is no longer built (~2 GB saved over 42 tracts).
                    tract_binary_arrays[atlas_name] = tract_arr > 0
                    tract_skeletons_arrays[atlas_name] = _skeletonize_cropped(tract_binary_arrays[atlas_name])

    # Validate that all loaded arrays share the same voxel grid shape as the reference T1.
    reference_shape = tuple(ref_img.shape)
    _assert_shape_compatibility(reference_shape, lesion_arrays, "lesion")
    _assert_shape_compatibility(reference_shape, map_arrays, "map")
    _assert_shape_compatibility(reference_shape, tract_arrays, "tract")

    return InputBundle(
        ref_img=ref_img,
        brain_array=brain_array,
        lesion_arrays=lesion_arrays,
        binary_lesion_arrays=binary_lesion_arrays,
        dilated_lesion_arrays=dilated_lesion_arrays,
        penumbras_arrays=penumbras_arrays,
        distance_maps=distance_maps,
        tract_arrays=tract_arrays,
        thresholded_tract_arrays=thresholded_tract_arrays,
        tract_binary_arrays=tract_binary_arrays,
        tract_skeletons_arrays=tract_skeletons_arrays,
        map_arrays=map_arrays,
        skeleton_array=skeleton_array,
        voxel_vol=voxel_vol,
        tract_atlas_files_dict=tract_atlas_files_dict,
        ventricle_array=ventricle_array,
        tract_trk_files=tract_trk_files,
        tract_source=tract_source,
    )

# ============================================================================
# CLI Argument Parsing
# ============================================================================

def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract brain imaging metrics and write a CSV row per subject."
    )
    parser.add_argument("--subject", required=False, default="", help="Subject ID")
    parser.add_argument("--t1-path", required=True, help="T1 file / dir / glob")
    parser.add_argument(
        "--lesion-path", required=False, default=None,
        help="Lesion mask spec: [name=]path,[name2=]path2 (default name: lesion)",
    )
    parser.add_argument(
        "--tracts-path", required=False, default=None,
        help="Tracts file / dir / glob",
    )
    parser.add_argument(
        "--t1-brain-path", required=False, default=None,
        help="Skull-stripped T1 brain image for TIV computation",
    )
    parser.add_argument(
        "--skip-skullstrip", action="store_true", default=False,
        help="Do not attempt TIV computation (no brain file available)",
    )
    parser.add_argument(
        "--maps", required=False, default="",
        help="Map spec: name;path,name2;path2",
    )

    parser.add_argument(
        "--skeleton-path", required=False, default=None,
        help="Skeleton mask file / dir / glob (optional)",
    )

    parser.add_argument(
        "--parcellation-path", required=False, default=None,
        help="Label image (GIF / NeuroMorph or FreeSurfer aparc+aseg) in any space; with --ventricle-labels it gives the "
             "ventricle mask for the periventricular and deep lesion volumes (optional).",
    )
    parser.add_argument(
        "--ventricle-labels", required=False, default=None,
        help="Comma-separated label values of the ventricles in --parcellation-path "
             "(GIF: 50,51,52,53; FreeSurfer: 4,5,43,44).",
    )

    parser.add_argument(
        "--tract-atlas-path", required=False, default=None,
        help="Tract atlas files / dir / glob (optional)",
    )

    parser.add_argument(
        "--output", required=False, default="./",
        help="Output path",
    )

    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose logging")

    parser.add_argument(
        "--tract-mode", required=False, default="tractography-atlas",
        choices=["tractography", "tractography-atlas", "atlas"],
        help="Tract analysis mode: tractography, tractography-atlas (default), or atlas",
    )

    parser.add_argument(
        "--qc-report", required=False, default=None,
        help="Path to tract_qc qc_report.csv. When provided, QC pass/fail is read "
             "from this file instead of being recomputed inline per tract.",
    )

    return parser.parse_args()

# ============================================================================
# CLI
# ============================================================================

def main() -> int:
    args = _parse_args()

    log_level = logging.INFO if args.verbose else logging.WARNING

    # Remove all handlers installed by library imports (we own the only handler).
    root_logger = logging.getLogger()
    for _h in root_logger.handlers[:]:
        root_logger.removeHandler(_h)
    root_logger.setLevel(log_level)
    _stream_handler = logging.StreamHandler(sys.stdout)
    _stream_handler.setLevel(log_level)
    _stream_handler.setFormatter(logging.Formatter("%(levelname)s %(message)s"))
    root_logger.addHandler(_stream_handler)

    # Nuclear fallback: disable INFO/DEBUG globally when not verbose.
    # logging.disable() is a process-wide threshold that no library code can override
    # by calling logging.getLogger().setLevel(logging.INFO) after this point.
    if args.verbose:
        logging.disable(logging.NOTSET)   # reset any prior disable
    else:
        logging.disable(logging.INFO)     # suppress INFO and DEBUG everywhere

    inputs = _load_inputs(args)

    # ── Write metrics to separate CSVs ───────────────────────────────────────
    output_path = Path(args.output)
    metrics_dir = output_path / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)

    # Optional extension hook: only runs if leukoquant/utils/extra/ is
    # present locally (see _load_optional_extension above). No-op otherwise.
    if _extension is not None and inputs.tract_arrays and hasattr(_extension, "on_tracts_ready"):
        _extension.on_tracts_ready(
            subject=args.subject,
            output_dir=str(output_path),
            tract_arrays=inputs.tract_arrays,
            tract_binary_arrays=inputs.tract_binary_arrays,
            tract_trk_files=inputs.tract_trk_files,
            voxel_sizes=inputs.ref_img.header.get_zooms()[:3],
            voxel_vol=inputs.voxel_vol,
        )

    # Per-lesion-set supports / labelling / distance transform, built once and
    # shared by every consumer below (whole-brain, tract-level).
    lesion_contexts: Dict[str, _LesionContext] = {
        _ln: _LesionContext(
            inputs.lesion_arrays[_ln], inputs.binary_lesion_arrays[_ln],
            inputs.dilated_lesion_arrays[_ln], inputs.penumbras_arrays[_ln],
            distance=inputs.distance_maps.get(_ln),
        )
        for _ln in inputs.lesion_arrays
    }

    # ── Metric family: whole-brain ───────────────────────────────────────────
    wb_df = _compute_whole_brain_metrics(
        subject=args.subject,
        brain_array=inputs.brain_array,
        lesion_dict=inputs.lesion_arrays,
        voxel_vol_mm3=inputs.voxel_vol,
    )
 
    if inputs.ventricle_array is not None:
        for _ln, _lesion in inputs.lesion_arrays.items():
            for _k, _v in _periventricular_deep_metrics(
                    _ln, _lesion, inputs.ventricle_array, inputs.ref_img.header.get_zooms()[:3], inputs.voxel_vol).items():
                wb_df[_k] = _v

    # Write whole-brain metrics.
    wb_csv = metrics_dir / "whole_brain_metrics.csv"
    print("✓ Whole Brain Metrics CSV saved")
    write_csv(wb_df, str(wb_csv))

    # ── Metric family: whole-brain lesion microstructure ─────────────────────
    # Deliberately NOT yet declared in metrics_workflow.smk's rule output: --
    # see future_plans/whole_brain_lesion_metrics.md in leukoquant-experiments.
    # Writing it unconditionally here lets every subject whose metrics rule
    # fires from now on pick it up as a side effect of work already scheduled,
    # with no extra Snakemake-triggered reruns; declaring it as a tracked
    # output later is a separate, deliberate step for backfilling subjects
    # that finish before that point.
    wb_lesion_df = _compute_whole_brain_lesion_microstructure_metrics_fast(
        subject=args.subject,
        map_arrays=inputs.map_arrays,
        lesion_arrays=inputs.lesion_arrays,
        binary_lesion_arrays=inputs.binary_lesion_arrays,
        penumbras_arrays=inputs.penumbras_arrays,
        dilated_lesion_arrays=inputs.dilated_lesion_arrays,
        lesion_contexts=lesion_contexts,
    )
    wb_lesion_csv = metrics_dir / "whole_brain_lesion_metrics.csv"
    print("✓ Whole Brain Lesion Microstructure Metrics CSV saved")
    write_csv(wb_lesion_df, str(wb_lesion_csv))

    # ── Peak width of Whole skeleton of tracts ───────────────────────────────────────────
    # Always written (empty DataFrame when no skeleton) so Snakemake can require the file.
    skeleton_metrics_df = pd.DataFrame()
    if inputs.skeleton_array is not None:
        skeleton_metrics_df = _skeleton_metrics_fast(args.subject, inputs.map_arrays, inputs.skeleton_array)

    skeleton_csv = metrics_dir / "skeleton_metrics.csv"
    print("✓ Skeleton Metrics CSV saved")
    write_csv(skeleton_metrics_df, str(skeleton_csv))

    # ── Lesion-level metrics ─────────────────────────────────────────────────
    if inputs.tract_arrays and inputs.map_arrays:
        from leukoquant.utils.lesion_metrics import _compute_lesion_metrics

        # Compute TIV volume (mm3) for lesion metrics
        tiv_volume = _compute_tiv_from_array(inputs.brain_array, inputs.voxel_vol) if inputs.brain_array is not None else 0.0

        # Sorted flat indices of every binary tract, computed once for all lesion sets
        tract_binary_idx = {n: np.flatnonzero(b > 0) for n, b in inputs.tract_binary_arrays.items()}

        # Run lesion-level metrics for each named lesion set and write outputs
        for lesion_name, lesion_array in inputs.lesion_arrays.items():
            lesion_df, tract_df = _compute_lesion_metrics(
                args.subject,
                lesion_array,
                tiv_volume,
                inputs.tract_binary_arrays,
                inputs.map_arrays,
                inputs.ref_img.header.get_zooms()[:3],
                tract_binary_idx=tract_binary_idx,
            )

            # Drop rows in tract_df where all numeric columns are NaN
            if not tract_df.empty:
                numeric_cols = tract_df.select_dtypes(include=[np.number]).columns.tolist()
                if numeric_cols:
                    tract_df = tract_df.dropna(how="all", subset=numeric_cols)

            lesion_csv = metrics_dir / f"{lesion_name}_lesion_metrics.csv"
            print(f"✓ {lesion_name.upper()} Lesion Metrics CSV saved")
            write_csv(lesion_df, str(lesion_csv))

            tract_csv = metrics_dir / f"{lesion_name}_tract_aggregated_lesion_metrics.csv"
            print(f"✓ {lesion_name.upper()} Tract Aggregated Lesion Metrics CSV saved")
            write_csv(tract_df, str(tract_csv))

            # Optional extension hook (see _load_optional_extension above).
            # Label with 26-connectivity (STRUCTURE_26) so lesion_label 1..n_components matches lesion_df.
            if _extension is not None and hasattr(_extension, "on_lesion_ready"):
                labeled_lesions, n_components = ndimage.label(
                    inputs.binary_lesion_arrays[lesion_name], structure=STRUCTURE_26
                )
                voxel_sizes = inputs.ref_img.header.get_zooms()[:3]
                _extension.on_lesion_ready(
                    subject=args.subject,
                    output_dir=str(output_path),
                    lesion_name=lesion_name,
                    labeled_lesions=labeled_lesions,
                    n_components=n_components,
                    lesion_df=lesion_df,
                    tract_arrays=inputs.tract_arrays,
                    tract_binary_arrays=inputs.tract_binary_arrays,
                    tract_trk_files=inputs.tract_trk_files,
                    affine=inputs.ref_img.affine,
                    voxel_sizes=voxel_sizes,
                    voxel_vol=inputs.voxel_vol,
                )

    # ── Tract-level metrics ──────────────────────────────────────────────────
    tract_level_metrics_df = _compute_tract_level_metrics_fast(
        subject=args.subject,
        tract_arrays=inputs.tract_arrays,
        tract_binary_arrays=inputs.tract_binary_arrays,
        tract_skeletons_arrays=inputs.tract_skeletons_arrays,
        tract_source=inputs.tract_source,
        lesion_arrays=inputs.lesion_arrays,
        binary_lesion_arrays=inputs.binary_lesion_arrays,
        dilated_lesion_arrays=inputs.dilated_lesion_arrays,
        penumbras_arrays=inputs.penumbras_arrays,
        map_arrays=inputs.map_arrays,
        voxel_vol=inputs.voxel_vol,
        voxel_sizes=inputs.ref_img.header.get_zooms()[:3],
        lesion_contexts=lesion_contexts,
    )

    # Write tract-level metrics.
    tract_csv = metrics_dir / "tract_level_metrics.csv"
    print("✓ Tract-Level Metrics CSV saved")
    write_csv(tract_level_metrics_df, str(tract_csv))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
