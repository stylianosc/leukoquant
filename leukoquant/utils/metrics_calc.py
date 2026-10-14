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
    return np.ascontiguousarray(np.asarray(_load(path).dataobj))

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


def _summary_nonzero(values: np.ndarray) -> Dict[str, float]:
    """Return summary stats over non-zero values (zeros if empty).

    Keys: min, max, mean, sum, median, std, 25th, 75th, IQR, peak_width (95th - 5th percentile),
    skew and kurtosis (excess kurtosis; see _shape_stats).
    """
    non_zero = values[np.nonzero(values)]
    if non_zero.size == 0:
        return {
            "min": 0.0,
            "max": 0.0,
            "mean": 0.0,
            "sum": 0.0,
            "median": 0.0,
            "std": 0.0,
            "25th": 0.0,
            "75th": 0.0,
            "IQR": 0.0,
            "peak_width": 0.0,
            "skew": 0.0,
            "kurtosis": 0.0,
        }

    # One np.percentile call (a single partition) instead of four; the values
    # are bit-identical to four separate calls (same linear interpolation).
    _p5, _p25, _p75, _p95 = np.percentile(non_zero, [5, 25, 75, 95])
    p25, p75, p95, p5 = float(_p25), float(_p75), float(_p95), float(_p5)
    # Float64 for the moments, so float32 maps do not lose precision in the 3rd and 4th powers.
    non_zero64 = non_zero.astype(np.float64, copy=False)
    mean = float(np.mean(non_zero64))
    skew, kurtosis = _shape_stats(non_zero64, mean)
    return {
        "min": float(np.min(non_zero)),
        "max": float(np.max(non_zero)),
        "mean": float(np.mean(non_zero)),
        "sum": float(np.sum(non_zero)),
        "median": float(np.median(non_zero)),
        "std": float(np.std(non_zero)),
        "25th": p25,
        "75th": p75,
        "IQR": p75 - p25,
        "peak_width": p95 - p5,
        "skew": skew,
        "kurtosis": kurtosis,
    }

def _summary_change_nonzero(lesion_map: np.ndarray, values_b: np.ndarray) -> Dict[str, float]:
    """Compute normalised symmetry-ratio statistics comparing inside-lesion vs. outside-lesion values.

    Parameters
    ----------
    lesion_map : np.ndarray
        Binary (or soft) lesion mask. Thresholded at 0.5 to produce a hard binary mask.
    values_b : np.ndarray
        Continuous value array (e.g. tract PD, FA map) with the same shape as ``lesion_map``.

    Returns
    -------
    Dict[str, float]
        Two keys per summary statistic S (min, max, mean, sum, median, std, 25th, 75th, IQR,
        peak_width; not skew or kurtosis, see SHAPE_STATS):
        - ``{S}_change`` : (inside - outside) / (inside + outside). Range [-1, 1].
            +1 means all values are inside the lesion; -1 means all values are outside.
            **Edge cases:**
            - If the lesion has zero overlap with the value array (inside=0, outside>0),
              returns NaN (metric undefined).
            - If both inside and outside are zero (denom=0), returns 0.0.
        - ``{S}_change_absolute`` : abs(inside - outside) / (inside + outside).
            Always >= 0. Returns NaN when inside=0 and denom!=0. Returns 0.0 when denom=0.
    """
    lesion_map = np.where(lesion_map < 0.5, 0, 1)
    non_a = lesion_map ^ 1
    overlap = lesion_map * values_b
    non_overlap = non_a * values_b

    overlap_stats = _summary_nonzero(overlap)
    non_overlap_stats = _summary_nonzero(non_overlap)

    stats: Dict[str, float] = {}
    for operation in overlap_stats.keys():
        if operation in SHAPE_STATS:
            # (a - b) / (a + b) needs two non-negative quantities; skewness and excess kurtosis can be
            # negative, so the ratio is undefined (and unbounded where a + b crosses zero).
            continue
        a_val = overlap_stats[operation]
        b_val = non_overlap_stats[operation]
        denom = a_val + b_val

        if denom == 0:
            change = 0.0
            change_abs = 0.0
        elif a_val == 0:
            # No values inside the lesion - metric is undefined, not -1
            change = float('nan')
            change_abs = float('nan')
        else:
            change = (a_val - b_val) / denom
            change_abs = abs(a_val - b_val) / denom

        stats[f"{operation}_change"] = change
        stats[f"{operation}_change_absolute"] = change_abs

    return stats

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
    {lesion_name}_volume_tiv_normalised_mm3
                       Lesion volume normalised by TIV (mm³) in whole brain (per lesion type).
    {lesion_name}_load_pct
                       Lesion load % in whole brain (per lesion type).
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
            row[f"{lesion_name}_volume_tiv_normalised_mm3"] = lesion_vol / tiv_mm3
        else:
            row[f"{lesion_name}_load_pct"] = None
            row[f"{lesion_name}_volume_tiv_normalised_mm3"] = None


    return pd.DataFrame([row])

# ============================================================================
# Metric family: whole-brain lesion microstructure
# ============================================================================

def _compute_whole_brain_lesion_microstructure_metrics(
    subject: str,
    map_arrays: Dict[str, np.ndarray],
    lesion_arrays: Dict[str, np.ndarray],
    binary_lesion_arrays: Dict[str, np.ndarray],
    penumbras_arrays: Dict[str, np.ndarray],
    dilated_lesion_arrays: Dict[str, np.ndarray],
) -> pd.DataFrame:
    """Compute whole-brain lesion-restricted microstructure metrics.

    Same idea as the tract-level "lesion microstructure" / "penumbra
    microstructure" families the tract loop below computes, just without the
    tract-mask restriction -- every map is summarised over the WHOLE
    lesion/penumbra/dilated-lesion region rather than its intersection with
    one of the 42 tracts. This is the metric family a whole-brain-scope score
    (e.g. Fazekas, a holistic visual rating) is actually comparable to; a
    tract-restricted metric never is, since Fazekas has no per-tract
    structure of its own.

    Runs once per subject (no tract loop), reusing the exact lesion/map
    arrays _load_inputs() already loaded for the tract-level pass below --
    no new data loading, no new CLI flags.

    Columns
    -------
    subject
    {map}_{lesion}_map_wb_lesion_binary_{stat}
        Map values inside the binary lesion mask (whole brain).
    {map}_{lesion}_map_wb_lesion_pd_{stat}
        Map values weighted by lesion probability density (whole brain).
    {map}_{lesion}_map_wb_penumbra_{stat}
        Map values in the penumbra ring around the lesion (whole brain).
    {map}_{lesion}_map_wb_dilated_lesion_{stat}
        Map values in the dilated-lesion (lesion + penumbra) region (whole brain).
    {stat} is one of the _summary_nonzero() keys: min, max, mean, sum,
    median, std, 25th, 75th, IQR, peak_width, skew, kurtosis.
    """
    row: Dict[str, Any] = {"subject": subject}

    for lesion_name, binary_array in binary_lesion_arrays.items():
        lesion_pd_array = lesion_arrays.get(lesion_name)
        penumbra_array = penumbras_arrays.get(lesion_name)
        dilated_array = dilated_lesion_arrays.get(lesion_name)

        for map_name, map_array in map_arrays.items():
            regions: Dict[str, np.ndarray] = {
                "map_wb_lesion_binary": map_array * binary_array,
            }
            if lesion_pd_array is not None:
                regions["map_wb_lesion_pd"] = map_array * lesion_pd_array
            if penumbra_array is not None:
                regions["map_wb_penumbra"] = map_array * penumbra_array
            if dilated_array is not None:
                regions["map_wb_dilated_lesion"] = map_array * dilated_array

            for suffix, region_values in regions.items():
                stats = _summary_nonzero(region_values)
                for stat_name, stat_value in stats.items():
                    row[f"{map_name}_{lesion_name}_{suffix}_{stat_name}"] = stat_value

    return pd.DataFrame([row])

# ============================================================================
# Metric family: tract-lesion load
# ============================================================================

def _compute_lesion_tract_intersection_details(
    lesion_array: np.ndarray,
    tract_array: np.ndarray,
    voxel_vol_mm3: float,
) -> Dict[str, Any]:
    """Compute detailed lesion-tract intersection metrics using connected components.

    For each individual lesion component, compute overlap with the tract mask.
    Aggregate across all intersecting components.

    Parameters
    ----------
    lesion_array : np.ndarray
        Lesion probability/intensity map.
    tract_array : np.ndarray
        Tract binary or probability map.
    voxel_vol_mm3 : float
        Single voxel volume in mm³.

    Returns
    -------
    Dict[str, Any]
        Dictionary with keys:
        - overall_percentage_overlap: 
        Percentage of volume of lesions intersecting tract to total volume of same
        lesions intersecting tract.
        - intersecting_lesions_total_volume: 
        sum of volumes of components intersecting tract
        - {mean/std/min/max/median/25th/75th/IQR}_percentage_tract_lesion_overlap: 
        overlap % stats across components
        - number_lesions_tract_overlap: count of components with any overlap
        - {mean/std/min/max/median/25th/75th/IQR}_lesions_overlap_volume: 
        stats on the overlap volume (mm³) of lesion within a tract.
    """
    # Create/reuse binary lesion and tract masks.
    lesion_binary = np.where(
            lesion_array < 0.5, 0, 1) if lesion_array.dtype in [
                np.float32, np.float64] else lesion_array > 0
    tract_binary = tract_array > 0

    # Identify/reuse connected components in lesion mask.
    labeled_lesions, num_components = ndimage.label(lesion_binary)  # type: ignore

    if num_components > 0:
        component_voxels = ndimage.sum(
            input=lesion_binary,
            labels=labeled_lesions,
            index=np.arange(1, num_components + 1),
        )
        component_volumes_mm3 = np.asarray(
            component_voxels, dtype=float) * voxel_vol_mm3
    else:
        component_volumes_mm3 = np.array([], dtype=float)

    result: Dict[str, Any] = {
        "in_tract_mm3": 0.0,
        "load_in_tract_pct": 0.0,
        "overall_percentage_overlap": 0.0,
        "intersecting_lesions_total_volume": 0.0,
        "number_lesions_tract_overlap": 0,
    }

    tract_volume_mm3 = float(np.sum(tract_binary)) * voxel_vol_mm3
    tract_lesion_load_mm3 = _compute_load(lesion_array, tract_array, voxel_vol_mm3)
    result["in_tract_mm3"] = tract_lesion_load_mm3
    if tract_volume_mm3 > 0:
        result["load_in_tract_pct"] = tract_lesion_load_mm3 / tract_volume_mm3 * 100.0

    lesion_sizes: List[float] = []
    overlap_percentages: List[float] = []
    overlap_volumes: List[float] = []  # Volume of overlap per lesion
    intersecting_total_volume = 0.0
    intersecting_components = 0

    # Process each connected component.
    for component_id in range(1, num_components + 1):
        component_mask = labeled_lesions == component_id
        component_volume = float(component_volumes_mm3[component_id - 1])

        # Compute overlap with tract.
        overlap_mask = np.logical_and(component_mask, tract_binary)
        overlap_volume_count = int(np.sum(overlap_mask)) * voxel_vol_mm3

        if overlap_volume_count > 0:
            lesion_sizes.append(component_volume)

            overlap_percentages.append(
                float(overlap_volume_count) / component_volume * 100.0)
            overlap_volumes.append(overlap_volume_count)
            intersecting_total_volume += component_volume
            intersecting_components += 1

    # Convert to arrays and compute statistics.
    lesion_sizes_array = np.array(lesion_sizes)
    overlap_percentages_array = np.array(overlap_percentages)
    overlap_volumes_array = np.array(overlap_volumes)

    # Total overlap volume tract-lesion (mm³).
    total_overlap_voxels = int(np.sum(np.logical_and(lesion_binary, tract_binary)))
    total_overlap_volume = total_overlap_voxels * voxel_vol_mm3
    if total_overlap_volume != overlap_volumes_array.sum():
        logger.warning(
            "Total overlap volume (%.2f mm³) does not match sum of component overlaps (%.2f mm³).",
            total_overlap_volume, overlap_volumes_array.sum()
        )

    # Update result with computed values.
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
        stats = _summary_nonzero(values)
        for stat_name, stat_value in stats.items():
            result[f"{stat_name}_{suffix}"] = stat_value

    return result

def _compute_tract_distance_metrics(
    tract_array: np.ndarray,
    lesion_array: np.ndarray,
    voxel_sizes: Tuple[float, float, float],
) -> Dict[str, float]:
    """Compute tract-to-lesion distance metrics for one tract array and one lesion array."""
    if tract_array.shape != lesion_array.shape:
        logger.warning(
            "Shapes mismatch: lesion %s vs tract %s; skipping tract distance.",
            lesion_array.shape,
            tract_array.shape,
        )
        return {
            "tract_nearest_lesion_distance": float("nan"),
            "tract_centre_of_mass_nearest_lesion_distance": float("nan"),
        }

    tract_binary = tract_array > 0

    tract_center_of_mass: Optional[Tuple[int, int, int]] = None
    if np.any(tract_binary):
        if tract_binary.ndim != 3:
            logger.warning("tract_array is not 3D (shape=%s); skipping centre-of-mass.", tract_binary.shape)
        else:
            com = ndimage.center_of_mass(tract_binary)
            tract_center_of_mass = (int(com[0]), int(com[1]), int(com[2]))

    lesion_binary = np.where(lesion_array < 0.5, 0, 1)

    overlap_voxels = int(np.sum(np.logical_and(tract_binary, lesion_binary > 0)))
    min_distance = 0.0 if overlap_voxels > 0 else np.inf
    min_distance_com = np.inf

    lesion_inverse = np.where(lesion_binary > 0, 0, 1)
    distance_map = cast(
        np.ndarray,
        ndimage.distance_transform_edt(lesion_inverse, sampling=voxel_sizes),
    )

    if np.any(tract_binary):
        masked_distance = np.where(tract_binary, distance_map, np.inf)
        min_distance = min(min_distance, float(np.min(masked_distance)))

    if tract_center_of_mass is not None:
        min_distance_com = min(min_distance_com, float(distance_map[tract_center_of_mass]))

    return {
        "tract_nearest_lesion_distance": min_distance if np.isfinite(min_distance) else np.nan,
        "tract_centre_of_mass_nearest_lesion_distance": (
            min_distance_com if np.isfinite(min_distance_com) else np.nan
        ),
    }

# ============================================================================
# FAST PATH: sparse, support-based implementations of the per-tract /
# per-lesion statistics.
#
# Why this is exactly equivalent to the dense code it replaces
# ------------------------------------------------------------
# Every statistic below is _summary_nonzero(product of a few dense arrays),
# i.e. it only ever looks at the NON-ZERO entries of an elementwise product,
# in C order. A product can only be non-zero where each factor is non-zero, so
# it is enough to evaluate the *same elementwise expression* on the voxels
# where the sparsest factor is non-zero (a tract's support, a lesion's support
# ...), taken in ascending C-order flat index. Because
#   * the same expression, on the same dtypes (we gather from the original
#     arrays, we never convert them), is applied to every voxel evaluated,
#   * zero products are dropped by _summary_nonzero exactly as before, and
#   * the surviving values reach np.sum/np.mean/np.percentile in the same
#     order and with the same dtype,
# the results are bit-identical to the dense implementation.
#
# The one situation where the argument breaks is a NaN / Inf anywhere in a map
# (NaN * 0 = NaN is "non-zero" in the dense code but outside every support in
# the sparse one). _inputs_are_finite() guards the fast path; if it ever
# returns False the legacy dense implementation is used instead.
# ============================================================================

def _flat(arr: np.ndarray) -> np.ndarray:
    """C-order 1-D view of an array (a copy only if it is not C-contiguous)."""
    return arr.reshape(-1)


def _inputs_are_finite(arrays) -> bool:
    """True if no floating-point array contains NaN/Inf (see the note above)."""
    for arr in arrays:
        if arr is not None and arr.dtype.kind in "fc" and not np.isfinite(arr).all():
            return False
    return True


class _LesionContext:
    """Per-lesion-set data computed once and reused for every tract and map."""

    def __init__(self, lesion_array, binary, dilated, penumbra):
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
        # 6-connected components of the binary lesion (the tract loop's labelling)
        labels, self.n_comp = ndimage.label(binary)
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
            dl, self.det_n = ndimage.label(det_binary.reshape(self.shape))
            self.det_labels = _flat(dl)
        self.det_comp_voxels = np.bincount(self.det_labels, minlength=self.det_n + 1)[1:]
        self._dist = None

    def distance_flat(self, voxel_sizes):
        """Euclidean distance (mm) to the nearest lesion voxel, once per lesion set."""
        if self._dist is None:
            lesion_inverse = np.where(self.lm, 0, 1).reshape(self.shape)
            self._dist = _flat(cast(
                np.ndarray,
                ndimage.distance_transform_edt(lesion_inverse, sampling=voxel_sizes),
            ))
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
        self.idx_skel = np.flatnonzero(_flat(tract_skeleton))
        self.has_negative = bool((self.vals < 0).any())
        # lesion-independent products with each map, filled lazily
        self.map_pd = {}

    def centre_of_mass_voxel(self):
        """int() of scipy's center_of_mass(tract_binary), computed exactly.

        scipy sums integer-valued coordinates in float64 (exact below 2**53)
        and divides by the voxel count, so an integer sum followed by one
        division gives the identical float, and hence the identical int().
        """
        n = len(self.idx_bin)
        coords = np.unravel_index(self.idx_bin, self.shape)
        return tuple(int(int(c.sum(dtype=np.int64)) / n) for c in coords)


def _change_from_stats(overlap_stats, non_overlap_stats):
    """Second half of _summary_change_nonzero (unchanged logic)."""
    stats: Dict[str, float] = {}
    for operation in overlap_stats.keys():
        if operation in SHAPE_STATS:
            # (a - b) / (a + b) needs two non-negative quantities; skewness and excess kurtosis can be
            # negative, so the ratio is undefined (and unbounded where a + b crosses zero).
            continue
        a_val = overlap_stats[operation]
        b_val = non_overlap_stats[operation]
        denom = a_val + b_val

        if denom == 0:
            change = 0.0
            change_abs = 0.0
        elif a_val == 0:
            change = float('nan')
            change_abs = float('nan')
        else:
            change = (a_val - b_val) / denom
            change_abs = abs(a_val - b_val) / denom

        stats[f"{operation}_change"] = change
        stats[f"{operation}_change_absolute"] = change_abs
    return stats


def _fast_summary_change(lc: _LesionContext, tc: _TractContext, values_at_support: np.ndarray):
    """_summary_change_nonzero(lesion_map=lesion, values_b=v) for a v whose
    non-zero voxels all lie in the tract's non-zero support.

    values_at_support: v evaluated at tc.idx_nz.
    """
    inside = lc.lm[tc.idx_nz]
    n_in = int(inside.sum())
    ones_in = np.ones(n_in, dtype=np.int64)                    # lesion_map      (int64)
    ones_out = np.ones(len(inside) - n_in, dtype=np.int64)     # lesion_map ^ 1  (int64)
    overlap_stats = _summary_nonzero(ones_in * values_at_support[inside])
    non_overlap_stats = _summary_nonzero(ones_out * values_at_support[~inside])
    return _change_from_stats(overlap_stats, non_overlap_stats)


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
        stats = _summary_nonzero(values)
        for stat_name, stat_value in stats.items():
            result[f"{stat_name}_{suffix}"] = stat_value

    return result


def _fast_tract_distance_metrics(lc: _LesionContext, tc: _TractContext, voxel_sizes):
    """Fast equivalent of _compute_tract_distance_metrics.

    The Euclidean distance transform depends only on the lesion mask, so it is
    computed once per lesion set (lc.distance_flat) instead of once per tract.
    """
    distance = lc.distance_flat(voxel_sizes)
    has_tract = len(tc.idx_bin) > 0

    overlap_voxels = int(np.count_nonzero(lc.lm[tc.idx_bin]))
    min_distance = 0.0 if overlap_voxels > 0 else np.inf
    min_distance_com = np.inf

    if has_tract:
        min_distance = min(min_distance, float(np.min(distance[tc.idx_bin])))
        com = tc.centre_of_mass_voxel()
        min_distance_com = min(min_distance_com, float(distance[np.ravel_multi_index(com, tc.shape)]))

    return {
        "tract_nearest_lesion_distance": min_distance if np.isfinite(min_distance) else np.nan,
        "tract_centre_of_mass_nearest_lesion_distance": (
            min_distance_com if np.isfinite(min_distance_com) else np.nan
        ),
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

        # map x tract PD over the tract's support (also reused by the whole-tract stats)
        map_tract_pd_at_support = {
            map_name: mflat[tc.idx_nz] * tc.vals for map_name, mflat in map_flat.items()
        }

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
            change_tract_pd_lesion = _fast_summary_change(lc, tc, tc.vals)

            tract_lesion_metrics = {
                "tract_pd_lesion_pd": stats_tract_pd_lesion_pd,
                "tract_pd_lesion_binary": stats_tract_pd_lesion_binary,
                "tract_pd_lesion_dilated": stats_tract_pd_lesion_dilated,
                "tract_pd_lesion_penumbra": stats_tract_pd_lesion_penumbra,
                "tract_pd_change_lesion_binary": change_tract_pd_lesion,
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
                stats_map_tract_pd_lesion_pd = _summary_nonzero(mflat[lc.idx_pd] * tpl_pd)
                stats_map_tract_binary_lesion_binary = _summary_nonzero(mflat[lc.idx_bin] * tbl_bin)
                stats_map_tract_pd_dilated_lesion = _summary_nonzero(mflat[lc.idx_dil] * tpl_dil)
                stats_map_tract_binary_dilated_lesion = _summary_nonzero(mflat[lc.idx_dil] * tbl_dil)
                stats_map_tract_pd_penumbra = _summary_nonzero(mflat[lc.idx_pen] * tpl_pen)
                stats_map_tract_binary_penumbra = _summary_nonzero(mflat[lc.idx_pen] * tbl_pen)
                stats_map_tract_binary_whole_overlapping_lesions = _summary_nonzero(mflat[idx_whole] * ones_whole)
                change_map_tract_pd_lesion = _fast_summary_change(lc, tc, map_tract_pd_at_support[map_name])

                tract_lesion_metrics_map = {
                    "map_tract_pd_lesion_pd": stats_map_tract_pd_lesion_pd,
                    "map_tract_binary_lesion_binary": stats_map_tract_binary_lesion_binary,
                    "map_tract_pd_dilated_lesion": stats_map_tract_pd_dilated_lesion,
                    "map_tract_binary_dilated_lesion": stats_map_tract_binary_dilated_lesion,
                    "map_tract_pd_penumbra": stats_map_tract_pd_penumbra,
                    "map_tract_binary_penumbra": stats_map_tract_binary_penumbra,
                    "map_tract_pd_change_lesion_binary": change_map_tract_pd_lesion,
                    "map_tract_binary_whole_overlapping_lesions": stats_map_tract_binary_whole_overlapping_lesions,
                }
                for suffix, metrics in tract_lesion_metrics_map.items():
                    for metric_name, metric_value in metrics.items():
                        tract_row[f"{map_name}_{lesion_name}_{suffix}_{metric_name}"] = metric_value

        for map_name, mflat in map_flat.items():
            if map_negative[map_name]:
                logger.warning(f"Warning: Map '{map_name}' contains negative values; check input data.")

            # Whole-tract map stats (tract PD / tract binary / tract skeleton)
            stats_map_tract_pd = _summary_nonzero(map_tract_pd_at_support[map_name])
            stats_map_tract_binary = _summary_nonzero(mflat[tc.idx_bin] * tbin_flat[tc.idx_bin])
            stats_map_tract_skeletonised = _summary_nonzero(
                _flat(tract_skeletons_arrays[tract_name])[tc.idx_skel] * mflat[tc.idx_skel]
            )
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
    """Support-based equivalent of _compute_whole_brain_lesion_microstructure_metrics."""
    row: Dict[str, Any] = {"subject": subject}
    map_flat = {name: _flat(arr) for name, arr in map_arrays.items()}

    for lesion_name, binary_array in binary_lesion_arrays.items():
        lc = (lesion_contexts or {}).get(lesion_name) or _LesionContext(
            lesion_arrays[lesion_name], binary_array,
            dilated_lesion_arrays[lesion_name], penumbras_arrays[lesion_name])

        for map_name, mflat in map_flat.items():
            regions: Dict[str, np.ndarray] = {
                "map_wb_lesion_binary": mflat[lc.idx_bin] * lc.binary[lc.idx_bin],
            }
            if lesion_name in lesion_arrays:
                regions["map_wb_lesion_pd"] = mflat[lc.idx_pd] * lc.lesion[lc.idx_pd]
            if lesion_name in penumbras_arrays:
                regions["map_wb_penumbra"] = mflat[lc.idx_pen] * lc.penumbra[lc.idx_pen]
            if lesion_name in dilated_lesion_arrays:
                regions["map_wb_dilated_lesion"] = mflat[lc.idx_dil] * lc.dilated[lc.idx_dil]

            for suffix, region_values in regions.items():
                stats = _summary_nonzero(region_values)
                for stat_name, stat_value in stats.items():
                    row[f"{map_name}_{lesion_name}_{suffix}_{stat_name}"] = stat_value

    return pd.DataFrame([row])


def _skeleton_metrics_fast(subject: str, map_arrays: Dict[str, np.ndarray], skeleton_array: np.ndarray) -> pd.DataFrame:
    """map x skeleton statistics over the skeleton's support only."""
    sk = _flat(skeleton_array)
    idx = np.flatnonzero(sk)
    sk_at = sk[idx]
    rows: List[Dict[str, Any]] = []
    for map_name, map_array in map_arrays.items():
        stats_skeleton_map = _summary_nonzero(_flat(map_array)[idx] * sk_at)
        row = {"subject": subject, "map_name": map_name}
        for stat_name, stat_value in stats_skeleton_map.items():
            row[f"skeleton_{stat_name}"] = stat_value
        rows.append(row)
    return pd.DataFrame(rows)


def _remove_small_or_flat_components(binary_array: np.ndarray) -> None:
    """In place: drop 6-connected components with < 4 voxels or flat in one axis.

    Equivalent to the legacy per-component loop (boolean mask + argwhere of the
    whole volume for every component) using one bincount and the components'
    bounding boxes: ptp(coords, axis) == 0 on some axis  <=>  bbox extent == 1.
    """
    labeled, n = ndimage.label(binary_array)
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


def _dilate_cropped(mask: np.ndarray, iterations: int) -> np.ndarray:
    """ndimage.binary_dilation(mask, iterations) inside the mask's bounding box.

    Dilation is local and the default border value is 0, so a bounding box with
    an (iterations + 1)-voxel margin gives the identical result for a fraction
    of the cost of a full-volume pass.
    """
    out = np.zeros(mask.shape, dtype=bool)
    nz = np.nonzero(mask)
    if len(nz[0]) == 0:
        return out
    m = iterations + 1
    box = tuple(slice(max(int(a.min()) - m, 0), min(int(a.max()) + m + 1, s)) for a, s in zip(nz, mask.shape))
    out[box] = ndimage.binary_dilation(input=mask[box], iterations=iterations)
    return out


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
    tract_arrays: Dict[str, np.ndarray]
    thresholded_tract_arrays: Dict[str, np.ndarray]
    tract_binary_arrays: Dict[str, np.ndarray]
    tract_skeletons_arrays: Dict[str, np.ndarray]
    map_arrays: Dict[str, np.ndarray]
    skeleton_array: Optional[np.ndarray]
    voxel_vol: float
    tract_atlas_files_dict: Dict[str, str]
    # Path to tracts.trk for each tract (None when atlas fallback was used)
    tract_trk_files: Dict[str, Optional[str]]
    # Origin of each tract: "tractography" or "atlas"
    tract_source: Dict[str, str] = field(default_factory=dict)

def _compute_tract_level_metrics_legacy(args, inputs: "InputBundle") -> pd.DataFrame:
    """Dense reference implementation of the tract-level metrics loop.

    Verbatim copy of the original loop from main(); used as the fallback when
    the fast sparse path cannot be proven equivalent (NaN/Inf in an input).
    """
    tract_level_metrics_df = pd.DataFrame()

    for tract_name, tract_array in inputs.tract_arrays.items():
        source = inputs.tract_source.get(tract_name, "tractography")
        tract_row = {
            "subject": args.subject,
            "tract_name": tract_name,
            "tract_source": source,
            "used_atlas": source == "atlas",
            "tract_volume_mm3": float(np.sum(tract_array > 0)) * inputs.voxel_vol,
        }

        # Sanity check if there are negative values in the tract array, which should not be the case for probability maps.
        if np.any(tract_array < 0):
            logger.warning(f"Warning: Tract '{tract_name}' contains negative values; check input data.")

        for lesion_name, lesion_array in inputs.lesion_arrays.items():
            # Sanity check if there are negative values in the lesion array, which should not be the case for probability maps.
            if np.any(lesion_array < 0):
                logger.warning(f"Warning: Lesion '{lesion_name}' contains negative values; check input data.")

            # Compute lesion load details for tract-lesion.
            tract_lesion_load_details = _compute_lesion_tract_intersection_details(
                lesion_array=lesion_array,
                tract_array=tract_array,
                voxel_vol_mm3=inputs.voxel_vol,
            )

            # Compute tract-lesion distance metrics.
            tract_lesion_distance_metrics = _compute_tract_distance_metrics(
                tract_array=tract_array,
                lesion_array=lesion_array,
                voxel_sizes=inputs.ref_img.header.get_zooms()[:3],
            )

            tract_pd_lesion_pd = tract_array * lesion_array
            tract_binary_lesion_binary = inputs.tract_binary_arrays[tract_name] * inputs.binary_lesion_arrays[lesion_name]
            tract_pd_lesion_binary = tract_array * inputs.binary_lesion_arrays[lesion_name]
            tract_pd_lesion_dilated = tract_array * inputs.dilated_lesion_arrays[lesion_name]
            tract_binary_lesion_dilated = inputs.tract_binary_arrays[tract_name] * inputs.dilated_lesion_arrays[lesion_name]
            tract_pd_lesion_penumbra = tract_array * inputs.penumbras_arrays[lesion_name]
            tract_binary_penumbra = inputs.tract_binary_arrays[tract_name] * inputs.penumbras_arrays[lesion_name]

            # Get lesions that overlap with the tract (including lesion areas outside the tract)
            # Connected components on lesion mask
            lesions_labels, _ = cast(
                Tuple[np.ndarray, int],
                ndimage.label(inputs.binary_lesion_arrays[lesion_name]),
            )
            # IDs of lesion components that overlap this tract
            overlap_component_ids = np.unique(lesions_labels * tract_binary_lesion_binary)
            overlap_component_ids = overlap_component_ids[overlap_component_ids != 0]
        
            whole_overlapping_lesions_binary = np.isin(lesions_labels, overlap_component_ids).astype(np.uint8)

            # Get summary stats for tract-lesion intersection.
            # Importance in tract (Tract PD) x Importance in lesion (Lesion PD)
            # Tract PD x Lesion PD
            stats_tract_pd_lesion_pd = _summary_nonzero(tract_pd_lesion_pd)
            # Tract PD x Lesion binary
            stats_tract_pd_lesion_binary = _summary_nonzero(tract_pd_lesion_binary)
            # Tract PD x Dilated lesion
            stats_tract_pd_lesion_dilated = _summary_nonzero(tract_pd_lesion_dilated)
            # Tract PD x Penumbra
            stats_tract_pd_lesion_penumbra = _summary_nonzero(tract_pd_lesion_penumbra)

            # Change: Tract PD x Lesion binary
            change_tract_pd_lesion = _summary_change_nonzero(
                lesion_map=lesion_array, values_b=tract_array)
            
            tract_lesion_metrics = {
                "tract_pd_lesion_pd": stats_tract_pd_lesion_pd,
                "tract_pd_lesion_binary": stats_tract_pd_lesion_binary,
                "tract_pd_lesion_dilated": stats_tract_pd_lesion_dilated,
                "tract_pd_lesion_penumbra": stats_tract_pd_lesion_penumbra,
                "tract_pd_change_lesion_binary": change_tract_pd_lesion,
            }

            # Add lesion load metrics to row
            for detail_key, detail_val in tract_lesion_load_details.items():
                tract_row[f"{lesion_name}_{detail_key}"] = detail_val
            
            # Add distance metrics to row
            for distance_key, distance_val in tract_lesion_distance_metrics.items():
                tract_row[f"{lesion_name}_{distance_key}"] = distance_val

            # Add tract-lesion metrics to row 
            for suffix, metrics in tract_lesion_metrics.items():
                for metric_name, metric_value in metrics.items():
                    tract_row[f"{lesion_name}_{suffix}_{metric_name}"] = metric_value
            

            for map_name, map_array in inputs.map_arrays.items():
                # Tract-lesion map stats.
                map_tract_pd = map_array * tract_array
                map_tract_pd_lesion_pd = map_array * tract_pd_lesion_pd
                map_tract_binary_lesion_binary = map_array * tract_binary_lesion_binary
                map_tract_pd_dilated_lesion = map_array * tract_pd_lesion_dilated
                map_tract_binary_dilated_lesion = map_array * tract_binary_lesion_dilated
                map_tract_pd_penumbra = map_array * tract_pd_lesion_penumbra
                map_tract_binary_penumbra = map_array * tract_binary_penumbra

                map_tract_whole_overlapping_lesions_binary = map_array * whole_overlapping_lesions_binary

                # Map x Tract PD x Lesion PD
                stats_map_tract_pd_lesion_pd = _summary_nonzero(map_tract_pd_lesion_pd)
                # Map x Tract Binary x Lesion binary
                stats_map_tract_binary_lesion_binary = _summary_nonzero(map_tract_binary_lesion_binary)
                # Map x Tract PD x Dilated lesion
                stats_map_tract_pd_dilated_lesion = _summary_nonzero(map_tract_pd_dilated_lesion)
                # Map x Tract Binary x Dilated lesion
                stats_map_tract_binary_dilated_lesion = _summary_nonzero(map_tract_binary_dilated_lesion)
                # Map x Tract PD x Penumbra
                stats_map_tract_pd_penumbra = _summary_nonzero(map_tract_pd_penumbra)
                # Map x Tract Binary x Penumbra
                stats_map_tract_binary_penumbra = _summary_nonzero(map_tract_binary_penumbra)
                # Map x Tract Binary x Whole overlapping lesions
                stats_map_tract_binary_whole_overlapping_lesions = _summary_nonzero(map_tract_whole_overlapping_lesions_binary)
                
                # Change: Map x Tract PD x Lesion binary
                change_map_tract_pd_lesion = _summary_change_nonzero(
                    lesion_map=lesion_array, values_b=map_tract_pd)
                
                tract_lesion_metrics_map = {
                    "map_tract_pd_lesion_pd": stats_map_tract_pd_lesion_pd,
                    "map_tract_binary_lesion_binary": stats_map_tract_binary_lesion_binary,
                    "map_tract_pd_dilated_lesion": stats_map_tract_pd_dilated_lesion,
                    "map_tract_binary_dilated_lesion": stats_map_tract_binary_dilated_lesion,
                    "map_tract_pd_penumbra": stats_map_tract_pd_penumbra,
                    "map_tract_binary_penumbra": stats_map_tract_binary_penumbra,
                    "map_tract_pd_change_lesion_binary": change_map_tract_pd_lesion,
                    "map_tract_binary_whole_overlapping_lesions": stats_map_tract_binary_whole_overlapping_lesions,
                }
                # Add map-tract-lesion metrics to row with lesion and map-specific prefix.
                for suffix, metrics in tract_lesion_metrics_map.items():
                    for metric_name, metric_value in metrics.items():
                        tract_row[f"{map_name}_{lesion_name}_{suffix}_{metric_name}"] = metric_value
                
        for map_name, map_array in inputs.map_arrays.items():
            # Sanity check if there are negative values in the map array, which should not be the case for probability maps.
            if np.any(map_array < 0):
                logger.warning(f"Warning: Map '{map_name}' contains negative values; check input data.")

            # Whole-tract map stats.
            map_tract_pd = map_array * tract_array
            map_tract_binary = map_array * inputs.tract_binary_arrays[tract_name]
            map_tract_skeletonised = inputs.tract_skeletons_arrays[tract_name] * map_array

            # Map x Tract PD
            stats_map_tract_pd = _summary_nonzero(map_tract_pd)
            # Map x Tract Binary
            stats_map_tract_binary = _summary_nonzero(map_tract_binary)
            # Map x Tract Skeleton
            stats_map_tract_skeletonised = _summary_nonzero(map_tract_skeletonised)

            tract_array_metrics_map = {
                "map_tract_pd": stats_map_tract_pd,
                "map_tract_binary": stats_map_tract_binary,
                "map_tract_skeletonised": stats_map_tract_skeletonised,
            }

            # Add map-tract metrics to row with map-specific prefix.
            for suffix, metrics in tract_array_metrics_map.items():
                for metric_name, metric_value in metrics.items():
                    tract_row[f"{map_name}_{suffix}_{metric_name}"] = metric_value

        tract_level_metrics_df = pd.concat(
            [tract_level_metrics_df, pd.DataFrame([tract_row])],
            ignore_index=True,
        )

    return tract_level_metrics_df

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

    # Load critical arrays once for reuse across metric families ───────────
    brain_array: Optional[np.ndarray] = None

    lesion_arrays: Dict[str, np.ndarray] = {}
    binary_lesion_arrays: Dict[str, np.ndarray] = {}
    dilated_lesion_arrays: Dict[str, np.ndarray] = {}
    penumbras_arrays: Dict[str, np.ndarray] = {}

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

        dilated_array = _dilate_cropped(binary_array, iterations=3)
        penumbra_array = np.where(dilated_array - binary_array > 0, 1, 0)

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
        tract_arrays=tract_arrays,
        thresholded_tract_arrays=thresholded_tract_arrays,
        tract_binary_arrays=tract_binary_arrays,
        tract_skeletons_arrays=tract_skeletons_arrays,
        map_arrays=map_arrays,
        skeleton_array=skeleton_array,
        voxel_vol=voxel_vol,
        tract_atlas_files_dict=tract_atlas_files_dict,
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

    # Fast (sparse) path guard: see the FAST PATH note above _flat(). Falls back
    # to the dense reference implementation if any input holds NaN/Inf.
    use_fast = _inputs_are_finite(
        list(inputs.map_arrays.values()) + list(inputs.tract_arrays.values())
        + list(inputs.lesion_arrays.values())
    )
    if not use_fast:
        logger.warning("Non-finite values found in an input array; using the dense (slow) metrics path.")
    # Per-lesion-set supports / labelling / distance transform, built once and
    # shared by every consumer below (whole-brain, tract-level).
    lesion_contexts: Dict[str, _LesionContext] = {}
    if use_fast:
        for _ln in inputs.lesion_arrays:
            lesion_contexts[_ln] = _LesionContext(
                inputs.lesion_arrays[_ln], inputs.binary_lesion_arrays[_ln],
                inputs.dilated_lesion_arrays[_ln], inputs.penumbras_arrays[_ln],
            )

    # ── Metric family: whole-brain ───────────────────────────────────────────
    wb_df = _compute_whole_brain_metrics(
        subject=args.subject,
        brain_array=inputs.brain_array,
        lesion_dict=inputs.lesion_arrays,
        voxel_vol_mm3=inputs.voxel_vol,
    )
 
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
    _wb_lesion_fn = (_compute_whole_brain_lesion_microstructure_metrics_fast if use_fast
                     else _compute_whole_brain_lesion_microstructure_metrics)
    _wb_kwargs = {"lesion_contexts": lesion_contexts} if use_fast else {}
    wb_lesion_df = _wb_lesion_fn(
        subject=args.subject,
        map_arrays=inputs.map_arrays,
        lesion_arrays=inputs.lesion_arrays,
        binary_lesion_arrays=inputs.binary_lesion_arrays,
        penumbras_arrays=inputs.penumbras_arrays,
        dilated_lesion_arrays=inputs.dilated_lesion_arrays,
        **_wb_kwargs,
    )
    wb_lesion_csv = metrics_dir / "whole_brain_lesion_metrics.csv"
    print("✓ Whole Brain Lesion Microstructure Metrics CSV saved")
    write_csv(wb_lesion_df, str(wb_lesion_csv))

    # ── Peak width of Whole skeleton of tracts ───────────────────────────────────────────
    # Always written (empty DataFrame when no skeleton) so Snakemake can require the file.
    skeleton_metrics_df = pd.DataFrame()
    if inputs.skeleton_array is not None and use_fast and _inputs_are_finite([inputs.skeleton_array]):
        skeleton_metrics_df = _skeleton_metrics_fast(args.subject, inputs.map_arrays, inputs.skeleton_array)
    elif inputs.skeleton_array is not None:
        rows: List[Dict[str, Any]] = []
        for map_name, map_array in inputs.map_arrays.items():
            skeleton_map = map_array * inputs.skeleton_array
            stats_skeleton_map = _summary_nonzero(skeleton_map)

            row = {"subject": args.subject, "map_name": map_name}
            for stat_name, stat_value in stats_skeleton_map.items():
                row[f"skeleton_{stat_name}"] = stat_value

            rows.append(row)

        skeleton_metrics_df = pd.DataFrame(rows)

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
            if _extension is not None and hasattr(_extension, "on_lesion_ready"):
                labeled_lesions, n_components = ndimage.label(inputs.binary_lesion_arrays[lesion_name])
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
    if use_fast:
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
    else:
        tract_level_metrics_df = _compute_tract_level_metrics_legacy(args, inputs)

    # Write tract-level metrics.
    tract_csv = metrics_dir / "tract_level_metrics.csv"
    print("✓ Tract-Level Metrics CSV saved")
    write_csv(tract_level_metrics_df, str(tract_csv))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
