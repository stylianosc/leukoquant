"""Choose which DWI archives of one ADNI session go into the pipeline.

ADNI sometimes ships several DWI archives for one session: the same acquisition twice (``ORIG_Axial_DTI`` and the corrected
``Axial_DTI``), a truncated first attempt followed by the full run, scans with a different field of view or slab, or single
b0 volumes acquired with the opposite phase-encoding direction (``b0_AP``/``b0_PA``). TRACULA concatenates all scans it is
given voxel by voxel, without registering them to each other, and then runs one eddy correction over the result. Scans
that are not on the same voxel grid would corrupt the data, and repeated scans only count the same measurements twice.
In every session of ADNI3 where this happens, the smaller scan is exactly the start of the larger scan's gradient table
(a truncated first attempt), so stacking added no new directions; such scans are dropped.

The names of the files inside the archives are not reliable for this (the time stamp is when ADNI processed the data,
not when it was scanned), so every decision uses what the converter produces: matrix, voxel size, position,
phase-encoding direction, acquisition time and the gradient table.

Rules, applied to the ADNI archives of one session (other inputs are never touched):
  1. Single b0 volumes (no diffusion weighting, e.g. ``b0_AP``/``b0_PA``) are dropped, unless ``topup="auto"`` and there
     are two of them with opposite phase-encoding directions on the same matrix and voxel size as the diffusion scan;
     then one of each direction is kept for TRACULA's reverse-polarity correction (``dob0 = 2``).
  2. Of the remaining scans one *main* scan is chosen, ranked by: most diffusion volumes; highest slab top (the brain
     rather than the neck, tops within ``SLAB_TOL_MM`` count as equal); isotropic voxels; not an "ORIG" copy; latest
     acquisition time; latest series number.
  3. Any other scan with the same matrix, voxel size, position (origin within ``ORIGIN_TOL_MM``) and phase-encoding
     direction as the main scan is stacked with it only if it adds at least one diffusion direction (same b-value shell,
     more than ``DIRECTION_TOL_DEG`` away from every direction already kept). Repeats, prefixes and ORIG copies add
     none and are dropped; scans on another grid, slab or phase-encoding direction are always dropped.
A session whose kept scan has fewer than six diffusion volumes cannot give a tensor fit; this is logged as a warning.
"""
from __future__ import annotations

import json
import logging
import re
import shutil
import tempfile
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# ADNI DICOM names, e.g. ADNI_021_S_7045_MR_Axial_DTI__br_raw_20220223120056490_2300_S1107498_I1547797.dcm
_DICOM_NAME = re.compile(r"MR_(?P<series>.+?)_br_raw_\d+_\d+_S(?P<series_no>\d+)_I")

# Tolerances (millimetres): voxel sizes that count as equal, and slab tops that count as the same height
VOXEL_TOL_MM = 0.01
SLAB_TOL_MM = 10.0
# Scans are stacked only when they sit at the same position; origins within this distance count as the same
ORIGIN_TOL_MM = 1.0
# Gradient directions closer than this (in degrees) count as the same direction, as do b-values within B_SHELL_TOL
DIRECTION_TOL_DEG = 5.0
B_SHELL_TOL = 60.0
# Fewest diffusion volumes from which a diffusion tensor can be fitted
MIN_DIFFUSION_VOLUMES = 6
# b-values below this count as b0 (no diffusion weighting)
B0_MAX = 50.0


@dataclass(frozen=True)
class ZipInfo:
    path: str
    series: str
    series_no: int


@dataclass(frozen=True)
class Scan:
    """What the converter produced for one archive."""
    shape: tuple         # full image shape, 3D for a single volume
    zooms: tuple         # voxel size in mm (x, y, z)
    origin: tuple        # world position of voxel (0, 0, 0) in mm
    bvals: tuple         # rounded b-values, one per volume
    bvecs: tuple         # rounded gradient directions, flattened
    pe_dir: str          # dcm2niix PhaseEncodingDirection ("j", "j-", "i", ...), "" if unknown
    acq_time: float = 0.0  # acquisition time in seconds (DICOM AcquisitionDateTime, else AcquisitionTime); 0 if unknown

    @property
    def n_volumes(self) -> int:
        return len(self.bvals)

    @property
    def n_diffusion(self) -> int:
        """Number of diffusion-weighted volumes (b above the b0 threshold)."""
        return int(sum(b > B0_MAX for b in self.bvals))

    def directions(self) -> List[Tuple[float, np.ndarray]]:
        """(b-value, unit gradient direction) of every diffusion-weighted volume."""
        vecs = np.asarray(self.bvecs, dtype=float).reshape(-1, 3)
        out = []
        for b, v in zip(self.bvals, vecs):
            norm = np.linalg.norm(v)
            if b > B0_MAX and norm > 0:
                out.append((float(b), v / norm))
        return out

    @property
    def isotropic(self) -> bool:
        return max(self.zooms) - min(self.zooms) < VOXEL_TOL_MM

    @property
    def slab_top(self) -> float:
        """World z of the top slice (assumes an axial grid, z increasing with the slice index)."""
        return self.origin[2] + self.zooms[2] * (self.shape[2] - 1)


def read_zip_info(zip_path: str) -> Optional[ZipInfo]:
    """Series name and series number from the DICOM names inside a zip.

    Returns None when the input is not a zip or its names do not follow the ADNI pattern; such inputs are never
    touched by the selection.
    """
    if not zip_path.lower().endswith(".zip"):
        return None
    try:
        with zipfile.ZipFile(zip_path) as zf:
            names = sorted(n for n in zf.namelist() if not n.endswith("/"))
    except (zipfile.BadZipFile, OSError):
        return None
    for name in names:
        match = _DICOM_NAME.search(name)
        if match:
            return ZipInfo(zip_path, match["series"], int(match["series_no"]))
    return None


def _per_volume_bvecs(array: np.ndarray, n_volumes: int) -> tuple:
    """Gradient directions as flat (x, y, z) triplets, one per volume, rounded.

    dcm2niix writes FSL layout (3 rows of n values); TRACULA-style files have n rows of 3. Both are accepted. A
    single volume loads as a flat array of 3.
    """
    arr = np.asarray(array, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1) if arr.size == 3 else arr.reshape(3, -1)
    if arr.shape[0] == 3 and arr.shape[1] == n_volumes:
        arr = arr.T
    return tuple(np.round(arr.reshape(-1), 3))


def _acquisition_seconds(sidecar: dict) -> float:
    """Acquisition time from a dcm2niix sidecar as seconds (0.0 when absent).

    dcm2niix writes ``AcquisitionTime`` without zero padding for some scanners (e.g. ``16:30:7.117500``), so the
    parts are parsed as numbers instead of with a fixed format. A full date-time is preferred when present so that
    scans on different days sort correctly.
    """
    text = str(sidecar.get("AcquisitionDateTime") or sidecar.get("AcquisitionTime") or "")
    match = re.search(r"(?:(\d{4})-(\d{1,2})-(\d{1,2})T)?(\d{1,2}):(\d{1,2}):(\d{1,2}(?:\.\d+)?)", text)
    if not match:
        return 0.0
    year, month, day, hour, minute, second = match.groups()
    days = 0
    if year:
        from datetime import date
        days = date(int(year), int(month), int(day)).toordinal()
    return days * 86400.0 + int(hour) * 3600 + int(minute) * 60 + float(second)


def convert_scan(zip_path: str) -> Scan:
    """Convert one archive in a scratch directory and describe the result."""
    # Imported here so that this module stays importable without the converter stack.
    import nibabel as nib
    from leukoquant.utils.dwi_utils import _convert_dicom_with_fallback

    work = Path(tempfile.mkdtemp())
    try:
        src, out = work / "src", work / "out"
        src.mkdir()
        out.mkdir()
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(src)
        dwi, bvec, bval = _convert_dicom_with_fallback(src, out)
        img = nib.load(str(dwi))
        sidecar = Path(str(dwi).replace(".nii.gz", ".json").replace(".nii", ".json"))
        pe_dir, acq_time = "", 0.0
        if sidecar.is_file():
            meta = json.loads(sidecar.read_text())
            pe_dir = str(meta.get("PhaseEncodingDirection", ""))
            acq_time = _acquisition_seconds(meta)
        bvals_arr = np.round(np.loadtxt(str(bval)).reshape(-1), 1) if bval else np.array([])
        return Scan(
            shape=tuple(img.shape),
            zooms=tuple(round(float(z), 3) for z in img.header.get_zooms()[:3]),
            origin=tuple(round(float(o), 2) for o in img.affine[:3, 3]),
            bvals=tuple(bvals_arr),
            bvecs=_per_volume_bvecs(np.loadtxt(str(bvec)), len(bvals_arr)) if bvec else (),
            pe_dir=pe_dir,
            acq_time=acq_time,
        )
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _same_matrix_and_voxels(a: Scan, b: Scan) -> bool:
    return a.shape[:3] == b.shape[:3] and np.allclose(a.zooms, b.zooms, atol=VOXEL_TOL_MM)


def _is_orig(info: ZipInfo) -> bool:
    return "orig" in info.series.lower()


def _same_position(a: Scan, b: Scan) -> bool:
    """Same matrix, voxel size, position and phase-encoding direction: the two can be stacked voxel by voxel."""
    return (
        _same_matrix_and_voxels(a, b)
        and np.allclose(a.origin, b.origin, atol=ORIGIN_TOL_MM)
        and a.pe_dir == b.pe_dir
    )


def _adds_directions(kept: List[Tuple[float, np.ndarray]], scan: Scan) -> bool:
    """True when ``scan`` has a diffusion direction (at its b-value) that ``kept`` does not already contain."""
    cos_tol = np.cos(np.radians(DIRECTION_TOL_DEG))
    for b, v in scan.directions():
        # +v and -v are the same axis of diffusion
        if not any(abs(b - kb) <= B_SHELL_TOL and abs(float(v @ kv)) >= cos_tol for kb, kv in kept):
            return True
    return False


def _opposite(pe_dir: str) -> str:
    """Opposite phase-encoding direction: ``j`` <-> ``j-``, ``i`` <-> ``i-``."""
    return pe_dir[:-1] if pe_dir.endswith("-") else pe_dir + "-"


def _best_scan(diffusion: List[Tuple[ZipInfo, Scan]]) -> Tuple[ZipInfo, Scan]:
    """The one diffusion scan to use; see rule 2 in the module docstring."""
    most = max(s.n_diffusion for _, s in diffusion)
    pool = [g for g in diffusion if g[1].n_diffusion == most]
    top = max(s.slab_top for _, s in pool)
    pool = [g for g in pool if top - g[1].slab_top <= SLAB_TOL_MM]
    return max(pool, key=lambda g: (g[1].isotropic, not _is_orig(g[0]), g[1].acq_time, g[0].series_no))


def _reverse_polarity_pair(b0_scans: List[Tuple[ZipInfo, Scan]], main: Scan) -> List[Tuple[ZipInfo, Scan]]:
    """One single-b0 scan of each of two opposite phase-encoding directions, on the matrix and voxel size of ``main``.

    Returns [] when there is no such pair. When a direction occurs more than once the latest acquisition is used.
    """
    usable = [g for g in b0_scans if g[1].pe_dir and _same_matrix_and_voxels(g[1], main)]
    latest: Dict[str, Tuple[ZipInfo, Scan]] = {}
    for g in sorted(usable, key=lambda g: (g[1].acq_time, g[0].series_no)):
        latest[g[1].pe_dir] = g
    for pe_dir, g in latest.items():
        if _opposite(pe_dir) in latest:
            return [g, latest[_opposite(pe_dir)]]
    return []


@dataclass
class Selection:
    """Archives to use (original order) and whether they include a reverse-polarity b0 pair for TRACULA's ``dob0 = 2``."""
    paths: List[str]
    reverse_polarity: bool = False
    dropped: Dict[str, str] = field(default_factory=dict)


def select_scans(
    candidates: List[str],
    scan_fn: Callable[[str], Scan] = convert_scan,
    topup: str = "none",
) -> Selection:
    """Choose the archives to use; see the module docstring for the rules.

    ``scan_fn`` converts an archive and describes it; it is a parameter so that tests do not need a DICOM converter.
    An archive that cannot be read or converted is kept untouched so that the pipeline reports the problem as before.
    ``topup`` is ``"none"`` (default) or ``"auto"``.
    """
    if topup not in ("none", "auto"):
        raise ValueError(f"topup must be 'none' or 'auto', got {topup!r}")

    infos: Dict[str, ZipInfo] = {}
    scans: Dict[str, Scan] = {}
    for path in candidates:
        info = read_zip_info(path)
        if not info:
            continue
        try:
            scans[path] = scan_fn(path)
        except Exception as exc:  # the pipeline will report an unusable archive itself
            logger.warning("Could not describe DWI archive %s (%s); leaving it in the list.", path, exc)
            continue
        infos[path] = info

    diffusion = [(infos[p], scans[p]) for p in scans if scans[p].n_diffusion > 0]
    b0_only = [(infos[p], scans[p]) for p in scans if scans[p].n_diffusion == 0]
    if len(scans) < 2 or not diffusion:  # nothing to choose from; leave the list alone
        return Selection(list(candidates))

    best_info, best = _best_scan(diffusion)
    keep = {best_info.path}
    # Rule 3: further scans on the same grid are stacked only when they add diffusion directions.
    known = best.directions()
    others = sorted((g for g in diffusion if g[0].path != best_info.path),
                    key=lambda g: (-g[1].n_diffusion, -g[1].acq_time, -g[0].series_no))
    for info, scan in others:
        if _same_position(best, scan) and _adds_directions(known, scan):
            keep.add(info.path)
            known = known + scan.directions()
    reverse_polarity = False
    if topup == "auto":
        pair = _reverse_polarity_pair(b0_only, best)
        keep.update(info.path for info, _ in pair)
        reverse_polarity = bool(pair)

    dropped: Dict[str, str] = {}
    for info, scan in diffusion + b0_only:
        if info.path in keep:
            continue
        if scan.n_diffusion == 0:
            dropped[info.path] = "single b0 volume without diffusion directions"
        else:
            dropped[info.path] = (
                f"adds no new diffusion directions or is on another grid, slab or phase-encoding direction "
                f"than {Path(best_info.path).name}"
            )
    for path, reason in dropped.items():
        logger.warning("DWI archive %s (%s) not used: %s.", Path(path).name, infos[path].series, reason)
    if best.n_diffusion < MIN_DIFFUSION_VOLUMES:
        logger.warning(
            "DWI archive %s has only %d diffusion volumes; a diffusion tensor needs at least %d.",
            Path(best_info.path).name, best.n_diffusion, MIN_DIFFUSION_VOLUMES,
        )
    return Selection([p for p in candidates if p not in dropped], reverse_polarity, dropped)
