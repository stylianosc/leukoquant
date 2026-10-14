"""DWI input preparation utilities.

Provides a single entrypoint `prepare_dwi_input` that accepts a DWI input
which can be:
- a .nii.gz file
- a .nii file (will be compressed to .nii.gz)
- a .zip archive containing a DICOM series
- a directory containing DICOM files

The function returns a dict with keys `dwi`, `bvecs`, `bvals`, and
`tempdir` (tempdir is None if no temporary workspace was created).

It raises informative exceptions when required outputs (bvec/bval/dwi)
cannot be produced.
"""

from __future__ import annotations

import collections
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path
from typing import Optional, Dict

import logging
import numpy as np
import argparse

logger = logging.getLogger(__name__)

def _find_sidecars(folder: Path, stem: Optional[str] = None) -> (Optional[Path], Optional[Path]):
    # dcm2niix emits singular "*.bvec"/"*.bval"; TRACULA's own preprocessed
    # DWI output uses the plural FSL convention "*.bvecs"/"*.bvals" (e.g.
    # dwi.bvecs, dwi.bvals) -- match both so sidecar auto-discovery works
    # for either source.
    #
    # When the caller knows the DWI's own filename stem (e.g. "dwi" for
    # dwi.nii.gz), match sidecars by exact stem instead of a bare glob --
    # TRACULA's own dmri/ output directory can contain intermediate
    # artifacts like dwi_orig.1.bvals (a pre-correction file from an
    # earlier trac-preproc step) alongside the real dwi.bvals, and a plain
    # "*.bvals" glob can't tell them apart even though only one is the
    # actual sidecar for dwi.nii.gz. A prefix match isn't enough either
    # ("dwi_orig.1.bvals" does start with "dwi") -- only an exact
    # "{stem}.bval(s)" name is unambiguous.
    if stem:
        bvecs = [p for p in (folder / f"{stem}.bvec", folder / f"{stem}.bvecs") if p.is_file()]
        bvals = [p for p in (folder / f"{stem}.bval", folder / f"{stem}.bvals") if p.is_file()]
    else:
        bvecs = list(folder.glob("*.bvec")) + list(folder.glob("*.bvecs"))
        bvals = list(folder.glob("*.bval")) + list(folder.glob("*.bvals"))

    if len(bvecs) > 1 or len(bvals) > 1:
        raise RuntimeError(f"Multiple bvec or bval files found in {folder}, cannot disambiguate.")

    return (bvecs[0] if bvecs else None), (bvals[0] if bvals else None)


# A diffusion tensor needs at least 6 directions; b0 volumes do not count.
MIN_DW_VOLUMES = 6


class InputDataError(ValueError):
    """The raw DWI itself is unusable (as opposed to a transient or environment
    failure). Only this error type leaves a failure marker, so a one-off disk
    or network problem can never get a healthy session skipped for good."""


def _find_converter_candidates(dcm2niix_primary: Optional[str] = None) -> list[tuple[str, str]]:
    """Return an ordered list of (converter_name, executable_path) to try.

    Discovers converters portably without machine-specific hardcoded paths:
      1. Primary dcm2niix (explicit path, module-relative external/noddi, or PATH)
      2. Rebuilt dcm2niix_nonspatial (module-relative external/noddi or PATH)
      3. FreeSurfer mri_convert (PATH or FREESURFER_HOME/bin)
    """
    candidates: list[tuple[str, str]] = []

    # Module-relative path to external/noddi
    module_noddi = Path(__file__).resolve().parent.parent / "external" / "noddi"

    # 1. Primary dcm2niix
    p = None
    if dcm2niix_primary and Path(dcm2niix_primary).is_file() and os.access(dcm2niix_primary, os.X_OK):
        p = dcm2niix_primary
    elif (module_noddi / "dcm2niix").is_file() and os.access(module_noddi / "dcm2niix", os.X_OK):
        p = str(module_noddi / "dcm2niix")
    elif shutil.which("dcm2niix"):
        p = shutil.which("dcm2niix")

    if p:
        candidates.append(("dcm2niix", p))

    # 2. Rebuilt dcm2niix_nonspatial
    nonspatial_candidates = []
    if dcm2niix_primary:
        nonspatial_candidates.append(str(Path(dcm2niix_primary).parent / "dcm2niix_nonspatial"))
    nonspatial_candidates.append(str(module_noddi / "dcm2niix_nonspatial"))
    if shutil.which("dcm2niix_nonspatial"):
        nonspatial_candidates.append(shutil.which("dcm2niix_nonspatial"))

    for nsc in nonspatial_candidates:
        if nsc and Path(nsc).is_file() and os.access(nsc, os.X_OK):
            candidates.append(("dcm2niix_nonspatial", nsc))
            break

    # 3. FreeSurfer mri_convert
    mc_candidates = []
    if shutil.which("mri_convert"):
        mc_candidates.append(shutil.which("mri_convert"))
    fs_home = os.environ.get("FREESURFER_HOME")
    if fs_home:
        mc_candidates.append(str(Path(fs_home) / "bin" / "mri_convert"))

    for mc in mc_candidates:
        if mc and Path(mc).is_file() and os.access(mc, os.X_OK):
            candidates.append(("mri_convert", mc))
            break

    return candidates


def _extract_dicom_slice_location(dcm_path: Path) -> Optional[float]:
    """Slice position of a DICOM file: SliceLocation, else the z of ImagePositionPatient.

    Reads the header only (no pixel data). Returns None for files that are not
    readable DICOMs or carry neither tag.
    """
    import pydicom

    try:
        ds = pydicom.dcmread(str(dcm_path), stop_before_pixels=True,
                             specific_tags=["SliceLocation", "ImagePositionPatient"])
        if "SliceLocation" in ds:
            return round(float(ds.SliceLocation), 3)
        if "ImagePositionPatient" in ds:
            return round(float(ds.ImagePositionPatient[2]), 3)
    except Exception:  # not a DICOM, or a malformed header: treat as unlocatable
        pass
    return None


def _filter_incomplete_slices(src_dir: Path) -> Path:
    """Detect and remove partial/trailing volume slices from an aborted acquisition.

    If an acquisition was interrupted midway through a volume, some slice locations
    will have fewer files than the complete volumes. dcm2niix will reject such datasets
    as inconsistent or non-spatial. Filtering out trailing slices restores volume consistency.

    Returns the directory to use (either src_dir or a temporary directory with symlinks).
    """
    # Use rglob to find all DICOMs recursively; convert to a list of unique
    # resolved paths (dict preserves insertion order, deduplication is via key lookup).
    seen: dict[Path, None] = {}
    for p in src_dir.rglob("*.dcm"):
        seen[p.resolve()] = None
    dcms: list[Path] = list(seen.keys())
    if not dcms:
        # No .dcm extension: treat all non-hidden files as potential DICOMs
        for p in src_dir.rglob("*"):
            if p.is_file() and not p.name.startswith("."):
                seen[p.resolve()] = None
        dcms = list(seen.keys())

    if len(dcms) < 20:
        return src_dir

    loc_counts: collections.Counter[float] = collections.Counter()
    file_locs: list[tuple[Path, Optional[float]]] = []
    for d in dcms:
        loc = _extract_dicom_slice_location(d)
        file_locs.append((d, loc))
        if loc is not None:
            loc_counts[loc] += 1

    if not loc_counts:
        return src_dir

    counts_freq = collections.Counter(loc_counts.values())
    expected_vol_count = counts_freq.most_common(1)[0][0]

    # Check if there are slice locations with partial volumes
    outliers = {loc for loc, count in loc_counts.items() if count < expected_vol_count}
    if not outliers:
        return src_dir

    logger.warning(
        "Detected incomplete volume slices in %s! %d slice locations have partial acquisitions. "
        "Creating filtered view with complete volumes only (%d full volumes).",
        src_dir,
        len(outliers),
        expected_vol_count,
    )
    filtered_dir = Path(tempfile.mkdtemp(prefix="filtered_dcm_"))
    kept = 0
    # Keep up to expected_vol_count files per valid slice location
    loc_seen: collections.Counter[float] = collections.Counter()
    for d, loc in file_locs:
        if loc is not None and loc not in outliers:
            if loc_seen[loc] < expected_vol_count:
                loc_seen[loc] += 1
                (filtered_dir / d.name).symlink_to(d.resolve())
                kept += 1
    logger.info("Filtered %d partial/excess slices; retained %d complete slices.", len(dcms) - kept, kept)
    return filtered_dir


def _convert_dicom_with_fallback(
    src_dir: Path,
    out_dir: Path,
    dcm2niix_path: Optional[str] = None,
    verbose: bool = False,
) -> tuple[Path, Optional[Path], Optional[Path]]:
    """Convert a DICOM directory to NIfTI with sidecars using a fallback chain.

    Fallback sequence:
      1. dcm2niix (standard)
      2. dcm2niix_nonspatial (built with -DmyInstanceNumberOrderIsNotSpatial)
      3. mri_convert (FreeSurfer)

    Also handles incomplete acquisitions (aborted scans with partial slice sets)
    by filtering trailing slices if standard conversion fails.

    Validates that the converted volume is not a corrupted 3D slice-dump when
    b-values indicate a multi-volume 4D DWI acquisition.

    Returns:
        (dwi_nii_path, bvec_path, bval_path)
    """
    candidates = _find_converter_candidates(dcm2niix_path)
    if not candidates:
        raise RuntimeError("No DICOM converter found (checked dcm2niix, dcm2niix_nonspatial, mri_convert).")

    errors: list[str] = []

    # First pass: try conversion directly on src_dir
    # Second pass: if dcm2niix failed with spatial/ordering or slice count errors, try filtering incomplete slices
    filtered_dir: Optional[Path] = None

    for attempt_dir in [src_dir, "FILTERED"]:
        current_src = attempt_dir
        if attempt_dir == "FILTERED":
            filtered_dir = _filter_incomplete_slices(src_dir)
            if filtered_dir == src_dir:
                break
            current_src = filtered_dir

        for name, exec_path in candidates:
            logger.info("Attempting DICOM conversion with %s (%s) on %s", name, exec_path, current_src)

            # Clean out_dir before each attempt to avoid mixing partial outputs
            for child in out_dir.iterdir():
                if child.is_file() or child.is_symlink():
                    child.unlink()
                elif child.is_dir():
                    shutil.rmtree(child, ignore_errors=True)

            try:
                if name in ("dcm2niix", "dcm2niix_nonspatial"):
                    cmd = [exec_path, "-z", "y", "-b", "y", "-o", str(out_dir), str(current_src)]
                    if verbose:
                        print(f"Running {name}: {' '.join(cmd)}")
                    proc = subprocess.run(cmd, capture_output=True, text=True)
                    if proc.returncode != 0:
                        raise RuntimeError(f"{name} exited with code {proc.returncode}. stderr: {proc.stderr.strip()}")

                elif name == "mri_convert":
                    dcms = sorted(current_src.glob("*.dcm")) + sorted(current_src.rglob("*.dcm"))
                    if not dcms:
                        dcms = [f for f in current_src.rglob("*") if f.is_file()]
                    if not dcms:
                        raise FileNotFoundError(f"No DICOM files found in {current_src} for mri_convert")
                    first_dcm = dcms[0]
                    out_nii = out_dir / "data.nii.gz"
                    cmd = [exec_path, str(first_dcm), str(out_nii)]
                    if verbose:
                        print(f"Running mri_convert: {' '.join(cmd)}")
                    env = os.environ.copy()
                    if "FREESURFER_HOME" not in env:
                        cand_fs_home = Path(exec_path).resolve().parent.parent
                        if (cand_fs_home / "bin" / "mri_convert").is_file():
                            env["FREESURFER_HOME"] = str(cand_fs_home)
                    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
                    if proc.returncode != 0:
                        raise RuntimeError(f"mri_convert exited with code {proc.returncode}. stderr: {proc.stderr.strip()}")

                # Find produced NIfTI files
                niis = list(out_dir.glob("*.nii.gz")) + list(out_dir.glob("*.nii"))
                if not niis:
                    raise FileNotFoundError(f"{name} produced no NIfTI files.")
                dwi_cand = niis[0]

                # Find sidecars in out_dir
                bvec_cand, bval_cand = _find_sidecars(out_dir)

                # Validate converted image structure
                import nibabel as nib
                img = nib.load(str(dwi_cand))

                # A DWI without gradient tables is unusable (mri_convert writes
                # the image only), so it must not count as a successful conversion.
                if not (bval_cand and bval_cand.exists() and bvec_cand and bvec_cand.exists()):
                    raise InputDataError(f"{name} produced no bval/bvec files.")

                if bval_cand and bval_cand.exists():
                    try:
                        bvals_arr = np.loadtxt(str(bval_cand)).reshape(-1)
                        n_bvals = len(bvals_arr)
                        if img.ndim < 4 and n_bvals > 1:
                            raise InputDataError(
                                f"Corrupted conversion: image is {img.ndim}D {img.shape} "
                                f"but {n_bvals} b-values were extracted (all slices stacked into 3D)."
                            )
                        if img.ndim >= 4 and img.shape[-1] != n_bvals:
                            raise InputDataError(
                                f"Volume mismatch: image has {img.shape[-1]} volumes but {n_bvals} b-values."
                            )
                    except Exception as val_err:
                        if isinstance(val_err, InputDataError):
                            raise
                        logger.debug("Could not parse bvals for validation: %s", val_err)

                logger.info("DICOM conversion succeeded using %s: %s", name, dwi_cand.name)
                if filtered_dir and filtered_dir != src_dir and filtered_dir.exists():
                    shutil.rmtree(filtered_dir, ignore_errors=True)
                return dwi_cand, bvec_cand, bval_cand

            except Exception as e:
                logger.warning("Conversion with %s failed: %s", name, e)
                errors.append(f"{name} (attempt on {current_src.name}): {e}")

    if filtered_dir and filtered_dir != src_dir and filtered_dir.exists():
        shutil.rmtree(filtered_dir, ignore_errors=True)

    raise InputDataError("All DICOM converters failed:\n" + "\n".join(f"  - {err}" for err in errors))


def _run_dcm2niix(dcm2niix_exec: str, src_dir: Path, out_dir: Path, verbose: bool = False) -> None:
    """Run DICOM conversion with fallback chain (retaining legacy function signature)."""
    _convert_dicom_with_fallback(src_dir, out_dir, dcm2niix_path=dcm2niix_exec, verbose=verbose)


def _standardize_bvals(bvals: np.ndarray) -> np.ndarray:
    """Ensure b-values are in s/mm² and round to the nearest 100.

    Clinical DWI b-values in s/mm² are typically 0, 1000, 2000, 3000.
    Some scanners or converters produce values in non-standard units:
      - Scaled ×1000 (ms/mm² or similar): 0, 1000000, 2000000 → divide by 1000
      - Floating-point drift within a shell: 995, 1005 → round to 1000

    Detection heuristic: if the maximum non-zero b-value exceeds 10 000,
    the values are assumed to be ×1000 too large and are divided by 1000
    before rounding.  The original array is never mutated.
    """
    bvals = bvals.copy().astype(float)
    non_zero = bvals[bvals > 0]
    if non_zero.size > 0 and non_zero.max() > 10_000:
        logger.warning(
            "b-values appear to be in non-standard units (max non-zero = %.1f). "
            "Dividing by 1000 to convert to s/mm².",
            non_zero.max(),
        )
        bvals /= 1000.0
    # Round to nearest 100 to force identical values within each shell,
    # which prevents FSL eddy's outlier detection from crashing with
    # "Mismatched DiffStatsVector" due to floating-point b-value spread.
    #bvals = np.round(bvals / 100.0) * 100.0
    #logger.info("b-values after standardisation: %s", np.unique(bvals))
    return bvals


def normalize_bvecs_unit(bvecs_path: str, out_path: Optional[str] = None) -> None:
    """Normalize bvec vectors to unit length, preserving zero vectors (b0).

    Operates on FSL-format bvecs (3 rows × N cols).  Each non-zero column is
    divided by its L2 norm so that dipy's gradient_table() accepts the result
    without raising a unit-vector tolerance error.

    Parameters
    ----------
    bvecs_path:
        Path to the bvecs file to normalize.
    out_path:
        Destination path.  When ``None`` the file is updated in-place.
    """
    bvecs = np.loadtxt(bvecs_path)

    # Work with (3, N) layout; transpose if loaded as (N, 3)
    if bvecs.ndim == 2 and bvecs.shape[1] == 3 and bvecs.shape[0] != 3:
        bvecs = bvecs.T

    if bvecs.ndim == 2 and bvecs.shape[0] == 3:
        norms = np.linalg.norm(bvecs, axis=0)
        norms[norms == 0] = 1.0  # leave b0 directions as zero vectors
        bvecs = bvecs / norms
    else:
        logger.warning(
            "normalize_bvecs_unit: unexpected bvecs shape %s; skipping.", bvecs.shape
        )
        return

    dest = out_path if out_path is not None else bvecs_path
    np.savetxt(dest, bvecs, fmt="%.15f")
    logger.debug("bvecs unit-normalized → %s", dest)


def normalize_bvec_bval_to_fsl(bvecs_path: str, bvals_path: str) -> None:
    """Normalize bvec and bval files to FSL/dcm2niix format in-place.

    FSL convention: bvecs → 3 rows × N cols (one row per axis X, Y, Z);
                    bvals → single row of N values.
    TRACULA convention: bvecs → N rows × 3 cols (one row per direction);
                        bvals → N rows of 1 value each.

    Detects format from array shape and rewrites only when a conversion is
    needed, so calling this on already-correct files is safe (idempotent).
    The N == 3 case (square 3×3 matrix) is logged as a warning and left
    unchanged because the format cannot be determined without extra context.

    After format transposing, vectors are also unit-normalized so downstream
    tools (dipy gradient_table, AMICO) never raise a unit-vector tolerance
    error regardless of scanner or conversion tool used.
    """
    bvecs = np.loadtxt(bvecs_path)

    if bvecs.ndim == 2 and bvecs.shape[0] == 3 and bvecs.shape[1] == 3:
        logger.warning(
            "bvecs is 3×3 - format is ambiguous (TRACULA vs FSL); "
            "skipping format normalization."
        )
    elif bvecs.ndim == 2 and bvecs.shape[1] == 3 and bvecs.shape[0] != 3:
        # TRACULA (N, 3) → FSL (3, N)
        logger.info(
            "bvecs detected in TRACULA column-per-direction format %s; "
            "transposing to FSL row-per-axis format.",
            bvecs.shape,
        )
        bvecs = bvecs.T
        np.savetxt(bvecs_path, bvecs, fmt="%.8f")
    elif bvecs.ndim == 2 and bvecs.shape[0] == 3:
        logger.debug("bvecs already in FSL format %s; no change.", bvecs.shape)
    else:
        logger.warning(
            "bvecs has unexpected shape %s; skipping normalization.", bvecs.shape
        )

    # Unit-normalize magnitudes so dipy/AMICO accept the vectors
    normalize_bvecs_unit(bvecs_path)

    # np.loadtxt returns a 1-D array for both FSL (one row) and TRACULA (one
    # column) bval files.  Standardise units and round before rewriting as a
    # single row to guarantee FSL format.
    bvals = np.loadtxt(bvals_path)
    bvals = _standardize_bvals(bvals)
    np.savetxt(bvals_path, bvals.reshape(1, -1), fmt="%g")
    logger.debug("bvals written as single-row FSL format.")


def normalize_bvec_bval_to_tracula(bvecs_path: str, bvals_path: str) -> None:
    """Normalize bvec and bval files to TRACULA format in-place.

    The inverse of normalize_bvec_bval_to_fsl(): TRACULA's trac-preproc reads
    bvecs as one row per DWI volume (N rows x 3 cols) and bvals as N rows of
    one value each -- the opposite of the FSL/dcm2niix row-per-axis
    convention (3 rows x N cols / 1 row x N cols) that dti/noddi need and
    that normalize_bvec_bval_to_fsl() produces.

    Feeding trac-preproc an FSL-format bvecs file makes it read the 3
    axis-rows as if they were 3 gradient vectors, regardless of the real
    direction count -- it fails with a hard "Found N b-values but 3
    gradient vectors" error for any subject with more than 3 real
    directions (confirmed against production logs and reproduced with
    synthetic data).

    Detects format from array shape and rewrites only when a conversion is
    needed, so calling this on already-TRACULA-format files is safe
    (idempotent). The N == 3 case (square 3x3 matrix) is logged as a
    warning and left unchanged, same ambiguity caveat as
    normalize_bvec_bval_to_fsl().
    """
    bvecs = np.loadtxt(bvecs_path)

    if bvecs.ndim == 2 and bvecs.shape[0] == 3 and bvecs.shape[1] == 3:
        logger.warning(
            "bvecs is 3×3 - format is ambiguous (TRACULA vs FSL); "
            "skipping format normalization."
        )
    elif bvecs.ndim == 2 and bvecs.shape[0] == 3 and bvecs.shape[1] != 3:
        # FSL (3, N) → TRACULA (N, 3)
        logger.info(
            "bvecs detected in FSL row-per-axis format %s; "
            "transposing to TRACULA column-per-direction format.",
            bvecs.shape,
        )
        bvecs = bvecs.T
        np.savetxt(bvecs_path, bvecs, fmt="%.8f")
    elif bvecs.ndim == 2 and bvecs.shape[1] == 3:
        logger.debug("bvecs already in TRACULA format %s; no change.", bvecs.shape)
    else:
        logger.warning(
            "bvecs has unexpected shape %s; skipping normalization.", bvecs.shape
        )

    # np.loadtxt returns a 1-D array for both FSL (one row) and TRACULA (one
    # column) bval files. Standardise units and rewrite as one value per
    # row to guarantee TRACULA format.
    bvals = np.loadtxt(bvals_path)
    bvals = _standardize_bvals(bvals)
    np.savetxt(bvals_path, bvals.reshape(-1, 1), fmt="%g")
    logger.debug("bvals written as one-value-per-row TRACULA format.")


def find_mismatched_bvec_volumes(bvecs_path: str, bvals_path: str) -> np.ndarray:
    """Return indices of volumes with a zero-length bvec but a nonzero bval.

    dipy's gradient_table() only enforces its unit-vector check on rows
    above the b0 threshold, so a direction that's genuinely missing
    (bvec == [0, 0, 0] - typically a corrupted or dropped gradient rather
    than a real b0) but still carries a nonzero bval raises a hard
    ValueError there instead of being caught at the source. A zero-length
    vector carries no directional information, so these volumes cannot be
    used for tensor/NODDI fitting regardless of their bval.

    Expects bvecs_path/bvals_path already in FSL format (3×N / 1×N), i.e.
    called after normalize_bvec_bval_to_fsl().
    """
    bvecs = np.loadtxt(bvecs_path)  # (3, N)
    bvals = np.loadtxt(bvals_path).reshape(-1)  # (N,)
    if bvecs.ndim != 2 or bvecs.shape[0] != 3 or bvecs.shape[1] != bvals.shape[0]:
        return np.array([], dtype=int)
    norms = np.linalg.norm(bvecs, axis=0)
    return np.where((norms < 1e-6) & (bvals > 0))[0]


def drop_volumes(dwi_path: str, bvecs_path: str, bvals_path: str, indices) -> None:
    """Remove the given volume indices from a DWI NIfTI and its bvecs/bvals, in place.

    Used to discard volumes whose gradient direction is missing/corrupted
    (see find_mismatched_bvec_volumes) - unlike zeroing their bval to treat
    them as b0, this keeps the b0/S0 reference free of diffusion-attenuated
    signal from a volume whose true direction is unknown, since dtifit and
    AMICO both average all b0 volumes together to build that reference.
    """
    import nibabel as nib

    indices = np.asarray(indices, dtype=int)
    if indices.size == 0:
        return

    img = nib.load(dwi_path)
    data = img.get_fdata()
    keep = np.ones(data.shape[-1], dtype=bool)
    keep[indices] = False

    nib.save(nib.Nifti1Image(data[..., keep], img.affine, img.header), dwi_path)

    bvecs = np.loadtxt(bvecs_path)[:, keep]
    bvals = np.loadtxt(bvals_path).reshape(-1)[keep]
    np.savetxt(bvecs_path, bvecs, fmt="%.15f")
    np.savetxt(bvals_path, bvals.reshape(1, -1), fmt="%g")

    logger.warning(
        "Dropped %d volume(s) from %s (indices %s) - zero-length bvec with "
        "nonzero bval, direction missing/corrupted. %d volumes remain.",
        indices.size, dwi_path, indices.tolist(), int(keep.sum()),
    )


def _record_input_failure(session_dir: Path, error_msg: str, source: str = "") -> None:
    """Atomically write <session_dir>/input_checks/dwi_failed.txt.

    session_dir must be the session's permanent output folder (the one
    process-all looks in), not node scratch, which is deleted when the job ends.
    Written via a temp file + os.replace so a reader sees either no marker or a
    complete one; concurrent writers (GIF and FreeSurfer TRACULA checking the
    same DWI) write the same content, so the last one landing is harmless.
    """
    try:
        checks_dir = Path(session_dir) / "input_checks"
        checks_dir.mkdir(parents=True, exist_ok=True)
        marker_path = checks_dir / "dwi_failed.txt"
        tmp_marker = checks_dir / f".tmp_{os.getpid()}_dwi_failed.txt"
        with open(tmp_marker, "w") as f:
            f.write(f"{error_msg.strip()}\n")
            if source:
                f.write(f"input: {source}\n")
        os.replace(tmp_marker, marker_path)
        logger.warning("Recorded DWI input failure marker at %s", marker_path)
    except OSError as exc:
        # The marker is an optimisation; failing to write it must not hide the real error.
        logger.warning("Could not record input check marker in %s: %s", session_dir, exc)


def check_input_failure_marker(
    session_dir: Path | str,
    raw_input_paths: Optional[list[Path | str]] = None,
    retry: bool = False,
    subject: str = "",
) -> bool:
    """Check if session should be skipped due to a previous DWI input check failure.

    Inspects `<session_dir>/input_checks/dwi_failed.txt`.

    Parameters
    ----------
    session_dir : Path or str
        Base output directory for the session/subject (contains `input_checks/`).
    raw_input_paths : list of (Path or str), optional
        List of raw DWI input candidate paths. If any existing file is newer than
        the failure marker, the marker is treated as stale, removed, and the session
        is retried.
    retry : bool, default False
        If True (e.g. `--retry-input-failures`), existing failure markers are removed
        and the session is not skipped.
    subject : str, optional
        Subject identifier for logging messages (defaults to directory name if empty).

    Returns
    -------
    bool
        True if the session should be skipped, False if processing should proceed.
    """
    session_path = Path(session_dir)
    subj_marker = session_path / "input_checks" / "dwi_failed.txt"
    if not subj_marker.is_file():
        return False

    subj_label = subject or session_path.name

    if retry:
        logger.info(
            "Subject '%s': retrying despite previous failure marker (--retry-input-failures)",
            subj_label,
        )
        try:
            subj_marker.unlink()
        except OSError:
            pass
        return False

    # Check if raw input files are newer than the marker
    try:
        marker_mtime = subj_marker.stat().st_mtime
        if raw_input_paths:
            for cand in raw_input_paths:
                p = Path(cand)
                if p.exists() and p.stat().st_mtime > marker_mtime:
                    logger.info(
                        "Subject '%s': DWI input is newer than failure marker; removing marker and retrying.",
                        subj_label,
                    )
                    try:
                        subj_marker.unlink()
                    except OSError:
                        pass
                    return False
    except Exception:
        pass

    reason = ""
    try:
        reason = subj_marker.read_text().strip()
    except Exception:
        pass

    logger.warning(
        "Skipping subject '%s': previously failed DWI input checks (%s). "
        "Use --retry-input-failures to reprocess.",
        subj_label,
        reason or "dwi_failed.txt",
    )
    return True


def prepare_dwi_input(dwi_input: str,
                      bvecs: Optional[str] = None,
                      bvals: Optional[str] = None,
                      output_dir: Optional[str] = None,
                      verbose: bool = False,
                      dcm2niix_path: Optional[str] = None,
                      session_dir: Optional[str] = None) -> Dict[str, Optional[str]]:
    """Prepare DWI, bvec and bval files from a variety of inputs.

    session_dir: the session's permanent output folder. When given and the raw
    DWI is unusable (InputDataError), a failure marker is written there so
    process-all skips the session next time. Other errors never leave a marker.

    Returns a dict: { 'dwi': <path>, 'bvecs': <path>, 'bvals': <path>, 'output_dir': <path> }

    Raises FileNotFoundError or RuntimeError on failure.
    """
    try:
        return _prepare_dwi_input_impl(dwi_input, bvecs, bvals, output_dir, verbose, dcm2niix_path=dcm2niix_path)
    except InputDataError as exc:
        if session_dir:
            _record_input_failure(Path(session_dir), str(exc), source=str(dwi_input))
        raise


def _prepare_dwi_input_impl(dwi_input: str,
                            bvecs: Optional[str] = None,
                            bvals: Optional[str] = None,
                            output_dir: Optional[str] = None,
                            verbose: bool = False,
                            dcm2niix_path: Optional[str] = None) -> Dict[str, Optional[str]]:
    input_p = Path(dwi_input)
    output_dir_path = Path(output_dir) if output_dir else input_p.parent

    # Only DICOM inputs (.zip archive or directory, both handled further
    # below) ever need dcm2niix -- .nii/.nii.gz inputs never call it.
    # Converter binaries are discovered portably in _convert_dicom_with_fallback.

    if not input_p.exists():
        raise FileNotFoundError(f"DWI input not found: {dwi_input}")

    # Case: .nii.gz -> copy to temp workspace
    if input_p.is_file() and input_p.name.endswith('.nii.gz'):
        logger.info("DWI input is .nii.gz file")

        temp_output_dir = Path(tempfile.mkdtemp())

        target = temp_output_dir / input_p.name
        shutil.copy2(input_p, target)
        dwi = target

        # Copy any JSON sidecar that lives next to the original input so the
        # later glob (which looks at the temp workspace) picks it up.
        for json_src in input_p.parent.glob("*.json"):
            shutil.copy2(json_src, temp_output_dir / json_src.name)

        # discover sidecars if not provided -- search the ORIGINAL input's
        # directory, not the temp workspace: only the .nii.gz and any .json
        # get copied there, never the bvec/bval sidecars themselves.
        if not bvecs or not bvals:
            bvec_found, bval_found = _find_sidecars(input_p.parent, stem=input_p.name[:-len(".nii.gz")])
            if not bvecs and bvec_found:
                bvecs = str(bvec_found)
            if not bvals and bval_found:
                bvals = str(bval_found)

    # Case: .nii -> compress into .nii.gz
    elif input_p.is_file() and input_p.suffix == '.nii':
        logger.info("DWI input is .nii file")

        temp_output_dir = Path(tempfile.mkdtemp())
        target = temp_output_dir / f"{input_p.stem}.nii.gz"

        with open(input_p, 'rb') as f_in, gzip.open(target, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        dwi = target

        for json_src in input_p.parent.glob("*.json"):
            shutil.copy2(json_src, temp_output_dir / json_src.name)

        if not bvecs or not bvals:
            bvec_found, bval_found = _find_sidecars(input_p.parent, stem=input_p.stem)
            if not bvecs and bvec_found:
                bvecs = str(bvec_found)
            if not bvals and bval_found:
                bvals = str(bval_found)


    # Case: zip -> extract and run converter with fallback on extracted folder
    elif input_p.is_file() and input_p.suffix == '.zip':
        logger.info("DWI input is .zip archive")
        
        # Create a temporary directory for extraction
        extract_dir = Path(tempfile.mkdtemp())

        #logger.info(f"Extracting zip archive {input_p} to {extract_dir}")
        with zipfile.ZipFile(input_p, 'r') as zf:
            zf.extractall(extract_dir)

        temp_output_dir = Path(tempfile.mkdtemp())

        # run converter with fallback on extracted folder
        dwi_found, bvec_cand, bval_cand = _convert_dicom_with_fallback(
            extract_dir, temp_output_dir, dcm2niix_path, verbose=verbose
        )
        
        # Delete extracted folder after use
        shutil.rmtree(extract_dir)
        dwi = dwi_found

        if not bvecs and bvec_cand:
            bvecs = str(bvec_cand)
        if not bvals and bval_cand:
            bvals = str(bval_cand)


    # Case: directory -> assume DICOMs and run converter with fallback
    elif input_p.is_dir():
        logger.info("DWI input is a directory (assuming DICOMs)")
        temp_output_dir = Path(tempfile.mkdtemp())
        dwi_found, bvec_cand, bval_cand = _convert_dicom_with_fallback(
            input_p, temp_output_dir, dcm2niix_path, verbose=verbose
        )
        dwi = dwi_found

        if not bvecs and bvec_cand:
            bvecs = str(bvec_cand)
        if not bvals and bval_cand:
            bvals = str(bval_cand)

    else:
        raise FileNotFoundError(f"Unsupported DWI input: {dwi_input}")

    # Final checks: ensure we have dwi, bvecs, bvals
    if not dwi or not Path(dwi).exists():
        raise FileNotFoundError(f"Prepared DWI not found after processing: {dwi}")
    if not bvecs or not Path(bvecs).exists():
        raise FileNotFoundError(f"bvecs file not found after processing. Expected at: {bvecs}")
    if not bvals or not Path(bvals).exists():
        raise FileNotFoundError(f"bvals file not found after processing. Expected at: {bvals}")

    # Rename outputs to standard names in output_dir
    final_dwi = output_dir_path / f"data.nii.gz"
    final_bvecs = output_dir_path / f"bvecs"
    final_bvals = output_dir_path / f"bvals"

    # Find json file and copy if exists
    json_files = list(Path(dwi).parent.glob("*.json"))
    if json_files:
        json_file = json_files[0]
        final_json = output_dir_path / "metadata.json"
        shutil.copy2(src=str(json_file), dst=final_json)
        logger.info(f"Copied associated JSON metadata: {json_file} -> {final_json}")

    logger.info(f"Copied prepared files dwi: {dwi} -> {final_dwi}, bvecs: {bvecs} -> {final_bvecs}, bvals: {bvals} -> {final_bvals}")

    shutil.copy2(src=str(dwi), dst=final_dwi)
    shutil.copy2(src=str(bvecs), dst=final_bvecs)
    shutil.copy2(src=str(bvals), dst=final_bvals)

    # The scratch workspace (temp_output_dir, set in every branch above)
    # served its purpose once the prepared files are copied to their
    # permanent location above -- without this cleanup, every call leaves a
    # full copy of the (often ~1GB) DWI volume behind in node-local /tmp.
    # z-score's healthy-cohort loop calls this once per healthy subject
    # within the same job, so an uncleaned temp dir here exhausted local
    # disk after ~30-odd subjects ("No space left on device") in practice.
    shutil.rmtree(temp_output_dir, ignore_errors=True)

    # Normalise bvecs/bvals to FSL format (3×N and 1×N) so that downstream
    # tools (dipy, dmipy) receive a consistent layout regardless of whether
    # the source files came from dcm2niix or TRACULA.
    normalize_bvec_bval_to_fsl(str(final_bvecs), str(final_bvals))

    # Validate final DWI dimensions and b-value consistency before dropping volumes:
    import nibabel as nib
    final_img = nib.load(str(final_dwi))
    if final_img.ndim < 4:
        raise InputDataError(f"Prepared DWI is not 4D: shape is {final_img.shape}.")
    if int(np.sum(np.loadtxt(str(final_bvals)).reshape(-1) > 50)) < MIN_DW_VOLUMES:
        raise InputDataError(
            f"Prepared DWI has invalid volume count: shape is {final_img.shape} "
            f"(minimum {MIN_DW_VOLUMES} diffusion-weighted volumes required for tensor reconstruction)."
        )

    bvals_arr = np.loadtxt(str(final_bvals)).reshape(-1)
    if len(bvals_arr) != final_img.shape[3]:
        raise InputDataError(
            f"Volume count mismatch: DWI image has {final_img.shape[3]} volumes "
            f"but {len(bvals_arr)} b-values were found."
        )

    # Drop any volume whose gradient direction is missing/corrupted (zero
    # bvec paired with a nonzero bval) - see drop_volumes() docstring for
    # why these are removed outright rather than relabelled as b0.
    mismatched = find_mismatched_bvec_volumes(str(final_bvecs), str(final_bvals))
    if mismatched.size > 0:
        drop_volumes(str(final_dwi), str(final_bvecs), str(final_bvals), mismatched)
        # Re-verify sufficient volumes remain after dropping corrupted ones
        final_img = nib.load(str(final_dwi))
        remaining_bvals = np.loadtxt(str(final_bvals)).reshape(-1)
        if int(np.sum(remaining_bvals > 50)) < MIN_DW_VOLUMES:
            raise InputDataError(
                f"Prepared DWI has insufficient valid volumes after dropping corrupted directions: "
                f"{final_img.shape[3]} remaining (minimum {MIN_DW_VOLUMES} diffusion-weighted required)."
            )

    logger.info(f"Final prepared files: DWI: {final_dwi}, bvecs: {final_bvecs}, bvals: {final_bvals}")

    return {
        'dwi': str(final_dwi.absolute()),
        'bvecs': str(final_bvecs.absolute()),
        'bvals': str(final_bvals.absolute()),
        'output_dir': str(output_dir_path.absolute())
    }

def main(argv=None):
    p = argparse.ArgumentParser(description='Prepare DWI input. Accepts .nii.gz, .nii, .zip (DICOM), or DICOM directory.')
    p.add_argument('--dwi', required=True, help='DWI NIfTI (.nii or .nii.gz)')
    p.add_argument('--bvecs', required=False, help='bvecs file')
    p.add_argument('--bvals', required=False, help='bvals file')
    p.add_argument('--outdir', required=True, help='Output directory')
    p.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    p.add_argument('--session-dir', required=False,
                   help="Session's permanent output folder; receives input_checks/dwi_failed.txt if the DWI is unusable")

    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        stream=sys.stdout)

    logging.info(f"Preparing DWI input from {args.dwi} with bvecs {args.bvecs} and bvals {args.bvals}")

    paths_dict = prepare_dwi_input(args.dwi, args.bvecs, args.bvals, args.outdir, session_dir=args.session_dir)

    return paths_dict["dwi"], paths_dict["bvecs"], paths_dict["bvals"], paths_dict["output_dir"]

if __name__ == '__main__':
    main()