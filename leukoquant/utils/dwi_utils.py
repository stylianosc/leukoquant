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


def _run_dcm2niix(dcm2niix_exec: str, src_dir: Path, out_dir: Path, verbose: bool = False) -> None:
    exec_path = Path(dcm2niix_exec)
    if not exec_path.exists():
        raise FileNotFoundError(f"dcm2niix executable not found at: {dcm2niix_exec}")
    if not exec_path.is_file():
        raise RuntimeError(f"dcm2niix path is not a file: {dcm2niix_exec}")
    if not os.access(str(exec_path), os.X_OK):
        raise RuntimeError(
            f"dcm2niix at {dcm2niix_exec} is not executable.\n"
            "Make it executable (chmod +x ...) or install dcm2niix in your PATH and omit dcm2niix_path."
        )

    cmd = [str(exec_path), '-z', 'y', '-b', 'y', '-o', str(out_dir), str(src_dir)]
    if verbose:
        print(f"Running dcm2niix: {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"dcm2niix failed (rc={proc.returncode}). stderr: {proc.stderr}")


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


def prepare_dwi_input(dwi_input: str,
                      bvecs: Optional[str] = None,
                      bvals: Optional[str] = None,
                      output_dir: Optional[str] = None,
                      verbose: bool = False) -> Dict[str, Optional[str]]:
    """Prepare DWI, bvec and bval files from a variety of inputs.

    Returns a dict: { 'dwi': <path>, 'bvecs': <path>, 'bvals': <path>, 'tempdir': <path or None> }

    Raises FileNotFoundError or RuntimeError on failure.
    """
    input_p = Path(dwi_input)
    output_dir_path = Path(output_dir) if output_dir else input_p.parent

    # Only DICOM inputs (.zip archive or directory, both handled further
    # below) ever need dcm2niix -- .nii/.nii.gz inputs never call it. Don't
    # require the binary to exist up front for every input type; leave the
    # existence/executable check to _run_dcm2niix itself, which only runs
    # when a DICOM branch is actually taken.
    dcm2niix_path = str(Path("/leukoquant/leukoquant/external/noddi/dcm2niix"))

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


    # Case: zip -> extract and run dcm2niix on extracted folder
    elif input_p.is_file() and input_p.suffix == '.zip':
        logger.info("DWI input is .zip archive")
        
        # Create a temporary directory for extraction
        extract_dir = Path(tempfile.mkdtemp())

        #logger.info(f"Extracting zip archive {input_p} to {extract_dir}")
        with zipfile.ZipFile(input_p, 'r') as zf:
            zf.extractall(extract_dir)

        temp_output_dir = Path(tempfile.mkdtemp())

        # run dcm2niix on extracted folder
        _run_dcm2niix(dcm2niix_path, extract_dir, temp_output_dir, verbose=verbose)
        
        # Delete extracted folder after use
        shutil.rmtree(extract_dir)
        
        # find first nii/nii.gz and sidecars
        niis = list(temp_output_dir.glob('*.nii.gz'))
        #logger.info(f"Files produced by dcm2niix: {list(tempdir.glob('*'))}")
        if not niis:
            raise FileNotFoundError(f"dcm2niix produced no NIfTI files from {input_p}")
        dwi = niis[0]

        bvec_candidates = list(temp_output_dir.glob(f"*.bvec"))
        bval_candidates = list(temp_output_dir.glob(f"*.bval"))

        if not bvecs and bvec_candidates:
            bvecs = str(bvec_candidates[0])
        if not bvals and bval_candidates:
            bvals = str(bval_candidates[0])


    # Case: directory -> assume DICOMs and run dcm2niix
    elif input_p.is_dir():
        logger.info("DWI input is a directory (assuming DICOMs)")
        temp_output_dir = Path(tempfile.mkdtemp())
        _run_dcm2niix(dcm2niix_path, input_p, temp_output_dir, verbose=verbose)

        niis = list(temp_output_dir.glob('*.nii.gz'))
        if not niis:
            raise FileNotFoundError(f"dcm2niix produced no NIfTI files from DICOM folder {input_p}")
        dwi = niis[0]

        bvec_candidates = list(temp_output_dir.glob(f"*.bvec"))
        bval_candidates = list(temp_output_dir.glob(f"*.bval"))

        if not bvecs and bvec_candidates:
            bvecs = str(bvec_candidates[0])
        if not bvals and bval_candidates:
            bvals = str(bval_candidates[0])

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

    # Drop any volume whose gradient direction is missing/corrupted (zero
    # bvec paired with a nonzero bval) - see drop_volumes() docstring for
    # why these are removed outright rather than relabelled as b0.
    mismatched = find_mismatched_bvec_volumes(str(final_bvecs), str(final_bvals))
    if mismatched.size > 0:
        drop_volumes(str(final_dwi), str(final_bvecs), str(final_bvals), mismatched)

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

    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        stream=sys.stdout)

    logging.info(f"Preparing DWI input from {args.dwi} with bvecs {args.bvecs} and bvals {args.bvals}")

    paths_dict = prepare_dwi_input(args.dwi, args.bvecs, args.bvals, args.outdir)

    return paths_dict["dwi"], paths_dict["bvecs"], paths_dict["bvals"], paths_dict["output_dir"]

if __name__ == '__main__':
    main()