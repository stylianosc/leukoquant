#!/usr/bin/env python3
"""
Simple NODDI-Watson wrapper based on the example notebook.

Usage:
  python run_noddi.py --dwi DWI.nii.gz --bvecs bvecs --bvals bvals --outdir ./noddi_out [--mask mask.nii.gz]

This script attempts to use dmipy (preferred) and dipy for scheme creation.
It fits a NODDI-Watson model and writes three volumes to `outdir`:
  - `odi.nii.gz`
  - `fiso.nii.gz`  (CSF / isotropic volume fraction)
  - `ficvf.nii.gz` (intra-cellular volume fraction estimate)

Notes:
 - Fitting NODDI to a full brain can be slow. Provide a mask to limit voxels.
 - This script is defensive: it tries to build a dmipy acquisition scheme and
   will raise a helpful message if required APIs are missing.
"""

from pathlib import Path
import argparse
import logging
import sys
import os

import numpy as np
import nibabel as nib
import warnings

# Allow imports of leukoquant when run from containers where the module is
# bind-mounted at /leukoquant but not installed in site-packages.
if "/leukoquant" not in sys.path and os.path.isdir("/leukoquant"):
    sys.path.insert(0, "/leukoquant")

from leukoquant.utils.dwi_utils import normalize_bvecs_unit

logger = logging.getLogger(__name__)
# Set root logger to WARNING by default (only show errors/warnings unless --verbose)
logging.getLogger().setLevel(logging.WARNING)

def load_bvecs_bvals(bvecs_path, bvals_path):
    bvecs = np.loadtxt(str(bvecs_path))
    bvals = np.loadtxt(str(bvals_path))
    # Ensure shape (N,3)
    if bvecs.ndim == 1:
        bvecs = bvecs[np.newaxis, :]
    if bvecs.shape[0] == 3 and bvecs.shape[1] != 3:
        bvecs = bvecs.T

    return bvecs, bvals


def create_acquisition_scheme(bvals, bvecs):
    """Create a dmipy-compatible acquisition scheme.

    Returns a scheme object usable by dmipy MultiCompartmentModel.fit
    or raises RuntimeError with guidance.
    """
    # Preferred: dmipy helper
    from dmipy.core.acquisition_scheme import acquisition_scheme_from_bvalues

    # Check if bvalues in s/mm^2 or s/m^2
    if np.max(bvals) < 1e7:
        logging.info("Detected b-values in s/mm^2; converting to s/m^2")
        bvals = bvals * 1e6  # convert to s/m^2

    try:
        return acquisition_scheme_from_bvalues(bvals, bvecs)
    except Exception as exc:
        raise RuntimeError(
            "Failed to build dmipy acquisition scheme. Ensure dmipy is up to date."
        ) from exc



def fit_noddi_and_write(dwi_path, bvecs_path, bvals_path, outdir, mask_path=None, nprocs=1):
    from numba.core.errors import NumbaPendingDeprecationWarning

    warnings.filterwarnings("ignore", category=NumbaPendingDeprecationWarning)
    
    logging.info("Loading DWI data...")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    img = nib.load(str(dwi_path))
    data = img.get_fdata()
    affine = img.affine

    bvecs, bvals = load_bvecs_bvals(bvecs_path, bvals_path)

    logging.info(f"Loaded bvals ({bvals.shape}) and bvecs ({bvecs.shape})")

    scheme = create_acquisition_scheme(bvals, bvecs)
    logging.info(f"Unique B0 Indices: {scheme.unique_b0_indices}")
    logging.info(f"Unique DWI indices: {scheme.unique_dwi_indices}")
    logging.info(f"Bvalues: {scheme.bvalues}")
    logging.info(f"Shell Indices: {scheme.shell_indices}")
    logging.info(f"B0 Mask shape: {scheme.b0_mask}")

    mask = None
    if mask_path is not None:
        mask = nib.load(str(mask_path)).get_fdata().astype(bool)
    else:
        b0_mask = scheme.b0_mask
        b0_indices = np.where(b0_mask)[0].tolist()
        if len(b0_indices) > 0:
            logging.info("No mask provided; using non-zero voxels in B0 images as mask.")
            first_b0 = b0_indices[0]
            mask = data[..., first_b0] > 0
        else:
            logging.info("No mask provided and no B0 images found; fitting on all voxels.")
            mask = data[..., 0] > 0  # fallback: non-zero in first volume
    
    logging.info(f"Fitting NODDI model on {np.sum(mask)} voxels")

    # Build the model (mirror of the notebook)
    try:
        from dmipy.signal_models import cylinder_models, gaussian_models
        from dmipy.distributions.distribute_models import SD1WatsonDistributed
        from dmipy.core.modeling_framework import MultiCompartmentModel
    except Exception as exc:
        raise RuntimeError("dmipy is required for this script. Install dmipy and retry.") from exc

    logging.info("Fitting NODDI-Watson model...")
    ball = gaussian_models.G1Ball()
    stick = cylinder_models.C1Stick()
    zeppelin = gaussian_models.G2Zeppelin()
    watson_dispersed_bundle = SD1WatsonDistributed(models=[stick, zeppelin])

    watson_dispersed_bundle.set_tortuous_parameter('G2Zeppelin_1_lambda_perp','C1Stick_1_lambda_par','partial_volume_0')
    watson_dispersed_bundle.set_equal_parameter('G2Zeppelin_1_lambda_par', 'C1Stick_1_lambda_par')
    watson_dispersed_bundle.set_fixed_parameter('G2Zeppelin_1_lambda_par', 1.7e-9)

    NODDI_mod = MultiCompartmentModel(models=[ball, watson_dispersed_bundle])
    NODDI_mod.set_fixed_parameter('G1Ball_1_lambda_iso', 3e-9)

    # Fit model. This may take long for full volumes.
    fitted = NODDI_mod.fit(scheme, data, mask=mask, use_parallel_processing=True)

    fitted_parameters = fitted.fitted_parameters

    logging.info(f" NODDI parameters: {list(fitted_parameters.keys())}")

    # Heuristics to extract desired outputs
    #  - fiso: the isotropic / CSF partial volume (usually 'partial_volume_0')
    #  - ficvf: intra-axonal fraction -> partial_volume_watson * SD1WatsonDistributed_..._partial_volume_0
    #  - odi: search for a key containing 'odi'

    csf_vf = fitted_parameters.get('partial_volume_0', None)
    if csf_vf is None:
        raise RuntimeError("Could not find CSF partial volume in fitted parameters.")
    
    extra_vf = fitted_parameters.get("partial_volume_1", None)
    if extra_vf is None:
        raise RuntimeError("Could not find extra-cellular partial volume in fitted parameters.")
    
    norm_vf_stick_watson = fitted_parameters.get("SD1WatsonDistributed_1_partial_volume_0", None)
    
    if norm_vf_stick_watson is None:
        raise RuntimeError("Could not find required partial volume in Watson distributed model.")

    vf_intra = extra_vf * norm_vf_stick_watson

    odi = fitted_parameters.get('SD1WatsonDistributed_1_SD1Watson_1_odi', None)
    # Find keys containing 'odi' if not found directly
    if odi is None:
        odi_keys = [k for k in fitted_parameters.keys() if 'odi' in k.lower()]
        if len(odi_keys) == 1:
            odi = fitted_parameters[odi_keys[0]]
        elif len(odi_keys) > 1:
            logging.warning(f"Multiple ODI keys found: {odi_keys}. Using the first one.")
            odi = fitted_parameters[odi_keys[0]]
        else:
            raise RuntimeError("Could not find ODI in fitted parameters.")
    
    # Prepare output arrays
    # FICVF: intra-cellular volume fraction
    ficvf_arr = np.zeros(data.shape[:3], dtype=np.float32)
    # FISO: isotropic volume fraction
    fiso_arr = np.zeros(data.shape[:3], dtype=np.float32)
    # ODI
    odi_arr = np.zeros(data.shape[:3], dtype=np.float32)

    # Fill output arrays
    if mask is not None:
        ficvf_arr[mask] = vf_intra[mask]
        fiso_arr[mask] = csf_vf[mask]
        odi_arr[mask] = odi[mask]

    else:
        ficvf_arr = vf_intra
        fiso_arr = csf_vf
        odi_arr = odi

    # Save results
    def save_array(arr, fname):
        if arr is None:
            logging.warning(f"No data for {fname}; skipping write.")
            return
        img_out = nib.Nifti1Image(arr.astype(np.float32), affine)
        file_path = outdir / fname
        nib.save(img_out, str(file_path))
        logging.info(f"Saved file: {file_path}")

    # Check number of unique DWI values and if there's only one 
    # then it's single-shell and we store only ODI.
    # Otherwise, it's multi-shell and we store all three.

    unique_dwi_indices = scheme.unique_dwi_indices.tolist()
    if len(unique_dwi_indices) == 1:
        logging.info("Single-shell DWI detected; writing only ODI output.")
        save_array(odi_arr, 'odi.nii.gz')
    else:
        logging.info("Multi-shell DWI detected; writing ODI, FISO, and FICVF outputs.")
        save_array(odi_arr, 'odi.nii.gz')
        save_array(fiso_arr, 'fiso.nii.gz')
        save_array(ficvf_arr, 'ficvf.nii.gz')
    

    logging.info(f"Wrote outputs to {outdir}")


def run_noddi_amico(dwi_path, bvecs_path=None, bvals_path=None, outdir: str = "", mask_path=None, scratch_dir=None):
    """Run NODDI fitting using AMICO framework.

    All computation and intermediate files are kept in scratch (when provided).
    Only the final output maps (odi, fwf, ndi) are copied to outdir at the end,
    keeping the output directory clean throughout the run.

    Intermediate files (scheme.txt, binarized_mask.nii.gz, dir.nii.gz) are
    written to scratch so the shell's keep_intermediate logic can retrieve them
    from there; when keep_intermediate is False the scratch trap discards them.
    """
    import amico
    import tempfile
    import shutil

    logging.info("Running NODDI fitting using AMICO...")
    amico.set_verbose(2)

    amico.setup()

    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)

    # All intermediate work lives in scratch when available; fall back to outdir
    # for standalone runs where no scratch directory is provided.
    scratch_path = Path(scratch_dir) if scratch_dir else outdir_path

    # Normalize bvecs to scratch to avoid dipy unit vector failures
    norm_bvecs_file = scratch_path / 'normalized_bvecs.bvec'
    normalize_bvecs_unit(bvecs_path, out_path=str(norm_bvecs_file))

    # scheme.txt and binarized_mask go to scratch so outdir stays clean.
    # The shell's keep_intermediate block reads them from scratch if needed.
    scheme_file = scratch_path / 'scheme.txt'
    amico.util.fsl2scheme(
        bvalsFilename=bvals_path,
        bvecsFilename=str(norm_bvecs_file),
        schemeFilename=str(scheme_file),
        bStep=100,
    )


    final_mask_path = mask_path
    if mask_path is not None:
        img = nib.load(mask_path)
        mask_data = (img.get_fdata() > 0).astype(np.uint8)
        binarized_mask_path = str(scratch_path / 'binarized_mask.nii.gz')
        nib.save(nib.Nifti1Image(mask_data, img.affine, img.header), binarized_mask_path)
        final_mask_path = binarized_mask_path

    # AMICO kernel generation and fitting run inside a private temp dir on scratch
    # to avoid NFS overhead and isolate per-job state.
    tmp_dir = Path(tempfile.mkdtemp(dir=scratch_path))
    os.environ["DIPY_HOME"] = str(tmp_dir.absolute())

    ae = amico.Evaluation(study_path=str(tmp_dir.absolute()), output_path=str(tmp_dir.absolute()))
    # b0_min_signal zeroes voxels whose b0 is <=1% of the mean in-brain b0 signal
    # *before* AMICO's normalization divides by it -- without this, a handful of
    # near-zero-signal voxels (edge/noise, not caught by the corrupted-bvec drop
    # in dwi_utils.py) blow up to Inf/NaN during the 1/norm_factor step and hard-crash
    # the whole fit. replace_bad_voxels is a safety net for any NaN/Inf that still
    # slips through from another source.
    ae.load_data(
        dwi_filename=dwi_path, scheme_filename=str(scheme_file), mask_filename=final_mask_path,
        b0_min_signal=0.01, replace_bad_voxels=0,
    )

    # AMICO scheme object attributes
    unique_bvals = sorted(set(ae.scheme.b.round(-2)))  # round to nearest 100 s/mm^2
    non_b0_shells = [b for b in unique_bvals if b > 0]

    logging.info(f"Unique b-values detected: {unique_bvals}")
    logging.info(f"DWI shells (non-b0): {non_b0_shells}")

    ae.set_model('NODDI')
    # regenerate=False: generate kernels only when absent (they always will be
    # in the fresh tmp_dir), avoiding AMICO's remove() call on partial files.
    ae.generate_kernels(regenerate=False)
    ae.load_kernels()
    ae.fit()

    # AMICO's save_results() tries to remove files but fails on directories.
    # The kernels/ subdirectory is created during generation, so remove it
    # before save_results() attempts cleanup to avoid IsADirectoryError.
    kernels_path = tmp_dir / 'kernels'
    if kernels_path.exists():
        shutil.rmtree(kernels_path, ignore_errors=True)

    ae.save_results()

    # Determine output vs intermediate maps based on shell count.
    # dir.nii.gz is an intermediate; it goes to scratch (not outdir) so the
    # shell's keep_intermediate block can optionally preserve it.
    if len(non_b0_shells) == 1:
        logging.info("Single-shell DWI detected; writing only ODI output.")
        output_maps      = ["ODI"]
        intermediate_maps = ["dir"]
    else:
        logging.info("Multi-shell DWI detected; writing ODI, FWF, and NDI outputs.")
        output_maps      = ["ODI", "FWF", "NDI"]
        intermediate_maps = ["dir"]

    # Copy only final output maps to outdir.
    for fname in output_maps:
        src_files = list(tmp_dir.glob(f'*{fname}*.nii.gz'))
        if src_files:
            dest_file = outdir_path / f'{fname.lower()}.nii.gz'
            shutil.copy2(src=str(src_files[0]), dst=str(dest_file))
            logging.info(f"Copied {src_files[0]} to {dest_file}")
        else:
            logging.warning(f"Output file for {fname} not found in results.")

    # Stage intermediate maps in scratch so the shell can handle keep_intermediate.
    for fname in intermediate_maps:
        src_files = list(tmp_dir.glob(f'*{fname}*.nii.gz'))
        if src_files:
            dest_file = scratch_path / f'{fname.lower()}.nii.gz'
            shutil.copy2(src=str(src_files[0]), dst=str(dest_file))
            logging.info(f"Staged intermediate {src_files[0]} to scratch at {dest_file}")
        else:
            logging.warning(f"Intermediate file for {fname} not found in results.")

    # Discard the inner temp dir; scratch_path itself is managed by the shell trap.
    shutil.rmtree(tmp_dir)



def main(argv=None):
    p = argparse.ArgumentParser(description='Fit NODDI-Watson and write ODI/FISO/FICVF')
    p.add_argument('--dwi', required=True, help='DWI NIfTI (.nii or .nii.gz)')
    p.add_argument('--bvecs', required=True, help='bvecs file')
    p.add_argument('--bvals', required=True, help='bvals file')
    p.add_argument('--outdir', required=True, help='Output directory')
    p.add_argument('--mask', required=False, help='Optional brain mask NIfTI')
    p.add_argument('--scratch', required=False, help='Scratch directory for AMICO temp files (defaults to outdir)')
    p.add_argument('--nprocs', type=int, default=1, help='Number of processes (if supported)')
    p.add_argument('--verbose', action='store_true')
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        stream=sys.stdout)

    logging.info(f"Running NODDI fit on {args.dwi} with bvecs {args.bvecs} and bvals {args.bvals}")

    run_noddi_amico(args.dwi, args.bvecs, args.bvals, args.outdir, mask_path=args.mask, scratch_dir=args.scratch)


if __name__ == '__main__':
    main()