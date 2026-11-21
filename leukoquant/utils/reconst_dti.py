import logging
import argparse
import os
import sys
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

# Allow imports of leukoquant when run from containers where the module is
# bind-mounted at /leukoquant but not installed in site-packages.
if "/leukoquant" not in sys.path and os.path.isdir("/leukoquant"):
    sys.path.insert(0, "/leukoquant")

from leukoquant.utils.dwi_utils import normalize_bvecs_unit

logger = logging.getLogger(__name__)
# Set root logger to WARNING by default (only show errors/warnings unless --verbose)
logging.getLogger().setLevel(logging.WARNING)

def compute_dti_metrics(nifti_fname,
                        bval_fname,
                        bvec_fname,
                        mask_fname=None,
                        out_dir=".",
                        save_intermediate=False,
                        ):
    """
    Load DWI, compute DTI (FA, MD, RGB, evecs) and save NIfTIs.
    Returns a dict with numpy arrays for fa, md, rgb, evecs and the mask.
    """
    import numpy as np
    import dipy.reconst.dti as dti
    from dipy.core.gradients import gradient_table
    from dipy.io.gradients import read_bvals_bvecs
    from dipy.io.image import load_nifti, save_nifti

    data, affine = load_nifti(nifti_fname)

    # Ensure bvecs are unit-length before dipy ingests them.  Some scanner
    # outputs have vectors that are close to but not exactly unit length,
    # which causes gradient_table() to raise a ValueError.
    normalize_bvecs_unit(bvec_fname)

    bvals, bvecs = read_bvals_bvecs(bval_fname, bvec_fname)
    gtab = gradient_table(bvals, bvecs=bvecs)

    # 3D brain mask
    if mask_fname is None:
        mask = np.ones(data.shape[:3], dtype=np.uint8)
        maskdata = data
    else:
        mask, _ = load_nifti(mask_fname)
        mask = (mask > 0).astype(np.uint8)
        maskdata = data * mask[..., np.newaxis]

    logger.info(f"Fitting Tensor model...")
    tenmodel = dti.TensorModel(gtab)
    tenfit = tenmodel.fit(maskdata)

    L1 = np.abs(tenfit.evals[..., 0])
    L2 = np.abs(tenfit.evals[..., 1])
    L3 = np.abs(tenfit.evals[..., 2])

    numerator_fa = np.sqrt((L1 - L2)**2 + (L2 - L3)**2 + (L3 - L1)**2)
    denominator_fa = np.sqrt(2 * (L1**2 + L2**2 + L3**2))

    FA = np.divide(numerator_fa, denominator_fa, 
                   out=np.zeros_like(numerator_fa), where=denominator_fa!=0)
    FA_INVERSE = 1 - FA
    RGB = dti.color_fa(FA, tenfit.evecs)
    MD = (L1 + L2 + L3) / 3.0
    RD = (L2 + L3) / 2.0
    AD = L1


    logger.info(f"Saving DTI metric NIfTIs to {out_dir}...")
    save_nifti(f"{out_dir}/fa.nii.gz", FA.astype(np.float32), affine)
    save_nifti(f"{out_dir}/fa_inverse.nii.gz", FA_INVERSE.astype(np.float32), affine)
    save_nifti(f"{out_dir}/md.nii.gz", MD.astype(np.float32), affine)
    save_nifti(f"{out_dir}/ad.nii.gz", AD.astype(np.float32), affine)
    save_nifti(f"{out_dir}/rd.nii.gz", RD.astype(np.float32), affine)
    
    save_nifti(f"{out_dir}/rgb.nii.gz", np.array(255 * RGB, "uint8"), affine)
    save_nifti(f"{out_dir}/l1.nii.gz", L1.astype(np.float32), affine)
    save_nifti(f"{out_dir}/l2.nii.gz", L2.astype(np.float32), affine)
    save_nifti(f"{out_dir}/l3.nii.gz", L3.astype(np.float32), affine)

    if save_intermediate:
        save_nifti(f"{out_dir}/masked_dwi.nii.gz", maskdata.astype(np.float32), affine)

    return {"fa": FA, "fa_inverse": FA_INVERSE, 
            "md": MD, "rgb": RGB, "rd": RD, "ad": AD, 
            "l1": L1, "l2": L2, "l3": L3, "mask": mask}

def main(argv=None):
    p = argparse.ArgumentParser(description='Fit DTI and write FA/MD/RD/AD/RGB/L1/L2/L3 maps')
    p.add_argument('--dwi', required=True, help='DWI NIfTI (.nii or .nii.gz)')
    p.add_argument('--bvecs', required=True, help='bvecs file')
    p.add_argument('--bvals', required=True, help='bvals file')
    p.add_argument('--outdir', required=True, help='Output directory')
    p.add_argument('--mask', required=False, help='Optional brain mask NIfTI')
    p.add_argument('--nprocs', type=int, default=1, help='Number of processes (if supported)')
    p.add_argument('--save-intermediate', action='store_true',
                   help='Save masked_dwi.nii.gz alongside metric maps (omitted by default to save ~100MB)')
    p.add_argument('--verbose', action='store_true')
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        stream=sys.stdout)

    logging.info(f"Running DTI fit on {args.dwi} with bvecs {args.bvecs} and bvals {args.bvals}")

    compute_dti_metrics(args.dwi, args.bvals, args.bvecs, args.mask,
                        out_dir=args.outdir, save_intermediate=args.save_intermediate)


if __name__ == '__main__':
    main()