import warnings
import logging
# Suppress verbose/repetitive nibabel warnings (e.g. qfac/pixdim[0] warnings)
logging.getLogger('nibabel').setLevel(logging.ERROR)
warnings.filterwarnings("ignore", category=UserWarning, module="nibabel")

import nibabel as nib
import numpy as np
import argparse
import sys
import os
import subprocess
from pathlib import Path
from dipy.io.streamline import load_tractogram, save_tractogram
from dipy.tracking import utils
from dipy.io.stateful_tractogram import StatefulTractogram, Space
from dipy.tracking.streamline import transform_streamlines
from scipy.ndimage import map_coordinates, center_of_mass

logger = logging.getLogger(__name__)


def _make_sft_safe(streamlines, reference, space, source_label=""):
    """Construct a StatefulTractogram, clamping any out-of-bounds streamline
    points to the volume boundary rather than discarding entire streamlines.

    Warping operations (affine or non-linear) can nudge a small number of
    streamline points just outside the target volume FOV - common when atlas
    tracts designed for a large template (MGH35, ~180 mm z) are warped to a
    dataset acquired with limited coverage (e.g. EPAD at 120 mm z). Removing
    entire streamlines because one endpoint clips 1 mm past the edge would bias
    density maps and connectivity metrics for those tracts. Clamping instead
    pulls those points to the nearest in-bounds coordinate, preserving the full
    tract for all intra-FOV measurements.
    """
    # DIPY >= 1.9 removed bbox_valid_check; the constructor no longer validates
    # bbox on creation, so pre-clamp coordinates are safe to pass in.
    sft = StatefulTractogram(streamlines, reference, space=space)
    sft.to_vox()

    dims = np.array(sft.dimensions, dtype=float)
    # StatefulTractogram.is_bbox_in_vox_valid() (called by save_tractogram)
    # does self.to_corner() before comparing against dims - i.e. it checks
    # center_coord + 0.5 <= dims, not center_coord <= dims, because we are in
    # the default Origin.NIFTI ("center") convention here. A point clamped to
    # `dims` (or even `dims - 1e-6`) in center-origin space is therefore still
    # `> dims` once shifted to corner-origin for the check, and the file we
    # just "fixed" fails to save with the exact same bbox error. The true
    # ceiling in center-origin space is dims - 0.5; np.nextafter in float32
    # then trims a hair further so the value also survives the float32
    # downcast when written to the .trk file.
    upper = np.nextafter((dims - 0.5).astype(np.float32), np.float32(0)).astype(np.float64)

    n_clamped = 0
    clamped = []
    for sl in sft.streamlines:
        clamped_sl = np.clip(sl, 0.0, upper)
        if not np.array_equal(clamped_sl, sl):
            n_clamped += 1
        clamped.append(clamped_sl)

    if n_clamped > 0:
        label = f" in '{source_label}'" if source_label else ""
        logger.warning(
            "Clamped %d/%d streamlines%s to FOV boundary (points from warp fell outside volume).",
            n_clamped, len(clamped), label,
        )

    return StatefulTractogram(clamped, reference, space=Space.VOX)


def get_transformed_tractogram(input_path, matrix, out_reference_path):
    """Apply an affine matrix to all streamlines and sync the header to a new reference."""
    sft = load_tractogram(filename=input_path, reference="same")
    transformed_points = transform_streamlines(sft.streamlines, matrix)
    return _make_sft_safe(transformed_points, out_reference_path, Space.RASMM, input_path)


def get_transformed_tractogram_nonlinear(
    input_path,
    deform_field_path,
    out_reference_path,
    init_affine_path=None,
    n_iter=20,
):
    """Warp streamlines from floating (MGH35) to reference (FA) space.

    Uses the inverse of a NiftyReg forward deformation field produced by:
        reg_transform -ref FA -cpp CPP -def deform_field.nii.gz

    That field is defined on the FA grid; each voxel stores the world-mm
    coordinate in MGH35 space that it maps to (forward: FA_world → MGH35_world).

    Inversion uses fixed-point iteration initialised from the affine-only
    estimate.  Update rule (subscript n = iteration):
        p_{n+1} = p_n + inv_A @ (q - def(p_n))
    where q is the MGH35 world-mm target, p is the current FA world-mm
    estimate, and A is the affine (FA_world → MGH35_world) used as the
    approximate Jacobian so the correction is expressed in FA mm.

    NiftyReg affine convention (-ref FA -flo MGH35 -aff AFF):
        AFF maps FA_world → MGH35_world  (reference → floating)
        inv(AFF) maps MGH35_world → FA_world  (used for init and Jacobian)
    """
    def_img = nib.load(deform_field_path)
    def_data = def_img.get_fdata()           # (X, Y, Z, 3), values = MGH35 world mm
    # Squeeze in case of singleton time dimension: (X, Y, Z, 1, 3) → (X, Y, Z, 3)
    def_data = np.squeeze(def_data)
    def_inv_affine = np.linalg.inv(def_img.affine)

    # Affine-only inverse for initialisation and as approximate Jacobian.
    if init_affine_path is not None:
        aff = np.loadtxt(init_affine_path)
        if aff.shape == (3, 4):
            aff = np.vstack([aff, [0, 0, 0, 1]])
        inv_aff = np.linalg.inv(aff)         # MGH35_world → FA_world
    else:
        # Fallback: identity (only valid if FA and MGH35 share world coords).
        inv_aff = np.eye(4)

    sft = load_tractogram(input_path, reference="same")
    sft.to_rasmm()

    new_streamlines = []
    for sl in sft.streamlines:
        N = len(sl)
        sl_hom = np.hstack([sl, np.ones((N, 1))]).T     # (4, N)

        # Affine-only initialisation: inv_aff @ q → FA world mm
        p = (inv_aff @ sl_hom)[:3].T                    # (N, 3) in FA world mm

        for _ in range(n_iter):
            # FA world mm → FA voxel coordinates
            p_hom = np.hstack([p, np.ones((N, 1))]).T
            vox = (def_inv_affine @ p_hom)[:3]           # (3, N)

            # Sample forward deformation field → predicted MGH35 world mm
            predicted = np.zeros((N, 3))
            for dim in range(3):
                # map_coordinates expects a tuple of coordinate arrays, one per dimension
                predicted[:, dim] = map_coordinates(
                    def_data[..., dim], tuple(vox), order=1, mode='nearest'
                )

            # Residual in MGH35 mm; map back to FA mm via inverse affine
            residual_mgh = sl - predicted                # (N, 3) in MGH35 mm
            residual_hom = np.hstack([residual_mgh, np.zeros((N, 1))]).T
            residual_fa = (inv_aff @ residual_hom)[:3].T  # (N, 3) in FA mm

            p = p + residual_fa

        new_streamlines.append(p)

    return _make_sft_safe(new_streamlines, out_reference_path, Space.RASMM, input_path)


def _densify_streamlines(streamlines, voxel_sizes, factor=0.5):
    """Interpolate points along streamlines so that no step is larger than factor * minimum voxel size."""
    step_size = factor * np.min(voxel_sizes)
    densified = []
    for sl in streamlines:
        diffs = np.diff(sl, axis=0)
        dists = np.linalg.norm(diffs, axis=1)
        if np.any(dists > step_size):
            new_sl = [sl[0]]
            for i, d in enumerate(dists):
                p0 = sl[i]
                p1 = sl[i+1]
                if d > step_size:
                    n_steps = int(np.ceil(d / step_size))
                    for step in range(1, n_steps):
                        new_sl.append(p0 + (step / n_steps) * (p1 - p0))
                new_sl.append(p1)
            densified.append(np.array(new_sl, dtype=np.float32))
        else:
            densified.append(sl)
    return densified


def get_density_map_nifti(input_path, normalise=True):
    """Return a NIfTI density map for the tractogram in its native space."""
    sft = load_tractogram(filename=input_path, reference="same")
    n = len(sft.streamlines)
    densified = _densify_streamlines(sft.streamlines, sft.voxel_sizes)
    dm = utils.density_map(densified, sft.affine, sft.dimensions).astype("float32")
    if normalise and n > 0:
        dm /= n
    return nib.Nifti1Image(dm, sft.affine)


def run_reg_resample(floating_path, reference_path, affine_file, output_path):
    """Run reg_resample to bring floating into reference space using an affine."""
    cmd = (
        f"reg_resample -ref {reference_path} -flo {floating_path} "
        f"-trans {affine_file} -inter 1 -res {output_path}"
    )
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)


def verify_density_maps(
    trk_path,
    ref_nifti_path,
    ref_image_path,
    com_threshold_mm=15.0,
    dice_threshold=0.3,
    output_verify_vol_path=None,
):
    """Verify a transformed trk by comparing its density map to a reference NIfTI.

    Both the trk and ref_nifti_path must already be in the same space
    (subject/T1/FA space).  A density map is generated from the trk
    streamlines in that space and compared to the reference NIfTI using
    Dice coefficient and CoM distance.  Prints a VERIFY OK / WARNING line.

    Used to cross-check:
    - tracula trks (affine-warped to T1) vs the reg_resample density map
      produced in the same command - both represent the same tract in T1 space
    - atlas trks (non-linearly warped to FA) vs the CPP-registered atlas
      density map - catches a wrong warp direction early
    """
    def _com_world(mask, affine):
        if not np.any(mask):
            return None
        vox = np.array(center_of_mass(mask))
        return (affine @ np.append(vox, 1.0))[:3]

    # Density map from the transformed trk in reference image space
    sft = load_tractogram(trk_path, reference=ref_image_path)
    ref_img = nib.load(ref_image_path)
    densified = _densify_streamlines(sft.streamlines, sft.voxel_sizes)
    dm = utils.density_map(densified, sft.affine, sft.dimensions)
    trk_mask = dm > 0

    if output_verify_vol_path:
        Path(output_verify_vol_path).parent.mkdir(parents=True, exist_ok=True)
        # Create a float32 NIfTI density map image normalised by total streamlines
        n = len(sft.streamlines)
        dm_float = dm.astype("float32")
        if n > 0:
            dm_float /= n
        dm_img = nib.Nifti1Image(dm_float, ref_img.affine)
        nib.save(dm_img, output_verify_vol_path)

    ref_nii = nib.load(ref_nifti_path)
    ref_mask = ref_nii.get_fdata() > 0

    # Dice on binary masks
    inter = np.logical_and(trk_mask, ref_mask).sum()
    denom = trk_mask.sum() + ref_mask.sum()
    dice = float(2 * inter / denom) if denom > 0 else 0.0

    # CoM distance in world mm
    trk_com = _com_world(trk_mask, ref_img.affine)
    ref_com = _com_world(ref_mask, ref_nii.affine)

    if trk_com is None or ref_com is None:
        print("VERIFY [WARNING]: one of the masks is empty - cannot compute metrics.")
        return False

    dist = float(np.linalg.norm(trk_com - ref_com))
    ok = dice >= dice_threshold and dist <= com_threshold_mm
    status = "OK" if ok else "WARNING"
    # print(
    #     f"VERIFY [{status}]: Dice={dice:.3f} (thr≥{dice_threshold}), "
    #     f"CoM dist={dist:.1f} mm (thr≤{com_threshold_mm} mm)"
    # )
    return ok


def main():
    parser = argparse.ArgumentParser(
        description="Transform TRK files and generate density maps."
    )
    parser.add_argument("--input-trk", required=True)
    parser.add_argument("--output-trk", required=False)
    parser.add_argument("--output-vol", required=False)
    parser.add_argument("--in-reference", required=False)
    parser.add_argument("--out-reference", required=False)
    parser.add_argument("--affine", required=False, help="4x4 affine matrix (.txt)")
    parser.add_argument(
        "--convert-affine",
        action="store_true",
        help="Re-header the trk to match --in-reference without moving streamlines",
    )
    parser.add_argument(
        "--warp",
        required=False,
        help="NiftyReg forward deformation field (.nii.gz) for non-linear trk warping "
             "(produced by reg_transform -ref FA -cpp CPP -def field.nii.gz)",
    )
    parser.add_argument(
        "--verify-nifti",
        required=False,
        help="NIfTI mask to compare the output trk CoM against (verification only)",
    )
    parser.add_argument(
        "--verify-output-vol",
        required=False,
        help="Optional path to output the generated density map from the verification TRK",
    )
    args = parser.parse_args()

    if not args.output_trk and not args.output_vol:
        print("Error: at least one of --output-trk or --output-vol must be specified.")
        sys.exit(1)

    # Load affine matrix once if provided
    mat = None
    if args.affine:
        try:
            mat = np.loadtxt(args.affine)
            if mat.shape == (3, 4):
                mat = np.vstack([mat, [0, 0, 0, 1]])
            elif mat.shape != (4, 4):
                print(f"Error: affine must be 3x4 or 4x4, got {mat.shape}")
                sys.exit(1)
        except Exception as exc:
            print(f"Error reading affine: {exc}")
            sys.exit(1)

    try:
        # 1. Affine trk transform
        if args.output_trk and mat is not None and args.out_reference and not args.warp:
            new_sft = get_transformed_tractogram(
                args.input_trk, np.linalg.inv(mat), args.out_reference
            )
            save_tractogram(new_sft, args.output_trk)

        # 2. Non-linear trk transform (fixed-point inversion of forward deformation field)
        if args.output_trk and args.warp and args.out_reference:
            new_sft = get_transformed_tractogram_nonlinear(
                input_path=args.input_trk,
                deform_field_path=args.warp,
                out_reference_path=args.out_reference,
                init_affine_path=args.affine,   # used for initialisation; optional
            )
            save_tractogram(new_sft, args.output_trk)

        # 3. Density map generation
        if args.output_vol:
            dm_img = get_density_map_nifti(args.input_trk)
            Path(args.output_vol).parent.mkdir(parents=True, exist_ok=True)
            nib.save(dm_img, args.output_vol)

            if args.out_reference and mat is not None:
                tmp = str(Path(args.output_vol).with_suffix("").with_suffix("")) + "_dwi_space.nii.gz"
                os.rename(args.output_vol, tmp)
                run_reg_resample(
                    floating_path=tmp,
                    reference_path=args.out_reference,
                    affine_file=args.affine,
                    output_path=args.output_vol,
                )
                os.remove(tmp)

        # 4. Header-only re-affine (no coordinate transform).
        #
        # This re-interprets the atlas trk's existing voxel-index coordinates
        # against --in-reference's affine/dims, on the assumption that they
        # already line up with that grid. When they don't quite (edge
        # rounding in the vendor atlas file), a point can land just outside
        # [0, dims) - invisible here since this used raw nibabel I/O with no
        # bbox check, but fatal on the next load_tractogram() downstream
        # (metrics_calc.sh's density-map step), which does validate. Route
        # through _make_sft_safe so that check happens now, with the same
        # clamp-to-boundary behaviour as the affine/nonlinear warp paths
        # above, rather than crashing one step later with a less useful
        # traceback.
        if args.convert_affine and args.in_reference and args.output_trk:
            raw = nib.streamlines.load(args.input_trk)
            new_sft = _make_sft_safe(
                raw.streamlines, args.in_reference, Space.VOX, args.input_trk
            )
            save_tractogram(new_sft, args.output_trk)

        # 5. Optional verification: generate density map from output trk and compare
        # to the reference NIfTI using Dice + CoM distance.
        if args.verify_nifti and args.output_trk and Path(args.output_trk).exists():
            if args.out_reference and Path(args.out_reference).exists():
                verify_density_maps(
                    args.output_trk,
                    args.verify_nifti,
                    args.out_reference,
                    output_verify_vol_path=args.verify_output_vol,
                )
            else:
                print("VERIFY [SKIP]: --out-reference required for density-map verification.")

    except Exception as exc:
        print(f"Processing failed: {exc}")
        sys.exit(1)


if __name__ == "__main__":
    main()
