"""Shared image utility functions for the leukoquant pipeline."""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

logger = logging.getLogger(__name__)


def reorient_to_ras(input_path: str, output_path: str, clip_negatives: bool = True) -> None:
    """Reorient a NIfTI image to RAS+ and (by default) clip negative voxels to 0.

    Uses nibabel.as_closest_canonical(), which flips/swaps data axes to reach
    the nearest RAS-aligned orientation without resampling or interpolation.

    Negative voxel values are never legitimate signal for a magnitude
    modality like T1/FLAIR (they're interpolation overshoot or bias-field
    correction artifacts) and are clipped to 0 - large-magnitude negatives
    (e.g. -32603) are known to crash mideface's mri_defacer with heap
    corruption.

    The file is only re-saved if reorientation or clipping actually changed
    something; otherwise it's copied verbatim (no re-save overhead, no
    floating-point precision loss).

    Args:
        input_path:  Path to the source NIfTI file (.nii or .nii.gz).
        output_path: Destination path for the reoriented, clipped image.
        clip_negatives: Clip negative voxels to 0 (default). Pass False to only
            reorient, leaving every voxel value untouched (used by the metrics
            stage, which must land on the same grid as the BaMoS/GIF outputs
            without altering the T1 intensities it was previously given).
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    img = nib.load(str(input_path))
    axcodes = nib.aff2axcodes(img.affine)
    needs_reorient = axcodes != ("R", "A", "S")

    img_ras = nib.as_closest_canonical(img) if needs_reorient else img
    # The voxel data is only read when negatives have to be looked for; a pure reorientation never needs it.
    data = np.asarray(img_ras.dataobj) if clip_negatives else None
    has_negatives = clip_negatives and bool(np.any(data < 0))

    is_same = Path(input_path).resolve() == Path(output_path).resolve()

    if not needs_reorient and not has_negatives:
        if not is_same:
            logger.info("Image already RAS+ with no negative voxels, copying as-is: %s", input_path)
            shutil.copy2(str(input_path), str(output_path))
        else:
            logger.info("Image already RAS+ with no negative voxels (in-place): %s", input_path)
        return

    if has_negatives:
        data = np.where(data < 0, 0, data)
        img_out = nib.Nifti1Image(data, img_ras.affine, img_ras.header)
    else:
        img_out = img_ras

    if is_same:
        tmp_out = Path(output_path).parent / f".tmp_{os.getpid()}_{Path(output_path).name}"
        nib.save(img_out, str(tmp_out))
        os.replace(str(tmp_out), str(output_path))
    else:
        nib.save(img_out, str(output_path))

    logger.info(
        "Processed %s -> %s (reoriented=%s from %s, clipped_negatives=%s)",
        input_path,
        output_path,
        needs_reorient,
        "".join(axcodes),
        has_negatives,
    )


def select_first_frame(input_path: str, output_path: str) -> None:
    """Collapse a 4D NIfTI with more than one volume down to its first volume.

    Some raw EPAD T1/FLAIR acquisitions are exported as a 4D image with a
    duplicate second frame (confirmed identical to the first, not a genuine
    second acquisition) instead of a normal 3D volume. Tools downstream
    (recon-all, GIF, BaMoS) all require a single 3D volume and either reject
    multi-frame input outright (recon-all: "input(s) cannot have multiple
    frames!") or silently misbehave.

    Only the first volume is kept. Uses `Nifti1Image.slicer` rather than
    manually rebuilding a header around a scaled float array, so any
    dtype/scl_slope/scl_inter encoding on the source file is handled the
    same way nibabel itself would for any other 4D-to-3D slice.

    If the image is already 3D (or 4D with exactly one volume), it's copied
    verbatim -- no re-save overhead, no requantization, and no assumption
    that a legitimately-single-volume 4D image needs "fixing".

    Args:
        input_path:  Path to the source NIfTI file (.nii or .nii.gz).
        output_path: Destination path for the single-volume image.
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    img = nib.load(str(input_path))

    is_same = Path(input_path).resolve() == Path(output_path).resolve()

    if img.ndim < 4 or img.shape[3] <= 1:
        if not is_same:
            logger.info("Image is already single-volume, copying as-is: %s", input_path)
            shutil.copy2(str(input_path), str(output_path))
        else:
            logger.info("Image is already single-volume (in-place): %s", input_path)
        return

    n_frames = img.shape[3]
    img_out = img.slicer[..., 0]
    if is_same:
        tmp_out = Path(output_path).parent / f".tmp_{os.getpid()}_{Path(output_path).name}"
        nib.save(img_out, str(tmp_out))
        os.replace(str(tmp_out), str(output_path))
    else:
        nib.save(img_out, str(output_path))

    logger.info(
        "Selected frame 0 of %d from %s -> %s",
        n_frames,
        input_path,
        output_path,
    )


def same_voxel_grid(path_a: str, path_b: str, atol_mm: float = 1e-2) -> bool:
    """True when two NIfTI images share one voxel grid: same 3-D shape and the same affine (within atol_mm).

    Header-only. Used to detect an image that is already on the target T1 grid before it would be resampled from DWI
    space with the DWI->T1 transform (which would be wrong for it).
    """
    a, b = nib.load(path_a), nib.load(path_b)
    return tuple(a.shape[:3]) == tuple(b.shape[:3]) and bool(np.allclose(a.affine, b.affine, rtol=0, atol=atol_mm))


def main(argv=None) -> None:
    p = argparse.ArgumentParser(
        description="Reorient a NIfTI image to RAS+, or collapse a multi-frame image to its first volume."
    )
    p.add_argument("--input",   required=True, help="Input NIfTI path (.nii or .nii.gz)")
    p.add_argument("--output",  help="Output NIfTI path (not used by same_grid)")
    p.add_argument("--ref",     help="Reference NIfTI for --op same_grid")
    p.add_argument(
        "--op",
        choices=["reorient", "canonical", "first_frame", "same_grid"],
        default="reorient",
        help="Operation to perform (default: reorient)",
    )
    p.add_argument("--verbose", action="store_true", help="Enable INFO-level logging")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s %(levelname)s: %(message)s",
        stream=sys.stdout,
    )

    if args.op == "same_grid":
        # Exit status 0 = same voxel grid, 1 = different (usable directly in a shell `if`).
        if not args.ref:
            p.error("--op same_grid needs --ref")
        sys.exit(0 if same_voxel_grid(args.input, args.ref) else 1)
    if not args.output:
        p.error("--output is required for this operation")
    if args.op == "first_frame":
        select_first_frame(args.input, args.output)
    elif args.op == "canonical":
        reorient_to_ras(args.input, args.output, clip_negatives=False)
    else:
        reorient_to_ras(args.input, args.output)


if __name__ == "__main__":
    main()
