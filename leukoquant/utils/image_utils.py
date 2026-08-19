"""Shared image utility functions for the leukoquant pipeline."""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

logger = logging.getLogger(__name__)


def reorient_to_ras(input_path: str, output_path: str) -> None:
    """Reorient a NIfTI image to RAS+ and clip negative voxels to 0.

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
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    img = nib.load(str(input_path))
    axcodes = nib.aff2axcodes(img.affine)
    needs_reorient = axcodes != ("R", "A", "S")

    img_ras = nib.as_closest_canonical(img) if needs_reorient else img
    data = np.asarray(img_ras.dataobj)
    has_negatives = bool(np.any(data < 0))

    if not needs_reorient and not has_negatives:
        logger.info("Image already RAS+ with no negative voxels, copying as-is: %s", input_path)
        shutil.copy2(str(input_path), str(output_path))
        return

    if has_negatives:
        data = np.where(data < 0, 0, data)
        img_out = nib.Nifti1Image(data, img_ras.affine, img_ras.header)
    else:
        img_out = img_ras

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

    if img.ndim < 4 or img.shape[3] <= 1:
        logger.info("Image is already single-volume, copying as-is: %s", input_path)
        shutil.copy2(str(input_path), str(output_path))
        return

    n_frames = img.shape[3]
    img_out = img.slicer[..., 0]
    nib.save(img_out, str(output_path))
    logger.info(
        "Selected frame 0 of %d from %s -> %s",
        n_frames,
        input_path,
        output_path,
    )


def main(argv=None) -> None:
    p = argparse.ArgumentParser(
        description="Reorient a NIfTI image to RAS+, or collapse a multi-frame image to its first volume."
    )
    p.add_argument("--input",   required=True, help="Input NIfTI path (.nii or .nii.gz)")
    p.add_argument("--output",  required=True, help="Output NIfTI path")
    p.add_argument(
        "--op",
        choices=["reorient", "first_frame"],
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

    if args.op == "first_frame":
        select_first_frame(args.input, args.output)
    else:
        reorient_to_ras(args.input, args.output)


if __name__ == "__main__":
    main()
