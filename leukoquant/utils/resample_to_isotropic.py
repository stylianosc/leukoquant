import nibabel as nib
from nibabel.processing import resample_to_output
import argparse
import os

def convert_to_isotropic(input_path, output_path, target_resolution=(1.0, 1.0, 1.0)):
    """
    Convert a NIfTI image to isotropic resolution using nibabel's resample_to_output.
    Args:
        input_path: Path to the input NIfTI file.
        output_path: Path where the isotropic NIfTI file will be saved.
        target_resolution: Tuple of (x, y, z) voxel sizes in mm for the output image.
    """
    # Load the input image
    img = nib.load(input_path)

    # Resample to isotropic resolution
    isotropic_img = resample_to_output(img, voxel_sizes=target_resolution)

    # Save the isotropic image
    nib.save(isotropic_img, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a NIfTI image to isotropic resolution.")
    parser.add_argument('--input', required=True, help='Path to the input NIfTI file')
    parser.add_argument('--output', required=False, help='Path to save the isotropic NIfTI file')
    parser.add_argument('--resolution', nargs=3, type=float, default=(1.0, 1.0, 1.0),
                        help='Target voxel size in mm (x y z), default is 1.0mm isotropic')

    args = parser.parse_args()

    if not args.output:
        input_dir = os.path.dirname(args.input)
        output_path = os.path.join(input_dir, "isotropic_" + os.path.basename(args.input))
    else:
        output_path = args.output

    convert_to_isotropic(args.input, output_path, tuple(args.resolution))