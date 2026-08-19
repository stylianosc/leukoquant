"""
Patch the FreeSurfer orientLAS script to handle PSL/PSR oriented DWI images.

orientLAS (a tcsh script called by trac-preproc) has two bugs for PSL/PSR images:

1. mri_info --dim, --cdc, --rdc, --sdc, --cras, --cres, --rres, --sres emit
   "INFO: ..." lines to stdout for non-standard orientations.  tcsh backtick
   substitution captures those tokens alongside the numeric output, corrupting
   $dims/$cdc/etc. and producing a garbled mri_convert call such as:
       mri_convert -oni is -onj INFO: -onk This ...

2. The orientation string is extracted with awk '{print $3}', which is fragile
   to varying field widths in the mri_info "Orientation" output line.

Usage:
    python patch_orientlas.py <path_to_orientLAS>

The file is modified in place.
"""

import sys


def patch(path: str) -> None:
    with open(path) as f:
        content = f.read()

    original = content

    # Fix 1: pipe every mri_info --xxx call through grep to strip INFO: lines
    # that contaminate tcsh backtick captures for PSL/PSR images.
    for flag in ("--dim", "--cdc", "--rdc", "--sdc", "--cras", "--cres", "--rres", "--sres"):
        content = content.replace(
            f"`mri_info {flag} $in`",
            f'`mri_info {flag} $in |& grep -v "^INFO:"`',
        )

    # Fix 2: use $NF (last field) instead of $3 to extract the orientation and
    # determinant strings so the extraction is robust to whitespace variations.
    for label in ("Orientation", "determinant"):
        content = content.replace(
            f"grep {label} | awk " + "'{print $3}'",
            f"grep {label} | awk " + "'{print $NF}'",
        )

    n_info = content.count('grep -v "^INFO:"')
    n_nf = content.count("print $NF")

    if content == original:
        print("WARNING: patch_orientlas.py made no changes - patterns not found", flush=True)
    else:
        print(f"patch_orientlas: {n_info}/8 INFO filters, {n_nf}/2 $NF fixes applied", flush=True)
        with open(path, "w") as f:
            f.write(content)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <path_to_orientLAS>", file=sys.stderr)
        sys.exit(1)
    patch(sys.argv[1])
