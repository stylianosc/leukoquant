"""
Validate FA-to-MGH35-template registration quality before trusting it for the
tract atlas fallback (metrics_calc.sh, "Preparing Tract Atlas" section).

transform_trk.py's --verify-nifti (used per atlas tract, downstream of this
check) only checks that two derivations of THIS SAME warp agree with each
other -- it cannot detect a genuinely bad reg_aladin/reg_f3d registration,
since both derivations would be consistently wrong together and still
"agree". This checks registration accuracy directly instead: correlate the
subject's own FA against the MGH35 template resampled into subject space
(produced by reg_f3d). A failed/degenerate registration (e.g. reg_f3d
converging to a bad local optimum) shows up here as a near-zero or negative
NCC (normalized cross-correlation), well before any tract gets warped from it.

Usage:
    python check_atlas_registration.py <subject_fa> <resampled_template_fa> [--threshold 0.4]

Prints "PASS <ncc>" or "FAIL <ncc>" to stdout and exits 0 on PASS, 1 on FAIL,
so the caller can act on either the exit code or the stdout text.
"""

import argparse
import sys
from typing import Optional

import nibabel as nib
import numpy as np

MIN_OVERLAP_VOXELS = 100


def compute_registration_ncc(subject_fa_path: str, resampled_template_path: str) -> Optional[float]:
    """Normalized cross-correlation between the subject's FA map and the
    template FA resampled into subject space. Returns None if the two images
    have too little non-background overlap to compute a meaningful NCC."""
    a = np.asarray(nib.load(subject_fa_path).dataobj).astype(float)
    b = np.asarray(nib.load(resampled_template_path).dataobj).astype(float)
    mask = (a > 0) | (b > 0)
    if mask.sum() < MIN_OVERLAP_VOXELS:
        return None
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def main() -> None:
    parser = argparse.ArgumentParser(description="Check FA→MGH35 template registration quality via NCC")
    parser.add_argument("subject_fa", help="Subject's own FA map")
    parser.add_argument("resampled_template", help="MGH35 template FA resampled into subject space")
    parser.add_argument("--threshold", type=float, default=0.4, help="Minimum NCC to pass (default: 0.4)")
    args = parser.parse_args()

    ncc = compute_registration_ncc(args.subject_fa, args.resampled_template)
    if ncc is None:
        print("FAIL 0.0")
        sys.exit(1)

    status = "PASS" if ncc >= args.threshold else "FAIL"
    print(f"{status} {ncc}")
    sys.exit(0 if status == "PASS" else 1)


if __name__ == "__main__":
    main()
