#!/usr/bin/env python3
"""
Inspect aparc+aseg label values to debug TraCULA issues.

This script checks:
1. What labels are present in the parcellation
2. Whether critical structures exist (cortex, thalamus, etc.)
3. Compares against expected FreeSurfer label ranges
"""

import argparse
import sys
from pathlib import Path

try:
    import nibabel as nib
    import numpy as np
except ImportError:
    print("ERROR: nibabel not found. Install with: pip install nibabel")
    sys.exit(1)


# Standard FreeSurfer label ranges
CRITICAL_LABELS = {
    "Left-Cerebral-Cortex": 3,
    "Right-Cerebral-Cortex": 42,
    "Left-Cerebral-White-Matter": 2,
    "Right-Cerebral-White-Matter": 41,
    "Left-Thalamus-Proper": 10,
    "Right-Thalamus-Proper": 49,
    "Left-Putamen": 12,
    "Right-Putamen": 51,
    "Left-Caudate": 11,
    "Right-Caudate": 50,
}

CORTICAL_LABELS = {
    "Left-Cerebral-Cortex": 3,
    "Right-Cerebral-Cortex": 42,
}

SUBCORTICAL_LABELS = {
    "Left-Thalamus": (10, 13),
    "Right-Thalamus": (49, 52),
}


def load_aparc(aparc_path):
    """Load aparc+aseg file and return data array."""
    try:
        img = nib.load(aparc_path)
        return img.get_fdata().astype(int)
    except Exception as e:
        print(f"ERROR: Failed to load {aparc_path}: {e}")
        sys.exit(1)


def analyze_labels(data):
    """Analyze label values in aparc+aseg."""
    unique_labels = np.unique(data)
    label_counts = {label: np.sum(data == label) for label in unique_labels}

    return unique_labels, label_counts


def check_critical_structures(unique_labels):
    """Check if critical structures are present."""
    missing = []
    present = []

    for name, label_val in CRITICAL_LABELS.items():
        if label_val in unique_labels:
            present.append(f"✓ {name} (label {label_val})")
        else:
            missing.append(f"✗ {name} (label {label_val})")

    return present, missing


def main():
    parser = argparse.ArgumentParser(
        description="Inspect aparc+aseg labels for TraCULA debugging"
    )
    parser.add_argument(
        "aparc_file",
        type=Path,
        help="Path to aparc+aseg.mgz or converted atlas"
    )
    parser.add_argument(
        "--show-all-labels",
        action="store_true",
        help="Show all labels found (not just critical ones)"
    )
    parser.add_argument(
        "--max-voxels",
        type=int,
        default=5,
        help="Show top N labels by voxel count (default: 5)"
    )

    args = parser.parse_args()

    if not args.aparc_file.exists():
        print(f"ERROR: File not found: {args.aparc_file}")
        sys.exit(1)

    print(f"\n📊 Analyzing: {args.aparc_file}")
    print("=" * 70)

    data = load_aparc(args.aparc_file)
    unique_labels, label_counts = analyze_labels(data)

    print(f"\n📈 Label Statistics:")
    print(f"   Total unique labels: {len(unique_labels)}")
    print(f"   Label range: {unique_labels.min()} - {unique_labels.max()}")
    print(f"   Total voxels: {data.size:,}")

    print(f"\n🔍 Critical Structures:")
    present, missing = check_critical_structures(unique_labels)

    if present:
        print("   Present:")
        for item in present:
            count = label_counts.get(CRITICAL_LABELS[[k for k in CRITICAL_LABELS if CRITICAL_LABELS[k] in [int(item.split()[-1])]][0]], 0)
            print(f"      {item} ({count:,} voxels)")

    if missing:
        print("\n   ⚠️  Missing:")
        for item in missing:
            print(f"      {item}")

    if args.show_all_labels:
        print(f"\n📋 All Labels (top {args.max_voxels} by voxel count):")
        sorted_labels = sorted(label_counts.items(), key=lambda x: x[1], reverse=True)
        for i, (label, count) in enumerate(sorted_labels[:args.max_voxels], 1):
            pct = 100 * count / data.size
            print(f"   {i:2d}. Label {label:3d}: {count:>10,} voxels ({pct:5.2f}%)")

    # Key warnings
    print(f"\n⚠️  Diagnostic Summary:")
    if len(unique_labels) < 10:
        print("   ⚠️  Very few labels found - parcellation may be incomplete")
    if 3 not in unique_labels and 42 not in unique_labels:
        print("   ⚠️  No cortical labels found - label conversion may have failed")
    if 10 not in unique_labels and 49 not in unique_labels:
        print("   ⚠️  No thalamus labels found - TraCULA may fail")

    print()


if __name__ == "__main__":
    main()
