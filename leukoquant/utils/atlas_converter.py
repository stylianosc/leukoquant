"""
Brain atlas converter with spatial clustering.

Converts parcellations between atlases using CSV mapping files.
Supports 1:1 direct mappings and 1:n spatial clustering for anatomical accuracy.
"""

import time
from pathlib import Path
from typing import Optional, Union, Dict, List, Tuple, Set
import nibabel as nib
import numpy as np
import pandas as pd


def convert_atlas_to_freesurfer(input_parcellation: Union[str, Path],
                               mapping_file: Union[str, Path],
                               output_path: Optional[Union[str, Path]] = None,
                               validate: bool = True,
                               verbose: bool = False) -> str:
    """
    Convert brain parcellation to FreeSurfer format.

    Supports 1:1 direct mappings and 1:n spatial clustering using nearest
    neighbor references for anatomically accurate distributions.

    Parameters
    ----------
    input_parcellation : str or Path
        Input parcellation NIfTI file
    mapping_file : str or Path
        CSV mapping file (Label, Source, Target, Nearest)
    output_path : str or Path, optional
        Output file path. Auto-generated if None
    validate : bool, default=True
        Run validation checks
    verbose : bool, default=False
        Show detailed progress

    Returns
    -------
    str
        Path to converted file
    """
    # Start timing
    start_time = time.time()

    # Load mapping
    if verbose:
        print(f"Loading mapping: {mapping_file}")
    mapping_df = pd.read_csv(mapping_file)

    # Load input parcellation
    if verbose:
        print(f"Loading input: {input_parcellation}")
    input_img = nib.load(str(input_parcellation))
    input_data = input_img.get_fdata().astype(np.int32)

    # Find unique labels in input (excluding background 0)
    input_labels = set(np.unique(input_data))
    input_labels.discard(0)  # Remove background
    
    if verbose:
        print(f"Input contains {len(input_labels)} unique labels: {sorted(input_labels)}")

    # Analyze mapping structure and classify mapping types
    if verbose:
        print("Analyzing mapping structure...")
    mapping_analysis = _analyze_mapping_structure(mapping_df, verbose)

    # Find unmapped labels (in input but not in mapping CSV)
    unmapped_labels = input_labels - mapping_analysis['all_source_labels']
    if unmapped_labels:
        print(f"⚠️  Warning: {len(unmapped_labels)} labels in input not found in mapping: {sorted(unmapped_labels)}")
        print("These labels will remain as background (0) in the output")

    # Initialize temporary array for safe processing
    temp_data = np.zeros_like(input_data, dtype=np.int32)

    if verbose:
        print("Applying mappings...")

    # Process 1:1 mappings (direct replacement)
    for source_label in mapping_analysis['one_to_one'].keys():
        if source_label in input_labels:
            target_label = mapping_analysis['one_to_one'][source_label]
            mask = (input_data == source_label)
            temp_data[mask] = target_label

    # Process 1:n mappings (spatial clustering)
    for source_label, mapping_info in mapping_analysis['one_to_many'].items():
        if source_label in input_labels:
            _apply_spatial_clustering(
                input_data, temp_data, source_label, mapping_info, verbose
            )

    # Safe copy to final output
    output_data = temp_data.copy()

    # Generate output path if not provided
    if output_path is None:
        input_path = Path(input_parcellation).absolute()
        input_filename = str(input_path).split('/')[-1].split('.', 1)[0]
        input_suffix = str(input_path).split('/')[-1].split('.', 1)[-1]
        
        output_path = input_path.parent / f"{input_filename}_freesurfer.{input_suffix}"

    output_path = Path(output_path).absolute()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save output
    output_img = nib.Nifti1Image(output_data, input_img.affine, input_img.header)
    nib.save(output_img, str(output_path))

    if verbose:
        print(f"Conversion completed: {output_path}")

    # Enhanced validation
    if validate:
        _validate_enhanced_conversion(input_data, output_data, verbose, unmapped_labels)

    # Calculate and print elapsed time with conversion details
    elapsed_time = time.time() - start_time

    # Detect atlas types from mapping file columns
    columns = list(mapping_df.columns)
    source_atlas = columns[1]  # Source column name
    target_atlas = columns[2]  # Target column name

    # Format conversion type
    conversion_type = f"{source_atlas} → {target_atlas}"

    print(f"⏱️  {conversion_type} parcellation conversion completed in {elapsed_time:.2f}s")
    print(f"Saved to: {output_path}")

    return str(output_path)


def _analyze_mapping_structure(mapping_df: pd.DataFrame, verbose: bool = False) -> Dict:
    """Analyze mapping structure and classify 1:1 vs 1:n mappings."""
    # Get column names
    columns = list(mapping_df.columns)
    source_col = columns[1]  # Source column
    target_col = columns[2]  # Target column
    nearest_col = columns[3]  # Nearest column

    # Group by source label to collect all targets for each source
    source_to_targets = {}  # source → list of all targets
    source_to_nearest = {}  # source → list of nearest refs

    for _, row in mapping_df.iterrows():
        # Parse source labels
        source_str = str(row[source_col])
        if not source_str or source_str.lower() == 'nan' or source_str.strip() == '':
            continue

        source_labels = _parse_labels(source_str)

        # Parse target labels
        target_str = str(row[target_col])
        if not target_str or target_str.lower() == 'nan' or target_str.strip() == '':
            continue

        target_labels = _parse_labels(target_str)

        # Parse nearest references
        nearest_str = str(row[nearest_col])
        nearest_refs = []
        if nearest_str and nearest_str.lower() != 'nan' and nearest_str.strip():
            nearest_refs = _parse_labels(nearest_str)

        # Collect targets for each source label
        for source_label in source_labels:
            if source_label not in source_to_targets:
                source_to_targets[source_label] = []
            if source_label not in source_to_nearest:
                source_to_nearest[source_label] = []

            source_to_targets[source_label].extend(target_labels)
            source_to_nearest[source_label].extend(nearest_refs)

    # Classify mapping types and collect all source labels
    all_source_labels = set(source_to_targets.keys())
    one_to_one = {}  # source → single target
    one_to_many = {}  # source → {targets, nearest_refs}

    for source_label in source_to_targets:
        unique_targets = list(set(source_to_targets[source_label]))  # Remove duplicates

        if len(unique_targets) == 1:
            # 1:1 mapping (even if appears in multiple CSV rows)
            one_to_one[source_label] = unique_targets[0]
        else:
            # 1:n mapping (multiple unique targets)
            one_to_many[source_label] = {
                'targets': unique_targets,
                'nearest': source_to_nearest[source_label]
            }

    if verbose:
        print(f"Mapping analysis: {len(one_to_one)} direct, {len(one_to_many)} spatial")
        print(f"All source labels in mapping: {sorted(all_source_labels)}")

        print(f"One-to-one mappings: {one_to_one}")
        print(f"One-to-many mappings: {one_to_many}")

    return {
        'all_source_labels': all_source_labels,
        'one_to_one': one_to_one,
        'one_to_many': one_to_many
    }


def _parse_labels(label_string: str) -> List[int]:
    """Parse label string supporting space and comma separation."""
    if not label_string or label_string.strip() == '':
        return []

    labels = []
    # Try space separation first
    if ' ' in label_string:
        parts = label_string.split()
        for part in parts:
            part = part.strip()
            if part.isdigit():
                labels.append(int(part))
    # Try comma separation
    elif ',' in label_string:
        parts = label_string.split(',')
        for part in parts:
            part = part.strip()
            if part.isdigit():
                labels.append(int(part))
    # Single number
    try:
        label = int(float(label_string.strip()))  # Convert to float first
    except ValueError:
        label = None
    if label is not None and str(label).isdigit():
        labels = [label]

    return labels


def _apply_spatial_clustering(input_data: np.ndarray, temp_data: np.ndarray,
                            source_label: int, mapping_info: Dict,
                            verbose: bool = False) -> None:
    """Apply spatial clustering for 1:n mappings using nearest references."""
    targets = mapping_info['targets']
    nearest_refs = mapping_info['nearest']

    if len(nearest_refs) != 2:
        if verbose:
            print(f"WARNING: Expected 2 nearest refs for {source_label}, got {len(nearest_refs)}")
        return

    # Calculate center of mass for nearest reference regions
    nearest_centers = {}
    for nearest_label in nearest_refs:
        mask = (input_data == nearest_label)
        if np.sum(mask) > 0:
            coords = np.where(mask)
            center = (
                np.mean(coords[0]),
                np.mean(coords[1]),
                np.mean(coords[2])
            )
            nearest_centers[nearest_label] = center

    if len(nearest_centers) != 2:
        if verbose:
            print(f"WARNING: Could not find both nearest regions for {source_label}")
        return

    # Get the two centers
    center1 = nearest_centers[nearest_refs[0]]
    center2 = nearest_centers[nearest_refs[1]]

    # Create vector between centers
    vector = np.array(center2) - np.array(center1)

    # Generate equidistant cluster centers
    cluster_centers = {}
    n_targets = len(targets)
    for i, target_label in enumerate(targets):
        t = (i + 1) / (n_targets + 1)  # Equidistant positions
        center = tuple(np.array(center1) + t * vector)
        cluster_centers[target_label] = center

    if verbose:
        print(f"SPATIAL: {source_label} → {targets} (nearest: {nearest_refs})")

    # Assign each source voxel to closest cluster
    source_mask = (input_data == source_label)
    source_coords = np.where(source_mask)

    for x, y, z in zip(source_coords[0], source_coords[1], source_coords[2]):
        # Find closest cluster center
        min_dist = float('inf')
        closest_target = targets[0]  # fallback

        voxel_pos = np.array([x, y, z])
        for target_label, center in cluster_centers.items():
            dist = np.linalg.norm(voxel_pos - np.array(center))
            if dist < min_dist:
                min_dist = dist
                closest_target = target_label

        # Assign to temporary array
        temp_data[x, y, z] = closest_target


def _validate_enhanced_conversion(input_data: np.ndarray, output_data: np.ndarray,
                                verbose: bool = False, unmapped_labels: Set[int] = None):
    """Enhanced validation of conversion results."""
    # Count unmapped voxels (labels not in CSV source column)
    unmapped_voxels = 0
    if unmapped_labels:
        for label in unmapped_labels:
            unmapped_voxels += np.sum(input_data == label)

    total_voxels = np.sum(input_data != 0)  # Non-background voxels
    unmapped_percentage = (unmapped_voxels / total_voxels * 100) if total_voxels > 0 else 0

    if verbose:
        print(f"Enhanced validation: {unmapped_percentage:.1f}% unmapped voxels")

        if unmapped_labels:
            print(f"Truly unmapped labels: {sorted(unmapped_labels)}")

        output_labels = sorted(np.unique(output_data))
        if output_labels:
            print(f"Output labels: {output_labels[:15]}{'...' if len(output_labels) > 15 else ''}")

    if unmapped_percentage > 10:
        print(f"⚠️  Warning: {unmapped_percentage:.1f}% of voxels unmapped")
        if unmapped_labels:
            print(f"   Unmapped labels: {sorted(unmapped_labels)}")

    success = unmapped_percentage < 10
    if verbose:
        status = "✅ SUCCESS" if success else "❌ FAILURE"
        print(f"   Validation result: {status}")

    return success


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="Convert brain parcellation to FreeSurfer format")
    parser.add_argument("--input-parcellation", "-i", required=True,
                       help="Input parcellation NIfTI file")
    parser.add_argument("--mapping-file", "-m", required=True,
                       help="CSV mapping file")
    parser.add_argument("--output-path", "-o",
                       help="Output file path (optional, auto-generated if not provided)")
    parser.add_argument("--no-validate", action="store_true",
                       help="Skip validation")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Verbose output")

    args = parser.parse_args()

    try:
        result = convert_atlas_to_freesurfer(
            input_parcellation=args.input_parcellation,
            mapping_file=args.mapping_file,
            output_path=args.output_path,
            validate=not args.no_validate,
            verbose=args.verbose
        )
        print(f"Conversion completed: {result}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
