import pandas as pd
import warnings
#!/usr/bin/env python3
"""
Tractography QC utility script for leukoquant.

This script performs QC on tractography outputs for a subject.
It is called from the Snakemake workflow and can also be run standalone.

Arguments:
    --subject            Subject ID
    --tractography-path  Path to tractography directory or file
    --output             Output CSV file for QC report
"""

import argparse
import os
from pathlib import Path
import glob
import numpy as np
from dipy.io.streamline import load_tractogram
from dipy.tracking import utils


def find_tract_files(subject: str, tractography_path: str, skip_subject_dir: bool = False):
    """
    Discover tract files and their tract names for a subject, supporting colon-separated base:glob paths.

    When skip_subject_dir=True the subject ID is NOT appended to the base path,
    meaning the base already points directly at the subject's tractography folder.

    Returns: (tract_files, tract_names, exists)
    """
    tract_files = []
    tract_names = []
    print(f"[tract_qc] find_tract_files: subject={subject!r}, tractography_path={tractography_path!r}, skip_subject_dir={skip_subject_dir}")

    if ":" in tractography_path:
        # Support base:glob[:space], where optional space (e.g., dwi/t1) should
        # not be treated as part of the filesystem glob.
        parts = tractography_path.split(":")
        base = parts[0]
        glob_pattern = parts[1] if len(parts) > 1 else ""
        base_path = Path(base)
        # When the base already includes the subject directory, don't append it again.
        search_root = base_path if skip_subject_dir else base_path / subject
        search_pattern = str(search_root / glob_pattern) if glob_pattern else str(search_root)
        print(f"[tract_qc] base={base!r}, glob_pattern={glob_pattern!r}")
        print(f"[tract_qc] base_path exists: {base_path.exists()} ({base_path})")
        print(f"[tract_qc] search_root exists: {search_root.exists()} ({search_root})")
        print(f"[tract_qc] search_pattern: {search_pattern!r}")
        tract_files = [str(f) for f in glob.glob(search_pattern)]
        print(f"[tract_qc] tract_files found ({len(tract_files)}): {tract_files}")
        tract_names = [Path(f).parent.name for f in tract_files]
        exists = len(tract_files) > 0
    else:
        tractography_path_obj = Path(tractography_path)
        print(f"[tract_qc] no colon - checking path directly: {tractography_path_obj}")
        print(f"[tract_qc] is_dir: {tractography_path_obj.is_dir()}, is_file: {tractography_path_obj.is_file()}")
        if tractography_path_obj.is_dir():
            for f in tractography_path_obj.rglob("*"):
                if f.is_file():
                    tract_files.append(str(f))
                    tract_names.append(f.parent.name)
            print(f"[tract_qc] tract_files found in dir ({len(tract_files)}): {tract_files[:5]}{'...' if len(tract_files) > 5 else ''}")
            exists = len(tract_files) > 0
        elif tractography_path_obj.is_file():
            tract_files.append(str(tractography_path_obj))
            tract_names.append(tractography_path_obj.parent.name)
            exists = True
        else:
            print(f"[tract_qc] path does not exist: {tractography_path_obj}")
            exists = False
    return tract_files, tract_names, exists

def qc_metrics(sft):
    """
    Compute QC metrics for a StatefulTractogram. Returns a dict of metric_name: value.
    Extend this function to add more metrics as needed.
    """
    metrics = {}
    streamlines = sft.streamlines
    voxel_dims = sft.voxel_sizes
    lengths_mm = np.array(list(utils.length(streamlines)), dtype=float) if streamlines else np.array([], dtype=float)
    # Metric 1: Streamline length STD
    streamlines_length_std = float(np.std(lengths_mm)) if lengths_mm.size > 0 else None
    # Metric 2: Volume to median streamline length ratio
    if lengths_mm.size > 0:
        median_length_mm = float(np.median(lengths_mm))
        # For volume, we calculate the density map, filter with 0.2 * max value,
        # and sum the remaining voxels. This is a proxy for the "volume" of the tract.
        dm = utils.density_map(streamlines, sft.affine, sft.dimensions)
        threshold = 0.2 * np.max(dm)
        volume_mm3 = np.sum(dm >= threshold) * np.prod(voxel_dims)

        ratio = volume_mm3 / median_length_mm if median_length_mm > 0 else None
    else:
        ratio = None
        volume_mm3 = None
        median_length_mm = None

    metrics = {
        "streamline_length_std": streamlines_length_std,
        "volume_median_streamline_length_ratio": ratio,
        "volume_mm3": volume_mm3,
        "median_streamline_length_mm": median_length_mm,
    }
    # Decision logic for pass/fail and reason using criteria
    criteria = {
        "streamline_length_std": 10**(-3),
        #"volume_median_streamline_length_ratio": 2
    }
    failed = []
    for key, threshold in criteria.items():
        value = metrics.get(key)
        if value is None or value < threshold:
            failed.append(f"{key}={value} (< {threshold})")
    if not failed:
        metrics["qc_pass"] = True
        metrics["reason"] = ""
    else:
        metrics["qc_pass"] = False
        metrics["reason"] = "; ".join(failed)
    return metrics

def qc_tract_file(tract_file: str) -> dict:
    """
    Perform QC on a single tract file. Returns a dict with QC results.
    For now, checks if file exists and is non-empty. Extend as needed.
    """
    qc_result = {
        "file": tract_file,
        "exists": False,
        "size": 0,
        "qc_pass": False,
        "reason": "",
        "streamline_length_std": None,
        "volume_median_streamline_length_ratio": None
    }

    try:
        if not os.path.isfile(tract_file):
            qc_result["reason"] = "File does not exist"
            qc_result["qc_pass"] = False
            return qc_result
        qc_result["exists"] = True
        qc_result["size"] = os.path.getsize(tract_file)
        if qc_result["size"] == 0:
            qc_result["reason"] = "File is empty"
            qc_result["qc_pass"] = False
            return qc_result
        if not tract_file.lower().endswith('.trk'):
            warnings.warn(f"File {tract_file} is not a .trk file. Only .trk files are supported for QC.")
            qc_result["reason"] = "Not a .trk file. Only .trk files are supported for QC."
            qc_result["qc_pass"] = False
            return qc_result
        try:
            sft = load_tractogram(tract_file, reference="same")
            metrics = qc_metrics(sft)
            qc_result.update(metrics)
            qc_result["qc_pass"] = metrics["qc_pass"]
            qc_result["reason"] = metrics["reason"]
        except Exception as e:
            qc_result["reason"] = f"Error in .trk QC: {e}"
            qc_result["qc_pass"] = False
    except Exception as e:
        qc_result["reason"] = f"Error: {e}"
        qc_result["qc_pass"] = False
    return qc_result

from typing import Optional

def run_tract_qc(subject: str, tractography_path: str, output_path: Optional[str] = None, skip_subject_dir: bool = False):
    tract_files, tract_names, exists = find_tract_files(subject, tractography_path, skip_subject_dir=skip_subject_dir)
    print(f"[tract_qc] run_tract_qc: {len(tract_files)} tract file(s) found for subject {subject!r}")
    for f in tract_files:
        print(f"[tract_qc]   -> {f}")
    qc_results = [qc_tract_file(f) for f in tract_files]
    report = {
        "subject": subject,
        "tractography_path": tractography_path,
        "exists": exists,
        "tract_files": tract_files,
        "tract_names": tract_names,
        "qc_results": qc_results,
        "notes": "QC passed" if exists else "Tractography file not found"
    }
    # Handle output_path as directory or file
    if not output_path:
        output_dir = os.getcwd()
        output_csv = os.path.join(output_dir, "qc_report.csv")
    elif os.path.isdir(output_path):
        output_dir = output_path.rstrip("/")
        output_csv = os.path.join(output_dir, "qc_report.csv")
    else:
        print(f"Output {output_path}, Is Dir: {os.path.isdir(output_path)}, \
              Is Dir Rstrip: {os.path.isdir(output_path.rstrip('/'))}")
        output_dir = os.path.dirname(output_path)
        output_csv = output_path

    os.makedirs(output_dir, exist_ok=True)
    print(f"Output Dir: {output_dir}, Output Path: {output_path}, Output CSV: {output_csv}")
    # Write all fields from all QC result dicts using pandas
    df = pd.DataFrame(report["qc_results"])
    df.to_csv(output_csv, index=False)
    print(f"QC results written to {output_csv}")

def main():
    parser = argparse.ArgumentParser(description="Tractography QC utility")
    parser.add_argument('--subject', required=True, help='Subject ID')
    parser.add_argument('--tractography-path', required=True, help='Path to tractography directory or file')
    parser.add_argument('--output', required=False, help='Output directory or file for QC report (CSV). If not set, uses current directory.')
    parser.add_argument('--skip-subject-dir', action='store_true', default=False,
                        help='Do not append the subject ID to the tractography base path (use when base already contains the subject folder).')
    args = parser.parse_args()
    run_tract_qc(args.subject, args.tractography_path, args.output, skip_subject_dir=args.skip_subject_dir)

if __name__ == "__main__":
    main()
