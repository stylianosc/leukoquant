#!/usr/bin/env python
"""Utilities for TRACULA workflows.

Provides a small helper to render a dmrirc template for a subject.
The script is intended to be called inside the Singularity container as:

python /leukoquant/leukoquant/utils/tracula_utils.py \
    --subject <subj> --subjects-dir <subjects_dir> --outdir /output/<subj>

It writes `/output/<subj>/dmrirc.config`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Dict, List, Optional, Tuple, Union
import shutil
import glob

import nibabel as nib
import numpy as np

from leukoquant.utils.dwi_utils import normalize_bvec_bval_to_tracula


def is_shelled_acquisition(
    bvals_path: Union[str, List[str]],
    b_range: int = 50,
    max_shells: int = 10,
    min_vols_per_shell: int = 10,
) -> bool:
    """Return True if every non-b0 shell has enough volumes for eddy's PEAS predictor.

    Accepts a single bvals file path or a list of paths; when a list is given the
    b-values are combined before grouping so multi-scan sessions are handled correctly.

    Uses the same b_range tolerance as FSL eddy (default 50 s/mm²). The max_shells
    parameter is retained for logging but no longer affects the return value - a
    dataset with many well-populated shells is still suitable for PEAS. Only
    under-populated shells (min_vols < min_vols_per_shell) warrant --dont_peas.

    Args:
        bvals_path: Path (or list of paths) to bvals file(s) (FSL single-row format).
        b_range: Maximum spread within a shell - matches eddy's --b_range default.
        max_shells: Informational only; no longer gates the return value.
        min_vols_per_shell: Minimum volumes required in every non-b0 shell.

    Returns:
        True if all non-b0 shells meet min_vols_per_shell; False otherwise.
    """
    # Normalize to list and merge all b-values.
    paths = [bvals_path] if isinstance(bvals_path, str) else list(bvals_path)
    all_bvals: list[float] = []
    for bp in paths:
        all_bvals.extend(np.loadtxt(bp).flatten().tolist())

    # Greedy grouping: assign each b-value (sorted) to the first existing shell
    # whose centre is within b_range, otherwise open a new shell.
    shell_centers: list[float] = []
    shell_counts: dict[float, int] = {}
    for bval in sorted(all_bvals):
        matched = next(
            (c for c in shell_centers if abs(bval - c) <= b_range), None
        )
        if matched is not None:
            shell_counts[matched] += 1
        else:
            shell_centers.append(bval)
            shell_counts[bval] = 1

    non_b0 = {c: n for c, n in shell_counts.items() if c > b_range}

    if not non_b0:
        return False  # No DWI volumes - not shelled

    n_shells = len(non_b0)
    min_vols = min(non_b0.values())
    # max_shells is logged but no longer gates; only min_vols matters for dont_peas.
    shelled = min_vols >= min_vols_per_shell

    print(f"Debug: b-value shells detected: {shell_counts}")
    print(f"Debug: non-b0 shells: {non_b0}")
    print(f"Debug: shelled: {shelled}")
    print(f"Debug: n_shells={n_shells}, min_vols_per_shell={min_vols}")
    print(f"Debug: Shell centers: {shell_centers}")

    return shelled


def prepare_freesurfer_dir_for_gif(
    subject_id: str,
    gif_parcellation: str,
    recon_all_dir: str,
    output_dir: str,
) -> str:
    """
    Prepares a FreeSurfer subject directory for TRACULA using a GIF parcellation.

    This function creates a new subject directory in the specified output directory.
    It then copies the necessary files from the original recon-all directory and
    the converted GIF parcellation into the new directory.

    Args:
        subject_id: The subject ID.
        gif_parcellation: Path to the converted GIF parcellation (aparc+aseg.mgz).
        recon_all_dir: Path to the original FreeSurfer recon-all directory.
        output_dir: Path to the output directory where the new subject directory will be created.

    Returns:
        Path to the new subject directory.
    """
    # Create the new subject directory
    new_subject_dir = Path(output_dir) / subject_id
    new_subject_dir.mkdir(parents=True, exist_ok=True)

    # Define the files and directories to copy from the original recon-all directory
    files_to_copy = [
        # Surface files
        "surf/lh.white", "surf/rh.white", "surf/lh.thickness", "surf/rh.thickness", "surf/rh.sphere.reg", "surf/lh.sphere.reg", 
        "surf/rh.pial", "surf/lh.pial",
        # Label files
        "label/lh.cortex.label", "label/rh.cortex.label",
        # Original T1 
        "mri/orig.mgz", 
        # Brain mask
        "mri/brain.mgz", "mri/brainmask.mgz", 
        # norm.mgz: Intensity-standardised volume (WM scaled to 110) created by 'mri_ca_normalize' using the Gaussian Classifier Array (GCA) model, nu.mgz, talairach.lta, and brainmask.mgz.
        "mri/norm.mgz", 
        # White matter segmentation
        "mri/wm.mgz",
        "mri/transforms/talairach.xfm", "mri/aseg.presurf.mgz",
        # Binary volume mask for cortical ribbon
        "mri/ribbon.mgz", 
    ]

    # Copy the files from the original recon-all directory
    for item in files_to_copy:
        source_path = Path(recon_all_dir) / subject_id / item
        destination_path = new_subject_dir / item
        destination_path.parent.mkdir(parents=True, exist_ok=True)

        for file in glob.glob(str(source_path)):
            if Path(file).is_dir():
                shutil.copytree(file, destination_path, dirs_exist_ok=True)
            else:
                shutil.copy(file, destination_path)

    # Copy the converted GIF parcellation
    shutil.copy(gif_parcellation, new_subject_dir / "mri" / "aparc+aseg.mgz")

    return str(new_subject_dir)


def pe_direction_converter(pe_dir: str) -> str:
    """Convert phase encoding direction from BIDS format to TRACULA format.
    If already in TRACULA format, returns as is. If unrecognized, return empty string.

    BIDS uses "i", "j", "k" with optional "-" suffix for reverse direction.
    TRACULA expects RL, LR, PA, AP, IS, SI.

    Examples:
        "i" -> "RL"
        "i-" -> "LR"
        "j" -> "PA"
        "j-" -> "AP"
        "k" -> "IS"
        "k-" -> "SI"
    """

    pe_dir = pe_dir.strip().lower()
    if pe_dir in ["i", "i-"]:
        return "RL" if pe_dir == "i" else "LR"
    elif pe_dir in ["j", "j-"]:
        return "PA" if pe_dir == "j" else "AP"
    elif pe_dir in ["k", "k-"]:
        return "IS" if pe_dir == "k" else "SI"
    elif pe_dir in ["rl", "lr", "pa", "ap", "is", "si"]:
        return pe_dir.upper()
    else:
        return ""  # Unrecognized format


def derive_epi_factor(total_readout_time: float, echo_spacing: float) -> int:
    """Derive the EPI echo train length from timing parameters.

    For a single-shot EPI acquisition the echo train length equals the full
    reconstruction matrix size in the phase-encode direction, which can be
    recovered as:

        epi_factor = round(TotalReadoutTime / EffectiveEchoSpacing + 1)

    The +1 accounts for the first k-space line carrying zero echo-spacing
    offset.

    Args:
        total_readout_time: Total EPI readout time in seconds.
        echo_spacing: Effective echo spacing in seconds.

    Returns:
        Derived EPI factor as an integer.

    Raises:
        ValueError: If echo_spacing is zero.
    """
    if echo_spacing == 0:
        raise ValueError("echo_spacing must be non-zero.")

    return round(total_readout_time / echo_spacing + 1)


_METADATA_KEYS: Dict[str, List[str]] = {
    "pe_dir": [
        "PhaseEncodingDirection",
        "PhaseEncoding",
        "InPlanePhaseEncodingDirectionDICOM",
        "pe_dir",
    ],
    "epi_factor": ["EchoTrainLength", "EPIFactor", "epi_factor"],
    "echo_spacing": [
        "EffectiveEchoSpacing",
        "EchoSpacing",
        "echo_spacing",
        "DerivedVendorReportedEchoSpacing",
    ],
    "total_readout_time": ["TotalReadoutTime", "total_readout_time"],
}


def _extract_dwi_metadata(dwi_path: str) -> Dict[str, str]:
    """Extract TRACULA-relevant metadata from the JSON sidecar of a single DWI file.

    Looks for a *.json file in the same directory as dwi_path (BIDS convention:
    same stem + .json).  Falls back to any *.json in that directory if the
    stem-matched file is absent.  Returns a dict with keys pe_dir, epi_factor,
    echo_spacing, and total_readout_time; empty strings for missing fields.

    epi_factor is derived from TotalReadoutTime / EffectiveEchoSpacing + 1 when
    no explicit EchoTrainLength / EPIFactor key is present.
    pe_dir is converted from BIDS (i/j/k) to TRACULA (RL/LR/PA/AP/IS/SI) format.
    """
    result = {k: "" for k in _METADATA_KEYS}

    dwi_p = Path(dwi_path)
    # Prefer the same-stem JSON (e.g. data_0.json or sub-01_dwi.json);
    # also accept metadata.json / metadata_N.json from dwi_utils;
    # fall back to any *.json in the same directory.
    stem = dwi_p.name.split(".")[0]
    stem_json = dwi_p.with_name(stem + ".json")
    metadata_json = dwi_p.with_name("metadata.json")
    # metadata_N.json produced by the workflow loop (e.g. metadata_0.json)
    scan_suffix = stem.rsplit("_", 1)[-1]  # e.g. "0" from "data_0"
    indexed_metadata_json = dwi_p.with_name(f"metadata_{scan_suffix}.json")

    if stem_json.exists():
        json_candidates = [stem_json]
    elif indexed_metadata_json.exists():
        json_candidates = [indexed_metadata_json]
    elif metadata_json.exists():
        json_candidates = [metadata_json]
    else:
        json_candidates = list(dwi_p.parent.glob("*.json"))

    if not json_candidates:
        print(f"Warning: no JSON sidecar found alongside {dwi_path}; "
              "pe_dir / epi_factor / echo_spacing will be left empty for this scan.")
        return result

    try:
        with json_candidates[0].open("r") as jf:
            meta = json.load(jf)
    except (OSError, json.JSONDecodeError) as exc:
        print(f"Warning: failed to read JSON sidecar {json_candidates[0]}: {exc}; "
              "pe_dir / epi_factor / echo_spacing will be left empty for this scan.")
        return result

    for var_name, keys in _METADATA_KEYS.items():
        for key in keys:
            if key in meta:
                result[var_name] = str(meta[key])
                break

    # Derive epi_factor from timing fields when no explicit key is present.
    if not result["epi_factor"]:
        if result["total_readout_time"] and result["echo_spacing"]:
            derived = derive_epi_factor(
                float(result["total_readout_time"]),
                float(result["echo_spacing"]),
            )
            result["epi_factor"] = str(derived)
            print(f"Info: epi_factor not found in JSON for {dwi_path}; "
                  f"derived from TotalReadoutTime / EffectiveEchoSpacing + 1 = {derived}")
        else:
            print(f"Warning: epi_factor not found in JSON for {dwi_path} and cannot be "
                  "derived (missing TotalReadoutTime or EffectiveEchoSpacing).")

    # Convert PE direction from BIDS to TRACULA format.
    if result["pe_dir"]:
        result["pe_dir"] = pe_direction_converter(result["pe_dir"])

    return result


def _filter_dwi_by_dimensions(
    dwi_list: List[str],
    bvecs_list: List[str],
    bvals_list: List[str],
) -> Tuple[List[str], List[str], List[str]]:
    """Return the largest subset of DWI scans that share identical spatial dimensions.

    trac-preproc calls mri_concat on all DWI inputs, which requires every volume
    to have the same first-3-axis shape.  When scans disagree, keep the biggest
    group; ties are broken by taking the group whose scan indices are lowest
    (i.e. the earliest scans in sorted order).

    Returns the filtered (dwi, bvecs, bvals) lists unchanged when all scans
    already agree on dimensions, or when only one scan is present.
    """
    if len(dwi_list) <= 1:
        return dwi_list, bvecs_list, bvals_list

    # Group scan indices by spatial shape (dim1, dim2, dim3).
    groups: Dict[Tuple[int, ...], List[int]] = {}
    for idx, dwi_path in enumerate(dwi_list):
        try:
            img = nib.load(dwi_path)
            key = tuple(int(d) for d in img.shape[:3])
        except Exception as exc:
            print(f"Warning: could not read dimensions of {dwi_path}: {exc}; "
                  "this scan will be assigned dimension key (0, 0, 0).")
            key = (0, 0, 0)
        groups.setdefault(key, []).append(idx)

    if len(groups) == 1:
        # All scans agree - nothing to filter.
        print(f"Info: all {len(dwi_list)} DWI scan(s) share dimensions "
              f"{next(iter(groups))} - no filtering needed.")
        return dwi_list, bvecs_list, bvals_list

    # Select the largest group; on tie, prefer the one whose first index is
    # smallest (i.e. the group that contains the earliest scan).
    best_key = max(groups, key=lambda k: (len(groups[k]), -groups[k][0]))
    best_indices = sorted(groups[best_key])

    dropped = sum(len(v) for k, v in groups.items() if k != best_key)
    print(
        f"Warning: DWI scans have mismatched spatial dimensions. "
        f"Keeping {len(best_indices)} scan(s) with dimensions {best_key}; "
        f"dropping {dropped} scan(s) with incompatible dimensions "
        f"({', '.join(str(k) for k in groups if k != best_key)}). "
        f"mri_concat requires all inputs to share the same first-3-axis shape."
    )

    dwi_out   = [dwi_list[i]   for i in best_indices]
    bvecs_out = [bvecs_list[i] for i in best_indices]
    bvals_out = [bvals_list[i] for i in best_indices]
    return dwi_out, bvecs_out, bvals_out


def write_dmrirc(
    subject: str,
    subjects_dir: str,
    outdir: str,
    template_path: Optional[str] = None,
    dwi: Union[str, List[str], None] = None,
    bvecs: Union[str, List[str], None] = None,
    bvals: Union[str, List[str], None] = None,
    usethalnuc: Optional[int] = None,
    doeddy: Optional[int] = 1,
    dob0: Optional[int] = 0,
    dmrirc_out: Optional[str] = None,
    verbose: bool = False,
) -> Tuple[Path, bool]:
    """Render a dmrirc template and write it to disk.

    Accepts one or more DWI scans (and their matching bvec/bval files).
    When multiple scans are provided, dcmlist / bveclist / bvallist are written
    as space-separated basenames inside a single dcmroot directory, which is the
    TRACULA multi-scan format documented in Example 3 of the dmrirc wiki page.

    doeddy is set to 2 (FSL eddy) only when every scan has a JSON sidecar that
    provides pe_dir, epi_factor, AND echo_spacing.  If any scan is missing any
    of these three fields, doeddy falls back to 1 (eddy_correct).

    After dimension filtering, is_shelled_acquisition is called on the filtered
    bval files. When any non-b0 shell has fewer than 10 volumes, doeddy is
    automatically lowered to 1 (eddy_correct) to avoid FSL eddy PEAS failure
    (GetGlobal2DWIIndexMapping).

    If `dmrirc_out` is provided it will be used as the full output path;
    otherwise a file named `dmrirc.config` is written inside `outdir`.

    Returns:
        Tuple of (dmrirc_path, is_shelled) where is_shelled indicates whether
        the filtered scans are suitable for PEAS (True) or need --dont_peas (False).
    """
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    if dmrirc_out:
        dmrirc_out_path = Path(dmrirc_out)
        dmrirc_out_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        dmrirc_out_path = outdir_p / "dmrirc.config"

    TEMPLATE_NAME = "dmrirc.template"
    tpl_path = Path(template_path) if template_path else Path(__file__).parent / TEMPLATE_NAME
    if not tpl_path.exists():
        raise FileNotFoundError(f"dmrirc template not found at: {tpl_path}")

    # Normalise inputs to lists so the rest of the function is uniform.
    dwi_list:   List[str] = ([dwi]   if isinstance(dwi,   str) else list(dwi   or []))
    bvecs_list: List[str] = ([bvecs] if isinstance(bvecs, str) else list(bvecs or []))
    bvals_list: List[str] = ([bvals] if isinstance(bvals, str) else list(bvals or []))

    # Filter to the largest subset of scans that share the same spatial
    # dimensions - mri_concat (called by trac-preproc) requires this.
    dwi_list, bvecs_list, bvals_list = _filter_dwi_by_dimensions(
        dwi_list, bvecs_list, bvals_list
    )

    # trac-preproc expects TRACULA-native bvecs/bvals (N rows x 3 cols / N
    # rows of 1 value each), but these files were prepared upstream by
    # dwi_utils.prepare_dwi_input() (shared with dti/noddi), which always
    # normalizes to FSL format (3 rows x N cols / 1 row x N cols) instead.
    # Left uncorrected, trac-preproc reads the FSL bvecs file's 3 axis-rows
    # as if they were 3 gradient vectors, regardless of the real direction
    # count, and fails with "Found N b-values but 3 gradient vectors" for
    # any subject with more than 3 real directions -- confirmed against
    # production logs (OASIS3/ADNI3/EPAD/DPUK) and reproduced with
    # synthetic data. Convert back to TRACULA format here, right before the
    # files are referenced in the dmrirc template, since this function is
    # the one place that owns TRACULA's format contract.
    for bvec_path, bval_path in zip(bvecs_list, bvals_list):
        if bvec_path and bval_path:
            normalize_bvec_bval_to_tracula(bvec_path, bval_path)

    # Check if --dont_peas is needed using only the dimension-filtered bval files.
    is_shelled = is_shelled_acquisition(bvals_list) if bvals_list else False

    # ── per-scan metadata extraction ─────────────────────────────────────────
    pe_dirs:       List[str] = []
    epi_factors:   List[str] = []
    echo_spacings: List[str] = []

    for dwi_path in dwi_list:
        meta = _extract_dwi_metadata(dwi_path)
        pe_dirs.append(meta["pe_dir"])
        epi_factors.append(meta["epi_factor"])
        echo_spacings.append(meta["echo_spacing"])
        print(f"Extracted metadata for {dwi_path}: {meta}")

    # ── doeddy decision ───────────────────────────────────────────────────────
    # Use FSL eddy (doeddy=2) only when every scan has all three required fields
    # AND all non-b0 shells have sufficient volumes for PEAS.
    if doeddy is None:
        doeddy = 1

    if doeddy == 2:
        if not is_shelled:
            # Under-populated shells: FSL eddy's PEAS GP predictor requires ≥10 volumes
            # per shell and crashes with GetGlobal2DWIIndexMapping when that is not met.
            # Fall back to eddy_correct (doeddy=1) which has no such constraint.
            print(
                "Warning: one or more non-b0 shells have fewer than 10 volumes; "
                "setting doeddy=1 (eddy_correct) to avoid FSL eddy PEAS failure."
            )
            doeddy = 1
        else:
            any_missing_meta = any(
                not pe_dirs[i] or not epi_factors[i] or not echo_spacings[i]
                for i in range(len(dwi_list))
            )
            if any_missing_meta:
                print(
                    "Warning: one or more DWI scans are missing JSON metadata "
                    "(pe_dir / epi_factor / echo_spacing); "
                    "setting doeddy=1 (eddy_correct) for all scans."
                )
                doeddy = 1

    if dob0 is None:
        dob0 = 0

    # ── template variable assembly ────────────────────────────────────────────
    # All prepared files share a common parent directory (dcmroot).
    # bvec/bval filenames are listed as space-separated basenames.
    dcmroot   = str(Path(dwi_list[0]).parent) if dwi_list else ""
    dcmlist   = " ".join(Path(d).name for d in dwi_list)
    bveclist  = " ".join(Path(b).name for b in bvecs_list if b)
    bvallist  = " ".join(Path(b).name for b in bvals_list if b)
    pe_dir_val      = " ".join(pe_dirs)
    epi_factor_val  = " ".join(epi_factors)
    echo_spacing_val = " ".join(echo_spacings)

    # subjlist must repeat the subject ID once per DWI scan - TRACULA validates
    # that len(subjlist) == len(dcmlist) and errors with "must specify as many
    # subject IDs as DWI scans" if they differ.
    subjlist_value = " ".join([subject] * len(dwi_list))

    kwargs = {
        "subject":          subject,
        "subjlist_value":   subjlist_value,
        "subjects_dir":     subjects_dir,
        "dcmroot_value":    dcmroot,
        "dcmlist_value":    dcmlist,
        "bveclist_value":   bveclist,
        "bvallist_value":   bvallist,
        "usethalnuc":       usethalnuc if usethalnuc is not None else "",
        "doeddy":           doeddy,
        "dob0":             dob0,
        "output_dir_value": str(outdir_p.absolute()),
        "pe_dir":           pe_dir_val,
        "epi_factor":       epi_factor_val,
        "echo_spacing":     echo_spacing_val,
    }

    try:
        from string import Template
        with tpl_path.open("r") as f:
            tpl = Template(f.read())
            content = tpl.substitute(**kwargs)

        print("Rendered dmrirc content:")
        print(content)
        print("Variables replaced in template:")
        for key in kwargs:
            if f"${key}" in tpl.template:
                print(f"  {key}: {kwargs[key]}")
    except KeyError as e:
        raise KeyError(f"Template placeholder missing: {e}") from e

    dmrirc_out_path.write_text(content)
    return dmrirc_out_path, is_shelled


def main(argv=None):
    p = argparse.ArgumentParser(description="Write dmrirc for TRACULA")
    p.add_argument("--subject", required=True)
    p.add_argument("--subjects-dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--template", required=False)
    # Accept one or more DWI / bvec / bval paths to support multi-scan sessions.
    p.add_argument("--dwi",   nargs="+", required=False)
    p.add_argument("--bvecs", nargs="+", required=False)
    p.add_argument("--bvals", nargs="+", required=False)
    p.add_argument("--usethalnuc", required=False, type=int)
    p.add_argument("--doeddy", required=False, type=int)
    p.add_argument("--dob0", required=False, type=int)
    p.add_argument("--dmrirc-out", required=False)
    p.add_argument(
        "--is-shelled-flag-file",
        required=False,
        help="If provided, writes '1' (acquisition is shelled) or '0' (not shelled) to this path.",
    )
    p.add_argument("--verbose", action="store_true")

    args = p.parse_args(argv)

    try:
        dmrirc, is_shelled = write_dmrirc(
            args.subject,
            args.subjects_dir,
            args.outdir,
            args.template,
            dwi=args.dwi,
            bvecs=args.bvecs,
            bvals=args.bvals,
            usethalnuc=args.usethalnuc,
            doeddy=args.doeddy,
            dob0=args.dob0,
            dmrirc_out=args.dmrirc_out,
            verbose=args.verbose,
        )
        print(f"Wrote dmrirc: {dmrirc}")

        # Write the is_shelled flag if requested.
        if args.is_shelled_flag_file:
            flag_path = Path(args.is_shelled_flag_file)
            flag_path.parent.mkdir(parents=True, exist_ok=True)
            # Write 1 = shelled (PEAS ok), 0 = not shelled (may need --dont_peas or other action).
            flag_path.write_text("1" if is_shelled else "0")
            print(f"is_shelled flag written to {flag_path}: {'1' if is_shelled else '0'}")
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == '__main__':
    main()
