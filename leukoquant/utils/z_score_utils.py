"""Utilities for the Z-score computation workflow.

Provides path-resolution helpers used by ``z_score_workflow.smk`` and a
design-matrix generator invoked via its CLI by ``z_score_calc.sh``.
"""

import glob
import sys
import logging
import numpy as np
import pandas as pd
import argparse

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Path helpers (used at Snakemake workflow parse time)
# ---------------------------------------------------------------------------

def resolve_single_match(pattern: str, label: str) -> str:
    """Return the single file matched by *pattern*, raising if 0 or >1 found."""
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise ValueError(f"No {label} file found for pattern: {pattern}")
    if len(matches) > 1:
        raise ValueError(
            f"Multiple {label} files found for pattern: {pattern}. Matches: {matches}"
        )
    return matches[0]


def get_t1_path(subject: str, config: dict) -> str:
    """Resolve and return the T1 image path for *subject*.

    Supports two config layouts:
    - ``t1_pattern``: a glob string with ``{subject}`` placeholder (standalone workflow).
    - ``t1_map``: a dict mapping subject ID → absolute path (process-all workflow,
      where each subject's T1 is already resolved to a concrete path).
    """
    if "t1_map" in config:
        path = config["t1_map"].get(subject)
        if not path:
            raise ValueError(f"Subject '{subject}' not found in config t1_map")
        return path
    pattern = config["t1_pattern"].format(subject=subject)
    return resolve_single_match(pattern, "T1")


def get_metric_path(subject: str, metric: str, config: dict) -> str:
    """Resolve and return the metric image path for *subject* and *metric*.

    Metric values support two formats:
    - Plain glob: ``/base/{subject}/fa.nii.gz`` - used by the standalone workflow,
      resolved via filesystem glob (file must exist at DAG-build time).
    - Colon-separated: ``/base/{subject}/dti/outputs:fa.nii.gz:space`` - used by
      process-all, where the path is fully deterministic so no glob is performed
      (the file will be produced by an upstream rule before z-score runs).
    """
    raw = config["metrics"][metric].format(subject=subject)
    parts = raw.split(":")
    if len(parts) >= 2:
        # Deterministic path from process-all: skip glob, just join base + filename.
        base = parts[0].rstrip("/")
        filename = parts[1].lstrip("/")
        return f"{base}/{filename}"
    return resolve_single_match(raw, f"metric '{metric}'")


def get_dwi_path(subject: str, config: dict) -> str:
    """Resolve and return the DWI image path for *subject*.

    Supports two config layouts, mirroring get_t1_path():
    - ``dwi_map``: a dict mapping subject ID -> absolute path (process-all workflow).
    - ``dwi_pattern``: a glob string with ``{subject}`` placeholder (standalone workflow).
    """
    if "dwi_map" in config:
        path = config["dwi_map"].get(subject)
        if not path:
            raise ValueError(f"Subject '{subject}' not found in config dwi_map")
        return path
    if "dwi_pattern" not in config:
        raise ValueError("dwi_pattern not found in config")
    pattern = config["dwi_pattern"].format(subject=subject)
    return resolve_single_match(pattern, "DWI")


def get_bval_path(subject: str, config: dict) -> str:
    """Resolve and return the bval file path for *subject*.

    Supports two config layouts, mirroring get_t1_path()/get_dwi_path():
    - ``bval_map``: a dict mapping subject ID -> absolute path (process-all workflow).
    - ``bval_pattern``: a glob string with ``{subject}`` placeholder (standalone workflow).
    """
    if "bval_map" in config:
        path = config["bval_map"].get(subject)
        if not path:
            raise ValueError(f"Subject '{subject}' not found in config bval_map")
        return path
    if "bval_pattern" not in config:
        raise ValueError("bval_pattern not found in config")
    pattern = config["bval_pattern"].format(subject=subject)
    return resolve_single_match(pattern, "bval")


def translate_path(path: str, bind_map: dict) -> str:
    """Translate an absolute host path to its container-mounted equivalent.

    Uses longest-prefix match so nested directories resolve correctly.
    Returns the original path unchanged if no bind entry covers it.
    """
    for host_dir, container_dir in sorted(bind_map.items(), key=lambda x: -len(x[0])):
        if path == host_dir:
            return container_dir
        if path.startswith(host_dir + "/"):
            return container_dir + path[len(host_dir):]
    return path


# ---------------------------------------------------------------------------
# Design-matrix builder (heavy deps imported lazily)
# ---------------------------------------------------------------------------

def _process_continuous_term(vals, t_val, power=1):
    """Process a continuous covariate with optional polynomial expansion.
    
    Applies centering, then raises to the given power.
    
    Args:
        vals: Series of healthy values (numeric)
        t_val: Target value (numeric)
        power: Exponent (1 for linear, 2 for quadratic, etc.)
    
    Returns:
        Tuple of (processed_healthy_values, processed_target_value)
    """
    mean_val = vals.mean()
    centered = vals - mean_val
    t_centered = t_val - mean_val
    
    if power == 1:
        return centered, t_centered
    else:
        # Center first, then raise to power: (x - mean)^power
        powered = centered ** power
        t_powered = t_centered ** power
        return powered, t_powered


def generate_design_matrix(
    demo_path: str,
    target_id: str,
    healthy_ids_path: str,
    covariates_str: str,
    poly_str: str,
    design_mat_path: str,
    target_vals_path: str,
    verbose: bool = False,
) -> None:
    """Build an FSL GLM design matrix from demographics and write to disk.

    Parameters
    ----------
    demo_path:
        Path to the demographics CSV.  First column must be the subject ID.
    target_id:
        Subject ID of the target (patient) subject.
    healthy_ids_path:
        Path to a plain-text file listing healthy subject IDs (one per line).
    covariates_str:
        Comma-separated covariate column names, e.g. ``"Age,Sex"``.
        Pass an empty string for no covariates.
    poly_str:
        Comma-separated polynomial terms in format ``covariate^power``,
        e.g. ``"Age^2,Age^3,BMI^2"``. For backward compatibility, plain
        covariate names (e.g. ``"Age"``) default to power=2. Pass empty string for none.
    design_mat_path:
        Output path for the FSL-format design matrix (space-separated, no header).
    target_vals_path:
        Output path for a space-separated file of target regressor values.
    verbose:
        If True, log detailed debug information about the design matrix construction.
    """


    covariates = [c.strip().lower() for c in covariates_str.split(",") if c.strip()] if covariates_str else []
    
    # Parse polynomial terms: format is "covariate_power" (e.g., "age_2", "age_3")
    # This comes pre-formatted from the processor
    poly_terms = []
    if poly_str:
        for p in poly_str.split(","):
            p = p.strip()
            if not p:
                continue
            # Format: covariate_power (e.g., "Age_2") - already formatted by processor
            parts = p.rsplit("_", 1)  # Split from right to handle names with underscores
            if len(parts) != 2:
                raise ValueError(f"Invalid polynomial term format '{p}'. Expected 'covariate_power' format (e.g., 'Age_2')")
            cov_name = parts[0].strip().lower()
            try:
                power = int(parts[1].strip())
                if power < 1:
                    raise ValueError(f"Polynomial power must be >= 1, got {power} for term '{p}'")
            except ValueError as e:
                if "invalid literal" in str(e):
                    raise ValueError(f"Invalid polynomial power in term '{p}': must be an integer")
                raise
            poly_terms.append((cov_name, power))

    df = pd.read_csv(demo_path)
    df.columns = df.columns.str.strip().str.lower()

    # Prefer an explicit "subject" column (leukoquant's harmonized
    # clinical_data_nn.csv schema puts "session" first, "subject" second) --
    # only fall back to the first column for demographics files that don't
    # follow that schema.
    id_col = "subject" if "subject" in df.columns else df.columns[0]
    df[id_col] = df[id_col].astype(str)

    with open(healthy_ids_path, "r") as fh:
        healthy_ids = [line.strip() for line in fh if line.strip()]

    # Match healthy IDs with exact match first, then substring match (case-insensitive)
    matched_healthy_ids = []
    for hid in healthy_ids:
        # Try exact match
        exact_matches = df[df[id_col] == hid]
        if not exact_matches.empty:
            matched_healthy_ids.append(hid)
        else:
            # Try substring match
            substring_matches = df[df[id_col].astype(str).str.contains(hid, case=False, na=False)]
            if not substring_matches.empty:
                # Use the first matched ID
                matched_healthy_ids.append(substring_matches[id_col].iloc[0])
            else:
                raise ValueError(f"Healthy subject '{hid}' not found in demographics (tried exact and substring match)")

    healthy_df = df[df[id_col].isin(matched_healthy_ids)].copy()
    healthy_df = healthy_df.set_index(id_col).reindex(matched_healthy_ids).reset_index()

    if len(healthy_df) != len(healthy_ids):
        missing = set(healthy_ids) - set(healthy_df[id_col].tolist())
        raise ValueError(f"Missing demographics for healthy subjects: {missing}")

    # Try exact match first, then substring match (case-insensitive)
    target_row = df[df[id_col] == target_id]
    if target_row.empty:
        target_row = df[df[id_col].astype(str).str.contains(target_id, case=False, na=False)]
    if target_row.empty:
        available_ids = df[id_col].tolist()[:10]  # Show first 10 for debugging
        raise ValueError(
            f"Target subject '{target_id}' not found in {demo_path} (tried exact match and substring match). "
            f"Available IDs (first 10): {available_ids}"
        )

    if verbose:
        logger.info(f"\nTarget subject '{target_id}' demographics:")
        logger.info(f"  {target_row.to_dict('records')[0]}")
        logger.info(f"\nHealthy subjects count: {len(healthy_df)}")
        #logger.info(f"Healthy subject IDs: {healthy_df[id_col].tolist()}")

    # Build design matrix for HEALTHY subjects only (target is NOT included).
    # The GLM fits a model to healthy data, then we use target_regressors
    # to predict what the target SHOULD look like, and compare to actual.
    design_matrix = pd.DataFrame()
    target_regressors = []  # Will hold target's covariate values (processed identically)
    column_headers = []    # Track column names for commenting
    
    # Track continuous variables for polynomial expansion
    continuous_vars = {}  # Store (vals, t_val, is_categorical) for each covariate

    if verbose:
        logger.info(f"\nBuilding design matrix with covariates: {covariates}")
        logger.info(f"Polynomial terms: {poly_terms}")

    for cov in covariates:
        if cov not in df.columns:
            raise ValueError(f"Covariate '{cov}' not found in demographics columns: {df.columns.tolist()}")
        
        col_healthy = healthy_df[cov]
        target_val_raw = target_row[cov].values[0]
        
        # Check for missing values in healthy cohort
        if col_healthy.isna().any():
            missing_ids = healthy_df[col_healthy.isna()][id_col].tolist()
            raise ValueError(f"Healthy subjects {missing_ids} have missing values for covariate '{cov}'")
        
        # Check for missing target value
        if pd.isna(target_val_raw):
            raise ValueError(f"Target subject '{target_id}' has missing value for covariate '{cov}'")

        is_string = col_healthy.dtype == object or not pd.api.types.is_numeric_dtype(col_healthy)
        is_categorical = False  # Will be set to True if variable is categorical

        if verbose:
            logger.info(f"\n  Covariate '{cov}':")
            logger.info(f"    Raw target value: {target_val_raw}")

        if is_string:
            # String/object categorical: sort unique values across both cohort and
            # target, then map to 0, 1, 2, ...
            all_unique = sorted(
                set(str(v) for v in col_healthy.tolist()) | {str(target_val_raw)}
            )
            mapping = {str(val): float(i) for i, val in enumerate(all_unique)}
            vals  = col_healthy.astype(str).map(mapping).astype(float)
            t_val = float(mapping[str(target_val_raw)])
            is_categorical = True
            if verbose:
                logger.info(f"    Categorical mapping: {mapping}")
        else:
            vals  = col_healthy.astype(float)
            t_val = float(target_val_raw)
            # Remap numeric binary columns that are not already 0/1
            # (e.g. Sex encoded as 1/2)
            numeric_unique = sorted(set(vals.tolist()) | {t_val})
            if len(numeric_unique) == 2 and set(numeric_unique) != {0.0, 1.0}:
                bin_map = {v: float(i) for i, v in enumerate(numeric_unique)}
                vals  = vals.map(bin_map)
                t_val = bin_map[t_val]
                is_categorical = True  # Binary variables are treated as categorical
                if verbose:
                    logger.info(f"    Binary remapping: {bin_map}")

        if verbose:
            logger.info(f"    Mapped target value: {t_val}")
            logger.info(f"    Is categorical: {is_categorical}")

        # Store for potential polynomial expansion later
        continuous_vars[cov] = (vals, t_val, is_categorical)

        # Only demean continuous variables, NOT categorical ones (like Sex)
        if is_categorical:
            # For categorical variables, use raw coded values (no demeaning)
            col_vals = vals
            t_col_val = t_val
            col_name = cov
            if verbose:
                logger.info(f"    Categorical - using raw coded values (no demeaning)")
                #logger.info(f"    Healthy values: {col_vals.tolist()}")
                logger.info(f"    Target value: {t_col_val}")
        else:
            # For continuous variables, demean using the helper function
            col_vals, t_col_val = _process_continuous_term(vals, t_val, power=1)
            col_name = f"{cov}_demeaned"
            if verbose:
                logger.info(f"    Continuous - demeaning")
                logger.info(f"    Mean (healthy): {vals.mean()}")
                logger.info(f"    Demeaned target value: {t_col_val}")

        design_matrix[col_name] = col_vals
        target_regressors.append(t_col_val)
        column_headers.append(col_name)

    # Now add polynomial terms (each entry is (covariate, power) tuple)
    for cov, power in poly_terms:
        if cov not in continuous_vars:
            raise ValueError(f"Polynomial term '{cov}^{power}' not found in covariates")
        
        vals, t_val, is_categorical = continuous_vars[cov]
        
        if is_categorical:
            raise ValueError(f"Cannot add polynomial term for categorical covariate '{cov}'. Polynomial terms only apply to continuous variables.")
        
        # Use helper function to apply polynomial transform: (x - mean)^power
        powered, t_powered = _process_continuous_term(vals, t_val, power=power)
        
        # Generate column name: Age^2, Age^3, etc.
        col_name = f"{cov}^{power}"
        
        if verbose:
            logger.info(f"  Polynomial term '{col_name}':")
            logger.info(f"    Powered centered target value: {t_powered}")
        
        design_matrix[col_name] = powered
        target_regressors.append(t_powered)
        column_headers.append(col_name)

    # Add intercept LAST (FSL convention)
    design_matrix["intercept"] = np.ones(len(healthy_df))
    target_regressors.append(1.0)
    column_headers.append("intercept")

    
    if verbose:
        logger.info(f"\nFinal design matrix shape: {design_matrix.shape}")
        logger.info(f"Design matrix columns: {design_matrix.columns.tolist()}")


    # Write design matrix in FSL format with column header comments
    num_waves = design_matrix.shape[1]
    num_points = design_matrix.shape[0]
    with open(design_mat_path, "w") as fh:
        fh.write(f"/NumWaves\t{num_waves}\n")
        fh.write(f"/NumPoints\t{num_points}\n")
        fh.write(f"# Columns: {' | '.join(column_headers)}\n")
        fh.write("/Matrix\n")
        for _, row in design_matrix.iterrows():
            fh.write(" ".join(map(str, row.values)) + "\n")

    # Write target regressors with column header comments
    with open(target_vals_path, "w") as fh:
        fh.write(f"# Columns: {' | '.join(column_headers)}\n")
        fh.write(" ".join(map(str, target_regressors)) + "\n")


# ---------------------------------------------------------------------------
# CLI entry point (called from z_score_calc.sh)
# ---------------------------------------------------------------------------

def _cli_generate_design_matrix(args) -> None:
    # Configure logging based on verbose flag
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format='%(message)s',
        stream=sys.stdout
    )
    
    generate_design_matrix(
        demo_path=args.demo,
        target_id=args.target_id,
        healthy_ids_path=args.healthy_ids,
        covariates_str=args.covariates,
        poly_str=args.poly_terms,
        design_mat_path=args.design_mat,
        target_vals_path=args.target_vals,
        verbose=args.verbose,
    )


def main() -> None:

    parser = argparse.ArgumentParser(description="Z-score utility CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    dm = subparsers.add_parser(
        "generate-design-matrix",
        help="Build an FSL GLM design matrix from a demographics CSV.",
    )
    dm.add_argument("--demo", required=True, help="Path to demographics CSV.")
    dm.add_argument("--target-id", required=True, dest="target_id", help="Target subject ID.")
    dm.add_argument("--healthy-ids", required=True, dest="healthy_ids", help="File listing healthy subject IDs.")
    dm.add_argument("--covariates", default="", help="Comma-separated covariate column names.")
    dm.add_argument("--poly-terms", default="", dest="poly_terms", help="Comma-separated polynomial terms formatted as 'covariate_power' (e.g., 'age_2,bmi_3'). Pre-formatted by the processor.")
    dm.add_argument("--design-mat", required=True, dest="design_mat", help="Output path for the design matrix.")
    dm.add_argument("--target-vals", required=True, dest="target_vals", help="Output path for the target regressor values.")
    dm.add_argument("--verbose", action="store_true", help="Enable verbose debug output.")

    args = parser.parse_args()
    if args.command == "generate-design-matrix":
        _cli_generate_design_matrix(args)


if __name__ == "__main__":
    main()
