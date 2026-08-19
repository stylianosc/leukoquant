"""
Snakemake workflow for atlas conversion (parcellation relabelling).

Converts one or more brain parcellation files from GIF label space to
FreeSurfer aparc+aseg label space using a CSV mapping file.

Config keys (all required):
  subjects                list of subject IDs
  input_parcellation_map  {subject: host_path}        – host-side input NIfTI per subject
  input_singularity_map   {subject: container_path}   – container-side input NIfTI per subject
  mapping_file            host path of the label-mapping CSV
  mapping_file_singularity  container path of the label-mapping CSV
  output_dir              host root output directory
  output_dir_singularity  container root output directory (default: /output)
  output_subdir           sub-path under {output_dir}/{subject}/ where
                          converted_atlas.mgz is written
                          (e.g. "tracula/intermediate" or "outputs")

Output path: {output_dir}/{subject}/{output_subdir}/converted_atlas.mgz
"""

import os
import re
import sys
from pathlib import Path

LEUKOQUANT_PARENT_DIR = config.get("leukoquant_parent_dir")
CONTAINER_NAME = config.get("container_name", "miniconda_unified_container")
CONTAINER_SIF = os.path.join(
    LEUKOQUANT_PARENT_DIR, "leukoquant/workflow/containers/", CONTAINER_NAME + ".sif")

sys.path.insert(0, LEUKOQUANT_PARENT_DIR)
from leukoquant.utils.container_utils import ensure_container
ensure_container(CONTAINER_SIF)

validate = config.get("validate", True)
verbose  = config.get("verbose",  False)

subjects = config.get("subjects", [])
if isinstance(subjects, str):
    subjects = [subjects]

if not subjects:
    raise ValueError(
        "atlas_converter_workflow.smk requires a non-empty 'subjects' list in config."
    )

input_parcellation_map   = config.get("input_parcellation_map",   {})
input_singularity_map    = config.get("input_singularity_map",    {})
mapping_file             = config.get("mapping_file", "")
mapping_file_singularity = config.get("mapping_file_singularity", "")

OUTPUT_DIR      = config.get("output_dir", "")
OUTPUT_DIR_SING = config.get("output_dir_singularity", "/output")

# Tool name used for standalone runs; callers (e.g. tracula_workflow) override
# output_subdir so that outputs and logs land under their own tool directory.
TOOL_NAME = "atlas_conversion"

# Subdirectory under {output_dir}/{subject}/ where the converted atlas lands.
# Callers set this so the output path matches what downstream tools expect.
# Examples:
#   "tracula/intermediate"    → {output_dir}/{subject}/tracula/intermediate/converted_atlas.mgz
#   "atlas_conversion/outputs"→ {output_dir}/{subject}/atlas_conversion/outputs/converted_atlas.mgz
OUTPUT_SUBDIR = config.get("output_subdir", f"{TOOL_NAME}/outputs")

# Derive the log tool name from the first segment of OUTPUT_SUBDIR so logs
# are grouped with the tool that triggered the conversion.
# Standalone: "atlas_conversion" → {subject}/atlas_conversion/logs/
# From tracula: "tracula"         → {subject}/tracula/logs/
_log_tool = OUTPUT_SUBDIR.split("/")[0]

wildcard_constraints:
    subject="|".join(re.escape(s) for s in subjects),

rule all:
    input:
        [f"{OUTPUT_DIR}/{s}/{OUTPUT_SUBDIR}/converted_atlas.mgz" for s in subjects]

rule convert_atlas:
    input:
        parcellation=lambda wildcards: input_parcellation_map[wildcards.subject],
        mapping=mapping_file,
    output:
        converted=f"{OUTPUT_DIR}/{{subject}}/{OUTPUT_SUBDIR}/converted_atlas.mgz",
    container:
        CONTAINER_SIF,
    params:
        input_singularity=lambda wildcards: input_singularity_map[wildcards.subject],
        mapping_file_singularity=mapping_file_singularity,
        output_singularity=lambda wildcards: (
            f"{OUTPUT_DIR_SING}/{wildcards.subject}/{OUTPUT_SUBDIR}/converted_atlas.mgz"
        ),
        validate_flag="--no-validate" if not validate else "",
        verbose_flag="--verbose" if verbose else "",
        log_file=lambda wildcards: (
            f"{OUTPUT_DIR_SING}/{wildcards.subject}/{_log_tool}/logs/atlas_conversion.log"
        ),
        error_file=lambda wildcards: (
            f"{OUTPUT_DIR_SING}/{wildcards.subject}/{_log_tool}/logs/atlas_conversion_error.log"
        ),
    resources:
        mem_mb=4*1024,
        time="01:00:00",
        name="atlas_conversion",
        workdir=lambda wildcards: f"{OUTPUT_DIR}/{wildcards.subject}/atlas_conversion",
    shell:
        """
        source /leukoquant/leukoquant/utils/bash_utils.sh
        mkdir -p "$(dirname "{params.log_file}")"
        exec > "{params.log_file}"
        exec 2> "{params.error_file}"

        echo "Date: $(date)"
        echo "Date: $(date)" >&2

        export LC_NUMERIC="en_US.UTF-8"
        export PYTHONUNBUFFERED=1

        python /leukoquant/leukoquant/utils/atlas_converter.py \
            --input-parcellation {params.input_singularity} \
            --mapping-file {params.mapping_file_singularity} \
            --output-path {params.output_singularity} \
            {params.validate_flag} \
            {params.verbose_flag}

        echo "Date: $(date)"
        echo "Date: $(date)" >&2
        echo "Atlas conversion completed successfully"
        """
