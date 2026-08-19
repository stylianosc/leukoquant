# leukoquant process-tracula

Run TRACULA (TRActs Constrained by UnderLying Anatomy) probabilistic tractography.

## Synopsis

```
leukoquant process-tracula [OPTIONS]
```

## Description

`process-tracula` runs the TRACULA probabilistic tractography pipeline, which
reconstructs major white matter pathways constrained by a brain parcellation atlas.

TRACULA requires:

1. **DWI data** (`--dwi`) - the diffusion-weighted images.
2. **A brain parcellation** - provided by either FreeSurfer `recon-all` or GIF,
   selected via `--parcellation`.
3. **A FreeSurfer reconstruction** - either pre-computed (`--freesurfer-recon-dir`)
   or run automatically when `--t1` is provided and `--parcellation` includes `freesurfer`.

### Parcellation modes

| `--parcellation` | Brain atlas used | FreeSurfer recon required |
|------------------|-----------------|--------------------------|
| `freesurfer` (default) | `aparc+aseg.mgz` from `recon-all` | Yes |
| `gif` | GIF multi-label parcellation, converted to FreeSurfer label space | Yes (for registration) |
| `freesurfer,gif` | Both parcellations; two independent tract sets produced | Yes |

For `--parcellation gif`, the GIF NIfTI is resolved from
`{output_dir}/{subject}/gif/outputs/*_NeuroMorph_Parcellation.nii.gz` automatically
(when `process-gif` output exists in `--output-dir`), or supplied explicitly via
`--brain-parcellation`.

Outputs are written to `{output_dir}/{subject}/tracula-{parcellation}/outputs/`.

## Parameters

| Flag | Short | Type | Default | Required | Description |
|------|-------|------|---------|----------|-------------|
| `--subject` | `-s` | `str` | - | **Yes**\*\* | Subject ID, or path to a `.txt` file listing one subject ID per line |
| `--dwi` | - | `str` | - | **Yes**\*\* | DWI NIfTI (or `.zip` DICOM archive) file or `{subject}` glob pattern |
| `--bvecs` | - | `str` | `None` | No | bvecs file or `{subject}` glob pattern |
| `--bvals` | - | `str` | `None` | No | bvals file or `{subject}` glob pattern |
| `--t1` | - | `str` | `None` | No* | T1 NIfTI file or `{subject}` glob pattern. Required when recon-all has not been run yet. |
| `--freesurfer-recon-dir` | - | `path` | `None` | No | Root directory containing per-subject recon-all outputs (`{dir}/{subject}/recon-all/outputs`). If omitted or missing, recon-all runs automatically (requires `--t1`). |
| `--parcellation` | - | `str` | `freesurfer` | No | Parcellation(s) to use. Comma-separated for multiple (e.g. `freesurfer,gif`). |
| `--brain-parcellation` | - | `path` | `None` | No | Explicit path to a brain parcellation NIfTI. For `gif`, omit to auto-resolve from GIF output directory. |
| `--output-dir` | `-o` | `path` | - | **Yes**\*\* | Root output directory for TRACULA results |
| `--scratch` | - | `path` | `None` | No | Optional scratch directory for bedpostx temporary files |
| `--scheduler` | - | `choice` | `local` | No | Execution scheduler: `local` or `sge` |
| `--cores` | `-c` | `int` | `1` | No | Number of Snakemake cores |
| `--config-yaml` | - | `path` | `None` | No | YAML file supplying defaults for any of the above - see [Config file](#config-file-config-yaml) below |
| `--force` | - | flag | `False` | No | Force-rerun TRACULA tractography (for `--parcellation`) even if outputs already exist |
| `--verbose` | `-v` | flag | `False` | No | Print Snakemake stdout to the console |

\* `--t1` is required when `recon-all` outputs are not already available.
\*\* Required overall, but may be supplied via `--config-yaml` instead of the CLI flag - see below.

### Parcellation options

| Value | Description |
|-------|-------------|
| `freesurfer` | FreeSurfer `aparc+aseg` atlas (default) |
| `gif` | GIF multi-label parcellation, auto-resolved or via `--brain-parcellation` |
| `freesurfer,gif` | Both parcellations in one run; produces two output subfolders |

### Scheduler options

| Value | Description |
|-------|-------------|
| `local` | Run all jobs on the current machine |
| `sge` | Submit jobs to an SGE cluster |

## Config file (`--config-yaml`)

Instead of (or alongside) CLI flags, `--subject`, `--dwi`, and `--output-dir`
may be supplied by a YAML file. **CLI flags always take priority** over the
equivalent YAML key when both are given - the YAML file only fills in
whatever the CLI left unset.

Recognized keys:

| YAML key | Equivalent flag |
|----------|------------------|
| `subject` (or `subject_input`) | `--subject` |
| `dwi` (or `dwi_pattern`) | `--dwi` |
| `bvecs` (or `bvecs_pattern`) | `--bvecs` |
| `bvals` (or `bvals_pattern`) | `--bvals` |
| `t1` (or `t1_pattern`) | `--t1` |
| `freesurfer_recon_dir` | `--freesurfer-recon-dir` |
| `parcellation` | `--parcellation` |
| `brain_parcellation` | `--brain-parcellation` |
| `output_dir` | `--output-dir` |
| `scheduler` | `--scheduler` |
| `cores` | `--cores` |
| `mapping_file` | *(no CLI equivalent - YAML-only; required for non-freesurfer, non-gif custom parcellations)* |
| `keep_intermediate` | *(no CLI equivalent - YAML-only)* |

```yaml
# tracula_config.yaml
subject: ./subjects.txt
dwi: "./data/{subject}/DWI/I*.zip"
freesurfer_recon_dir: "./recon_outputs/{subject}/recon-all/outputs"
parcellation: freesurfer
output_dir: ./tracula_outputs
scheduler: sge
```

```bash
leukoquant process-tracula --config-yaml tracula_config.yaml
```

## Returns

`dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `success` | `bool` | `True` if the workflow completed without error |
| `results_dir` | `str` | Absolute path to the output directory |
| `error` | `str` | Error message (only present when `success` is `False`) |

## Notes

- Each parcellation produces a separate output subfolder:
  `tracula-freesurfer/` and/or `tracula-gif/`.
- When `--parcellation gif` is used for multiple subjects, set `--output-dir`
  to the same directory where `process-gif` wrote its outputs so the GIF NIfTI
  is auto-resolved per subject without needing `--brain-parcellation`.
- For GIF parcellation, the GIF brain extraction mask (`*_NeuroMorph_Brain.nii.gz`)
  is automatically applied to `T1.mgz` to produce `brain.mgz` and `brainmask.mgz`,
  replacing the recon-all versions. Both files are identical (the GIF-masked T1).
- TRACULA is computationally intensive (bedpostx + path sampling).
  Use `--scheduler sge` for any dataset with more than one subject.

## Output structure

```
{output_dir}/
└── {subject}/
    ├── tracula-freesurfer/
    │   └── outputs/
    │       ├── dpath/
    │       │   ├── lh.cst_AS/
    │       │   │   └── path.pd.nii.gz
    │       │   └── ...
    │       └── merged_avg16_syn_bbr.mgz
    └── tracula-gif/           ← only when --parcellation gif or freesurfer,gif
        └── outputs/
            └── ...
```

## Examples

```bash
# Single subject, auto recon-all, freesurfer parcellation (default)
leukoquant process-tracula \
  --subject sub-001 \
  --dwi ./data/sub-001/DWI/data.nii.gz \
  --t1  ./data/sub-001/T1/scan.nii.gz \
  --output-dir ./outputs

# Multiple subjects, reuse existing recon-all, explicit freesurfer parcellation
leukoquant process-tracula \
  --subject ./subjects.txt \
  --dwi "./data/{subject}/DWI/I*.zip" \
  --freesurfer-recon-dir "./recon_outputs/{subject}/recon-all/outputs" \
  --parcellation freesurfer \
  --output-dir ./tracula_outputs \
  --scheduler sge

# Single subject, GIF parcellation (explicit brain parcellation path)
leukoquant process-tracula \
  --subject sub-001 \
  --dwi "./data/sub-001/DWI/I*.zip" \
  --freesurfer-recon-dir "./recon_outputs/sub-001/recon-all/outputs" \
  --brain-parcellation "./gif_outputs/sub-001/gif/outputs/sub001_NeuroMorph_Parcellation.nii.gz" \
  --parcellation gif \
  --output-dir ./gif_outputs \
  --scheduler sge

# Multiple subjects, GIF parcellation (auto-resolved from gif_outputs)
leukoquant process-tracula \
  --subject ./subjects.txt \
  --dwi "./data/{subject}/DWI/I*.zip" \
  --freesurfer-recon-dir "./recon_outputs/{subject}/recon-all/outputs" \
  --parcellation gif \
  --output-dir ./gif_outputs \
  --scheduler sge

# Multiple subjects, both freesurfer and gif in one run
leukoquant process-tracula \
  --subject ./subjects.txt \
  --dwi "./data/{subject}/DWI/I*.zip" \
  --freesurfer-recon-dir "./recon_outputs/{subject}/recon-all/outputs" \
  --parcellation freesurfer,gif \
  --output-dir ./gif_outputs \
  --scheduler sge
```

## See also

- [`process-gif`](process_gif.md) - produces the GIF parcellation for `--parcellation gif`
- [`process-recon`](process_recon.md) - produces the recon-all outputs for `--freesurfer-recon-dir`
- [`process-metrics`](process_metrics.md) - computes metrics on TRACULA tract outputs
- [`process-tract-qc`](process_tract_qc.md) - quality-controls TRACULA tract outputs
- [`process-all`](process_all.md) - runs TRACULA as part of the full pipeline

## Source

[`leukoquant/cli/main.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/cli/main.py) · [`leukoquant/core/tracula_processor.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/core/tracula_processor.py)
