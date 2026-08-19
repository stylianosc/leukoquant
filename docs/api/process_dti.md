# leukoquant process-dti

Run DTI (Diffusion Tensor Imaging) model fitting.

## Synopsis

```
leukoquant process-dti [OPTIONS]
```

## Description

`process-dti` fits the diffusion tensor model to multi-shell DWI data,
producing voxel-wise parametric maps including fractional anisotropy (FA),
mean diffusivity (MD), radial diffusivity (RD), and axial diffusivity (AD).

The command accepts a single DWI file or a glob pattern (e.g. a `.zip` archive
of DICOM images).  Gradient directions (`bvecs`) and b-values (`bvals`) can
be provided explicitly; if omitted, they are expected to be co-located with
the DWI file.

An optional brain mask can be provided via `--mask`.  If no mask is given,
fitting is performed on the whole image volume.  Alternatively, `--skull-strip`
auto-generates a mask using FreeSurfer `mri_synthstrip` before fitting.

Outputs are written to `{output_dir}/{subject}/dti/outputs/`.

## Parameters

| Flag | Short | Type | Default | Required | Description |
|------|-------|------|---------|----------|-------------|
| `--subject` | `-i` | `str` | - | **Yes**\* | Subject ID, or path to a `.txt` file listing one subject ID per line |
| `--dwi` | - | `str` | - | **Yes**\* | DWI NIfTI (or `.zip` DICOM archive) file or `{subject}` glob pattern |
| `--bvecs` | - | `str` | `None` | No | bvecs file or `{subject}` glob pattern |
| `--bvals` | - | `str` | `None` | No | bvals file or `{subject}` glob pattern |
| `--mask` | - | `str` | `None` | No | Brain mask NIfTI or `{subject}` glob pattern. If omitted, fits on the whole image. |
| `--skull-strip` | - | flag | `False` | No | Auto-generate a brain mask using `mri_synthstrip` before fitting. Ignored when `--mask` is provided. |
| `--output-dir` | `-o` | `path` | - | **Yes**\* | Root output directory for DTI results |
| `--scheduler` | `-s` | `choice` | `local` | No | Execution scheduler: `local` or `sge` |
| `--cores` | `-c` | `int` | `1` | No | Number of Snakemake cores |
| `--config-yaml` | - | `path` | `None` | No | YAML file supplying defaults for any of the above - see [Config file](#config-file-config-yaml) below |
| `--force` | - | flag | `False` | No | Force-rerun DTI fitting even if outputs already exist |
| `--verbose` | `-v` | flag | `False` | No | Print Snakemake stdout to the console |

\* Required overall, but may be supplied via `--config-yaml` instead of the CLI flag - see below.

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
| `mask_pattern` | `--mask` |
| `skull_strip` | `--skull-strip` |
| `output_dir` | `--output-dir` |
| `scheduler` | `--scheduler` |
| `cores` | `--cores` |
| `keep_intermediate` | *(no CLI equivalent - YAML-only)* |

```yaml
# dti_config.yaml
subject: ./subjects.txt
dwi: "./data/{subject}/DWI/I*.zip"
output_dir: ./outputs
scheduler: sge
```

```bash
leukoquant process-dti --config-yaml dti_config.yaml
```

## Returns

`dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `success` | `bool` | `True` if the workflow completed without error |
| `results_dir` | `str` | Absolute path to the output directory |
| `error` | `str` | Error message (only present when `success` is `False`) |

## Notes

- DTI maps (FA, MD, etc.) are in DWI space.  When used as `--maps` inputs to
  `process-metrics`, set `space` to `dwi` in the path spec.
- If `--skull-strip` is used, `mri_synthstrip` must be available inside the
  container image.

## Output structure

```
{output_dir}/
└── {subject}/
    └── dti/
        └── outputs/
            ├── fa.nii.gz    ← fractional anisotropy
            ├── md.nii.gz    ← mean diffusivity
            ├── rd.nii.gz    ← radial diffusivity
            ├── ad.nii.gz    ← axial diffusivity
            └── ...
```

## Examples

```bash
# Single subject, no mask
leukoquant process-dti \
  --subject sub-001 \
  --dwi ./data/sub-001/DWI/data.nii.gz \
  --output-dir ./outputs

# Multiple subjects, auto skull-strip, SGE
leukoquant process-dti \
  --subject ./subjects.txt \
  --dwi "./data/{subject}/DWI/I*.zip" \
  --skull-strip \
  --output-dir ./outputs \
  --scheduler sge

# Explicit brain mask
leukoquant process-dti \
  --subject sub-001 \
  --dwi ./data/sub-001/DWI/data.nii.gz \
  --mask ./data/sub-001/mask.nii.gz \
  --output-dir ./outputs
```

## See also

- [`process-noddi`](process_noddi.md) - NODDI fitting on the same DWI data
- [`process-metrics`](process_metrics.md) - uses DTI maps via `--maps`
- [`process-all`](process_all.md) - runs DTI fitting as part of the full pipeline

## Source

[`leukoquant/cli/main.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/cli/main.py) · [`leukoquant/core/dti_processor.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/core/dti_processor.py)
