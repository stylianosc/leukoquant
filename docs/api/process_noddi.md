# leukoquant process-noddi

Run NODDI (Neurite Orientation Dispersion and Density Imaging) model fitting.

## Synopsis

```
leukoquant process-noddi [OPTIONS]
```

## Description

`process-noddi` fits the NODDI microstructure model to multi-shell DWI data,
producing voxel-wise parametric maps including:

- **ODI** - Orientation Dispersion Index (neurite orientation variability)
- **ICVF** - Intra-Cellular Volume Fraction (neurite density)
- **ISOVF** - Isotropic Volume Fraction (free water fraction)

The command accepts the same DWI input conventions as `process-dti`.
A brain mask suppresses fitting outside brain tissue and substantially reduces
compute time; use `--skull-strip` if no mask is available.

Outputs are written to `{output_dir}/{subject}/noddi/outputs/`.

## Parameters

| Flag | Short | Type | Default | Required | Description |
|------|-------|------|---------|----------|-------------|
| `--subject` | `-i` | `str` | - | **Yes**\* | Subject ID, or path to a `.txt` file listing one subject ID per line |
| `--dwi` | - | `str` | - | **Yes**\* | DWI NIfTI (or `.zip` DICOM archive) file or `{subject}` glob pattern |
| `--bvecs` | - | `str` | `None` | No | bvecs file or `{subject}` glob pattern |
| `--bvals` | - | `str` | `None` | No | bvals file or `{subject}` glob pattern |
| `--mask` | - | `str` | `None` | No | Brain mask NIfTI or `{subject}` glob pattern. If omitted, fits on the whole image. |
| `--skull-strip` | - | flag | `False` | No | Auto-generate a brain mask using `mri_synthstrip` before fitting. Ignored when `--mask` is provided. |
| `--output-dir` | `-o` | `path` | - | **Yes**\* | Root output directory for NODDI results |
| `--scheduler` | `-s` | `choice` | `local` | No | Execution scheduler: `local` or `sge` |
| `--cores` | `-c` | `int` | `1` | No | Number of Snakemake cores |
| `--config-yaml` | - | `path` | `None` | No | YAML file supplying defaults for any of the above - see [Config file](#config-file-config-yaml) below |
| `--force` | - | flag | `False` | No | Force-rerun NODDI fitting even if outputs already exist |
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
# noddi_config.yaml
subject: ./subjects.txt
dwi: "./data/{subject}/DWI/I*.zip"
output_dir: ./outputs
scheduler: sge
cores: 4
```

```bash
leukoquant process-noddi --config-yaml noddi_config.yaml
```

## Returns

`dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `success` | `bool` | `True` if the workflow completed without error |
| `results_dir` | `str` | Absolute path to the output directory |
| `error` | `str` | Error message (only present when `success` is `False`) |

## Notes

- NODDI fitting is computationally intensive.  Using `--scheduler sge` and
  increasing `--cores` is strongly recommended for multi-subject datasets.
- Outputs are in DWI space.  Specify `space=dwi` when passing these maps to
  `process-metrics` via `--maps`.
- Multi-shell DWI data (at least two non-zero b-values) is required for NODDI
  fitting.

## Output structure

```
{output_dir}/
└── {subject}/
    └── noddi/
        └── outputs/
            ├── odi.nii.gz     ← orientation dispersion index
            ├── icvf.nii.gz    ← intra-cellular volume fraction
            ├── isovf.nii.gz   ← isotropic volume fraction
            └── ...
```

## Examples

```bash
# Single subject, auto skull-strip
leukoquant process-noddi \
  --subject sub-001 \
  --dwi ./data/sub-001/DWI/data.nii.gz \
  --skull-strip \
  --output-dir ./outputs

# Multiple subjects, explicit mask, SGE cluster
leukoquant process-noddi \
  --subject ./subjects.txt \
  --dwi  "./data/{subject}/DWI/I*.zip" \
  --mask "./data/{subject}/mask.nii.gz" \
  --output-dir ./outputs \
  --scheduler sge --cores 4
```

## See also

- [`process-dti`](process_dti.md) - DTI fitting on the same DWI data
- [`process-metrics`](process_metrics.md) - uses NODDI maps via `--maps`
- [`process-all`](process_all.md) - runs NODDI fitting as part of the full pipeline

## Source

[`leukoquant/cli/main.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/cli/main.py) · [`leukoquant/core/noddi_processor.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/core/noddi_processor.py)
