# API Reference

Complete parameter reference for all LeukoQuant CLI commands.

Each command wraps a Snakemake workflow that runs inside a Singularity/Apptainer container.
Commands can be executed locally or dispatched to an SGE cluster via `--scheduler sge`.

---

## Command index

Commands are listed in typical pipeline order.  Most can also be run as
standalone tools independent of the rest of the pipeline.

### Preprocessing

| Command | Description |
|---------|-------------|
| [`process-gif`](process_gif.md) | GIF brain segmentation and multi-label parcellation |
| [`process-recon`](process_recon.md) | FreeSurfer recon-all cortical reconstruction |
| [`process-bamos`](process_bamos.md) | BaMoS white matter hyperintensity (WMH) segmentation |
| [`process-dti`](process_dti.md) | DTI model fitting - FA, MD, RD, AD maps |
| [`process-noddi`](process_noddi.md) | NODDI model fitting - ODI, ICVF, ISOVF maps |
| [`process-atlas-conversion`](process_atlas_conversion.md) | Convert GIF parcellation labels to FreeSurfer label space |

### Tractography

| Command | Description |
|---------|-------------|
| [`process-tracula`](process_tracula.md) | TRACULA probabilistic tractography (freesurfer and/or gif parcellation) |

### Metrics & QC

| Command | Description |
|---------|-------------|
| [`process-metrics`](process_metrics.md) | Extract lesion-informed metrics along white matter tracts |
| [`process-tract-qc`](process_tract_qc.md) | Tractography quality control report |
| [`process-zscore`](process_zscore.md) | Voxel-wise Z-score maps relative to a healthy cohort |

### Full pipeline

| Command | Description |
|---------|-------------|
| [`process-all`](process_all.md) | Run all pipeline stages in the correct order |

---

## Common conventions

### Subject input

Every command accepts `--subject` as either:

- A single subject ID string (e.g. `sub-001`)
- A path to a plain-text file with one subject ID per line

### Path patterns

File path arguments support the `{subject}` placeholder, which is substituted
per-subject at runtime:

```
./data/{subject}/T1/I*.nii.gz   →   ./data/sub-001/T1/I001.nii.gz
                                     ./data/sub-002/T1/I002.nii.gz
                                     ...
```

Glob wildcards (`*`, `?`) are resolved on the host filesystem before the
pattern is passed to the container.

### Composite path format

Commands that take tractography or metric maps use a three-part colon-separated
spec:

```
base_dir:relative_glob:space
```

| Part | Description |
|------|-------------|
| `base_dir` | Root directory (may contain `{subject}`) |
| `relative_glob` | Glob relative to `base_dir` (or `base_dir/{subject}`) |
| `space` | Coordinate space: `t1`, `dwi`, or `atlas` |

Example: `./tracula_out/{subject}/tracula-freesurfer/outputs:dpath/*/path.pd.nii.gz:dwi`

### Config file (`--config-yaml`)

Every command accepts an optional `--config-yaml <path>` pointing to a YAML
file that supplies default values for its arguments. This is convenient for
complex, multi-flag invocations you want to save and re-run, or for
CI/scripted usage where writing a file is cleaner than a long command line.

**Precedence: CLI flags always win.** A `--config-yaml` value is only used
for an argument that was *not* also given on the command line - it never
overrides an explicit CLI flag. This lets you keep a shared config file for a
dataset and override just one or two values per invocation:

```bash
leukoquant process-metrics --config-yaml my_config.yaml --output-dir ./run2
```

YAML keys generally match the long-form flag name with dashes replaced by
underscores (e.g. `--tractography-path` → `tractography_path`). Most
commands also accept a shorter alias matching the flag itself (e.g. `t1` as
well as `t1_pattern`) - each command's page lists its exact recognized keys
under its own "Config file" section.

A `--config-yaml` path that doesn't exist raises an error immediately; a
command is never silently run with an unintended file.

### Returns

All commands return a Python `dict` with at minimum:

| Key | Type | Description |
|-----|------|-------------|
| `success` | `bool` | `True` if the workflow completed without error |
| `results_dir` | `str` | Absolute path to the output directory |
| `error` | `str` | Error message (only present when `success` is `False`) |
