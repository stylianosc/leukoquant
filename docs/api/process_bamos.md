# leukoquant process-bamos

Run BaMoS (Bayesian Model Selection) white matter hyperintensity (WMH) segmentation.

## Synopsis

```
leukoquant process-bamos [OPTIONS]
```

## Description

`process-bamos` segments white matter hyperintensities from FLAIR (and T1)
images using the BaMoS algorithm.  BaMoS requires a brain parcellation from GIF
to guide the segmentation.

There are two modes of operation, and they can be mixed within the same
batch:

1. **Standalone** - omit `--gif-results-dir` (or a given subject's
   `{subject}` pattern doesn't match anything). GIF is run internally on
   that subject's `--t1` / `--flair` images before BaMoS executes. Requires
   a local GIF installation (bundled, or via `GIF_HOME`) - see
   [`process-gif`](process_gif.md).

2. **Reuse existing GIF output** - provide `--gif-results-dir` pointing to a
   directory containing four files, matched loosely by substring:
   `*_NeuroMorph_Parcellation.nii.gz`, `*_NeuroMorph_Segmentation.nii.gz`,
   `*_NeuroMorph_prior.nii.gz`, and `*TIV*`. These don't need to come from
   this tool's own GIF run - any GIF output with files matching those
   patterns works, including one produced independently by someone else.
   GIF is skipped for that subject and only BaMoS runs. If the directory's
   parcellation file doesn't match the expected naming, `process-bamos`
   fails immediately with a clear error rather than a confusing downstream
   failure.

**Mixed batches**: for a multi-subject run, each subject is resolved
independently - some may have matching `--gif-results-dir` results while
others fall back to standalone mode. A subject with no matching directory
and no available GIF installation is skipped with a warning (not treated as
a fatal error for the whole batch); the rest of the batch still runs. If
*every* subject ends up unprocessable this way, the command fails with a
clear error before submitting anything.

Outputs are written to `{output_dir}/{subject}/bamos/outputs/`.

## Parameters

| Flag | Short | Type | Default | Required | Description |
|------|-------|------|---------|----------|-------------|
| `--subject` | `-i` | `str` | - | **Yes**\* | Subject ID, or path to a `.txt` file listing one subject ID per line |
| `--flair` | `-f` | `str` | - | **Yes**\* | FLAIR NIfTI file or `{subject}` glob pattern |
| `--t1` | `-t` | `str` | - | **Yes**\* | T1 NIfTI file or `{subject}` glob pattern |
| `--gif-results-dir` | `-g` | `str` | `None` | No | Directory containing a prior GIF output, or `{subject}` glob pattern. Subjects without a match fall back to running GIF internally (or are skipped with a warning if no GIF installation is available). |
| `--output-dir` | `-o` | `path` | - | **Yes**\* | Root output directory |
| `--scheduler` | `-s` | `choice` | `local` | No | Execution scheduler: `local` or `sge` |
| `--cores` | `-c` | `int` | `1` | No | Number of Snakemake cores |
| `--config-yaml` | - | `path` | `None` | No | YAML file supplying defaults for any of the above - see [Config file](#config-file-config-yaml) below |
| `--force` | - | flag | `False` | No | Force-rerun BaMoS lesion detection even if outputs already exist |
| `--verbose` | `-v` | flag | `False` | No | Print Snakemake stdout to the console |

\* Required overall, but may be supplied via `--config-yaml` instead of the CLI flag - see below.

### Scheduler options

| Value | Description |
|-------|-------------|
| `local` | Run all jobs on the current machine |
| `sge` | Submit jobs to an SGE cluster |

## Config file (`--config-yaml`)

Instead of (or alongside) CLI flags, `--subject`, `--flair`, `--t1`, and
`--output-dir` may be supplied by a YAML file. **CLI flags always take
priority** over the equivalent YAML key when both are given - the YAML file
only fills in whatever the CLI left unset.

Recognized keys:

| YAML key | Equivalent flag |
|----------|------------------|
| `subject` (or `subject_input`) | `--subject` |
| `flair` (or `flair_pattern`) | `--flair` |
| `t1` (or `t1_pattern`) | `--t1` |
| `gif_results_pattern` | `--gif-results-dir` |
| `output_dir` | `--output-dir` |
| `scheduler` | `--scheduler` |
| `cores` | `--cores` |
| `keep_intermediate` | *(no CLI equivalent - YAML-only)* |

```yaml
# bamos_config.yaml
subject: ./subjects.txt
flair: "./data/{subject}/FLAIR/I*.nii.gz"
t1: "./data/{subject}/T1/I*.nii.gz"
gif_results_pattern: "./gif_outputs/{subject}/gif/outputs"
output_dir: ./outputs
scheduler: sge
```

```bash
leukoquant process-bamos --config-yaml bamos_config.yaml
```

## Returns

`dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `success` | `bool` | `True` if the workflow completed without error |
| `results_dir` | `str` | Absolute path to the output directory |
| `error` | `str` | Error message (only present when `success` is `False`) |

## Notes

- `--gif-results-dir` can be a pattern such as
  `./gif_outputs/{subject}/gif/outputs` - the `{subject}` placeholder is
  substituted per-subject.
- A GIF installation (bundled, or `GIF_HOME`) is only required if at least
  one subject in the batch actually needs GIF run internally. A batch
  where every subject supplies matching `--gif-results-dir` results never
  needs one.
- The corrected lesion mask (`CorrectLesion_*.nii.gz`) is in T1 space and is
  the primary input for `process-metrics` `--lesion-path`.

## Output structure

```
{output_dir}/
└── {subject}/
    └── bamos/
        └── outputs/
            ├── CorrectLesion_*.nii.gz   ← corrected WMH mask (T1 space)
            ├── Lesion_*.nii.gz
            └── ...
```

## Examples

```bash
# Standalone - GIF runs internally
leukoquant process-bamos \
  --subject sub-001 \
  --flair ./data/sub-001/FLAIR/scan.nii.gz \
  --t1    ./data/sub-001/T1/scan.nii.gz \
  --output-dir ./outputs

# Reuse existing GIF output (multi-subject, SGE)
leukoquant process-bamos \
  --subject ./subjects.txt \
  --flair "./data/{subject}/FLAIR/I*.nii.gz" \
  --t1    "./data/{subject}/T1/I*.nii.gz" \
  --gif-results-dir "./gif_outputs/{subject}/gif/outputs" \
  --output-dir ./outputs \
  --scheduler sge
```

## See also

- [`process-gif`](process_gif.md) - produces the GIF parcellation consumed here
- [`process-metrics`](process_metrics.md) - uses the BaMoS lesion mask via `--lesion-path`
- [`process-all`](process_all.md) - runs BaMoS as part of the full pipeline

## Source

[`leukoquant/cli/main.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/cli/main.py) · [`leukoquant/core/bamos_processor.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/core/bamos_processor.py)
