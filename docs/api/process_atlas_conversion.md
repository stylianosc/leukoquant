# leukoquant process-atlas-conversion

Convert a GIF brain parcellation to FreeSurfer label space.

## Synopsis

```
leukoquant process-atlas-conversion [OPTIONS]
```

## Description

`process-atlas-conversion` re-labels a GIF multi-label parcellation NIfTI
to match the FreeSurfer `aparc+aseg` label scheme using a user-supplied CSV
mapping file.  The converted atlas (`converted_atlas.mgz`) can then be used
as the `--brain-parcellation` input to `process-tracula` for non-standard
parcellation schemes.

The mapping CSV must have at least two columns: the source GIF label integer
and the target FreeSurfer label integer.  Optional post-conversion label
validation checks that all expected labels are present in the output.

Outputs are written to `{output_dir}/{subject}/atlas_conversion/outputs/`.

## Parameters

| Flag | Short | Type | Default | Required | Description |
|------|-------|------|---------|----------|-------------|
| `--subject` | `-i` | `str` | - | **Yes**\* | Subject ID (used as the output subdirectory name) |
| `--input-parcellation` | `-p` | `path` | - | **Yes**\* | Input GIF parcellation NIfTI file (must exist) |
| `--mapping-file` | `-m` | `path` | - | **Yes**\* | CSV label mapping file (must exist). Columns: source label, target label. |
| `--output-dir` | `-o` | `path` | `None` | No | Root output directory. Defaults to the parent directory of `--input-parcellation`. |
| `--no-validate` | - | flag | `False` | No | Skip post-conversion label validation |
| `--scheduler` | `-s` | `choice` | `local` | No | Execution scheduler: `local` or `sge` |
| `--cores` | `-c` | `int` | `1` | No | Number of Snakemake cores |
| `--config-yaml` | - | `path` | `None` | No | YAML file supplying defaults for any of the above - see [Config file](#config-file-config-yaml) below |
| `--force` | - | flag | `False` | No | Force-rerun atlas conversion even if outputs already exist |
| `--verbose` | `-v` | flag | `False` | No | Print Snakemake stdout to the console |

\* Required overall, but may be supplied via `--config-yaml` instead of the CLI flag - see below.

### Scheduler options

| Value | Description |
|-------|-------------|
| `local` | Run all jobs on the current machine |
| `sge` | Submit jobs to an SGE cluster |

## Config file (`--config-yaml`)

Instead of (or alongside) CLI flags, `--subject`, `--input-parcellation`, and
`--mapping-file` may be supplied by a YAML file. **CLI flags always take
priority** over the equivalent YAML key when both are given - the YAML file
only fills in whatever the CLI left unset.

Recognized keys:

| YAML key | Equivalent flag |
|----------|------------------|
| `subject` | `--subject` |
| `input_parcellation` | `--input-parcellation` |
| `mapping_file` | `--mapping-file` |
| `output_dir` | `--output-dir` |
| `scheduler` | `--scheduler` |
| `cores` | `--cores` |

```yaml
# atlas_conversion_config.yaml
subject: sub-001
input_parcellation: ./gif_outputs/sub-001/gif/outputs/sub001_NeuroMorph_Parcellation.nii.gz
mapping_file: ./leukoquant/mappings/gif_to_freesurfer.csv
output_dir: ./atlas_conversion_outputs
```

```bash
leukoquant process-atlas-conversion --config-yaml atlas_conversion_config.yaml
```

## Returns

`dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `success` | `bool` | `True` if the workflow completed without error |
| `results_dir` | `str` | Absolute path to the output directory |
| `error` | `str` | Error message (only present when `success` is `False`) |

## Notes

- `--subject` is used only as the output subdirectory name; the input
  parcellation file is specified directly via `--input-parcellation` and is not
  glob-expanded.
- The built-in GIF→FreeSurfer mapping is provided at
  `leukoquant/mappings/gif_to_freesurfer.csv`.
- When `--output-dir` is omitted, the converted atlas is written alongside
  the input parcellation file.

## Output structure

```
{output_dir}/
└── {subject}/
    └── atlas_conversion/
        └── outputs/
            └── converted_atlas.mgz
```

## Examples

```bash
# Convert GIF parcellation to FreeSurfer labels using the built-in mapping
leukoquant process-atlas-conversion \
  --subject sub-001 \
  --input-parcellation ./gif_outputs/sub-001/gif/outputs/sub001_NeuroMorph_Parcellation.nii.gz \
  --mapping-file ./leukoquant/mappings/gif_to_freesurfer.csv \
  --output-dir ./atlas_conversion_outputs

# Skip label validation
leukoquant process-atlas-conversion \
  --subject sub-001 \
  --input-parcellation ./gif_outputs/sub-001/gif/outputs/sub001_NeuroMorph_Parcellation.nii.gz \
  --mapping-file ./leukoquant/mappings/gif_to_freesurfer.csv \
  --no-validate \
  --output-dir ./atlas_conversion_outputs
```

## See also

- [`process-gif`](process_gif.md) - produces the input parcellation NIfTI
- [`process-tracula`](process_tracula.md) - uses the converted atlas via `--brain-parcellation`
- [`process-all`](process_all.md) - runs atlas conversion as part of the full pipeline

## Source

[`leukoquant/cli/main.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/cli/main.py) · [`leukoquant/core/atlas_converter_processor.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/core/atlas_converter_processor.py)
