# Changelog

All notable changes to LeukoQuant are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.0.0] - 2026-10-14

### Changed (metric values differ from 1.x)
- Binary-tract x dilated-lesion statistics summarise the map once; both the fast and the legacy loop multiplied the map in twice.
- Lesion coalescence: a lesion that vanishes on erosion scores 0 (the next size test always overwrote it).
- Atlas streamlines are kept in RAS mm when re-headered to the MGH35 grid (atlas density maps were empty after resampling).
- Stray streamline points up to 0.01 voxel outside the field of view are clamped instead of failing the density map.

### Added
- Skewness and excess kurtosis for every region summary (tract, skeleton, lesion, penumbra, whole brain, Z-score maps); not used in the lesion-versus-outside change ratios.
- Faster metrics: sparse-support tract statistics, voxel-based fractal dimension, vectorised surface area, batch streamline transform.
- Memory requests sized from the data: metrics (from the T1 headers, with `LEUKOQUANT_METRICS_MEM_GB` override), z-score GLM and TRACULA; peak memory and scratch use reported by every step; out-of-memory and other step failures reported loudly.
- DWI conversion fallback chain (dcm2niix, non-spatial dcm2niix, mri_convert) with input validation; sessions with unusable input are skipped by `process-all` (`--retry-input-failures` to retry).
- TRACULA brainstem field-of-view guard and recovery; missing tracts fall back to the atlas in the metrics.
- `process-all --force` takes a list of stages; each run gets its own Snakemake config, so concurrent runs in one output folder no longer collide.
- 4-D T1 inputs are collapsed to their first frame.

### Fixed
- Containers get a real writable `/tmp`; one metrics scratch folder per task.
- `--rerun-triggers` no longer swallows the positional target in standalone processors.
- Requires snakemake-executor-plugin-sge 0.6.31.

## [1.0.0] - 2026-05-29

### Added

**Core pipeline commands**
- `process-gif` - GIF brain segmentation and parcellation (T1 and FLAIR inputs)
- `process-recon` - FreeSurfer `recon-all` cortical surface reconstruction
- `process-bamos` - BaMoS white matter hyperintensity (WMH) segmentation
- `process-dti` - DTI model fitting (FA, MD, AD, RD maps)
- `process-noddi` - NODDI model fitting (NDI, ODI, ISOVF maps) via AMICO
- `process-tracula` - TRACULA probabilistic tractography (FreeSurfer and GIF atlas variants)
- `process-atlas-conversion` - GIF label atlas conversion to FreeSurfer space
- `process-metrics` - Lesion-informed tract metrics and WMH region statistics
- `process-tract-qc` - Tractography quality control report
- `process-zscore` - Z-score normalisation of diffusion metrics against a healthy cohort
- `process-all` - Orchestrated end-to-end pipeline combining all modules above

**Infrastructure**
- Snakemake-based workflow engine with SGE cluster and local execution modes
- Singularity/Apptainer container support for all compute steps
- `setup_platform.sh` - automated Apptainer installation (Linux and macOS via Lima VM)
- Automated container download from Hugging Face Hub on first run
- Multi-subject batch execution via subject list files with `{subject}` path templating

**Metrics**
- Whole-brain WMH metrics CSV (volume, count, lesion-level DTI/NODDI statistics)
- Tract-level metrics CSV (per-tract lesion burden, FA/MD/NDI/ODI within tracts)

**Documentation**
- Full CLI reference for all 11 commands (`--help` and API docs)
- Installation and platform setup guide
- CLI usage examples (Wave 1 and Wave 2 workflows)
- Output directory structure reference
- GitHub Pages site at https://stylianosc.github.io/leukoquant/

[1.0.0]: https://github.com/stylianosc/leukoquant/releases/tag/v1.0.0
