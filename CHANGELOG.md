# Changelog

All notable changes to LeukoQuant are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.1.1] - 2026-10-28

### Changed
- Tract-to-lesion distance of a subject without any lesion: the nearest and the mean distance are now NaN (2.1.0 reported the 100 mm cap, which clipped every distance at 100 mm). The fraction within 5 mm stays 0, and the NAWM and the outer penumbra rings of a lesion-free subject are unchanged (the whole tract and empty). No distance was clipped and no participant was lesion-free in the three clinical cohorts of the paper, so none of its values change.

---

## [2.1.0] - 2026-10-22

### Changed (metric values and columns differ from 2.0.0)
- **Skeleton statistics** (`skeleton_*`, including PSMD and the mean skeleton FA / MD) are statistics of the map on a binary skeleton. Before, the thin binary skeleton was resampled to T1 space with linear interpolation and its fractional values multiplied the map. The skeleton is now resampled with nearest neighbour and thresholded at 0.5.
- **Lesion-versus-outside change**: the `(inside - outside) / (inside + outside)` ratios (`*_change`, `*_change_absolute`) are replaced by `*_diff` (inside minus reference, 10 statistics) and `auc` (probability that a lesion voxel exceeds a reference voxel) per map, computed against the tract's normal-appearing white matter (NAWM). The ratio is unbounded for signed maps such as Z-score maps (11-18% of values lay outside [-1, 1]).
- Region summaries report `p5` and `p95` instead of `min` and `max`, plus `n_voxels`; empty regions stay 0.0.
- Statistics of a map weighted by tract or lesion probability are probability-weighted (weighted mean, standard deviation, percentiles, shape) instead of statistics of map x probability.
- Penumbra and dilated lesion are defined by distance (3 mm), not by three voxel dilations.
- Lesion components are 26-connected everywhere (the tract metrics used 6-connectivity).
- Per-lesion: `merged_after_dilation` counts distinct lesions (it was always true); coalescence size bands are percentages of TIV (only the scores 0 and 1 could occur); fractal dimension uses box sizes that fit the lesion; surface area (and compactness, sphericity) from a mesh instead of exposed voxel faces; the three best-matching tracts are reported.
- Tract-to-lesion distance: nearest, mean and fraction within 5 mm replace the distance of the tract's centre of mass; the distance of a subject without lesion is the 100 mm cap (it was a meaningless border distance).
- `volume_tiv_normalised_mm3` removed (it equals `load_pct` / 100).
- NaN / Inf voxels in inputs are treated as missing data (set to 0, logged); the dense fallback implementations are removed.

### Added
- Penumbra rings 3-6 and 6-9 mm and the NAWM region with their statistics.
- Per-lesion `n_lesions_merged_after_dilation`, `max_inscribed_radius_mm`.
- Tract: share of the tract's streamline density inside lesions.
- Whole brain: lesion count, mean lesion volume; periventricular and deep lesion volumes and share (from `--parcellation-path` and `--ventricle-labels`).
- `metrics_output_dir` (process-all configuration): the metrics step writes its folders there, nothing existing is overwritten.
- Force guard for `process-all --force`: the run is planned first and refused when a stage upstream of the forced ones would also rerun (typically through newer modification times after a fresh clone); `--dry-run` prints the per-rule job counts and writes Snakemake's full plan with the reason of every job (`dry_run_plan.txt`), `--allow-upstream-rerun` overrides the refusal, and `--only-forced` plans and runs only the forced stages whatever the modification times say.
- The metrics output root is mounted in the container (`/metrics_output`) when `metrics_output_dir` is set.

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
