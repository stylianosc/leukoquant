# Changelog

All notable changes to LeukoQuant are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

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
