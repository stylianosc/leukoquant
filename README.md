# LeukoQuant

[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C%203.10-blue.svg)](https://www.python.org/)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue.svg)](https://stylianosc.github.io/leukoquant/)
[![Cite](https://img.shields.io/badge/cite-CITATION.cff-green.svg)](CITATION.cff)

**Lesion-informed white matter damage metrics toolkit for cerebral small vessel disease.**

LeukoQuant provides a complete, reproducible pipeline from raw MRI scans to quantitative white matter metrics. It wraps Snakemake workflows inside a clean CLI, supports local and SGE cluster execution, and handles all container-based software dependencies automatically.

> **Documentation:** [https://stylianosc.github.io/leukoquant/](https://stylianosc.github.io/leukoquant/)

> **Status:** LeukoQuant is currently available as a pre-release to support open science initiatives. A manuscript detailing the core methodology and our comprehensive clinical validation across the ADNI, OASIS Brain, and EPAD cohorts is currently in preparation.
>
> **Citation:** If you use this software in the meantime, please cite it via our Zenodo DOI: [Insert DOI]. If you are planning a large-scale clinical application, please contact the authors as our primary validation manuscript is imminent.

---

## Pipeline overview

```
Raw MRI (T1, FLAIR, DWI)
        │
        ├─ process-gif          → GIF brain segmentation & parcellation
        ├─ process-recon        → FreeSurfer recon-all
        ├─ process-bamos        → BaMoS white matter hyperintensity segmentation
        ├─ process-dti          → DTI fitting (FA, MD, AD, RD)
        ├─ process-noddi        → NODDI fitting (NDI, ODI, ISOVF)
        ├─ process-atlas-conversion → Convert GIF labels → FreeSurfer space
        │
        ├─ process-tracula      → TRACULA probabilistic tractography
        │
        ├─ process-metrics      → Lesion-informed tract metrics CSV
        ├─ process-tract-qc     → Tractography quality control report
        ├─ process-zscore       → Z-score maps (target vs. healthy cohort)
        │
        └─ process-all          → Full pipeline (all of the above, ordered)
```

---

## Prerequisites

| Requirement             | Check                   | Notes                                                                                                                                       |
| ----------------------- | ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| Python ≥ 3.8           | `python3 --version`   | Required                                                                                                                                    |
| Git                     | `git --version`       | Required                                                                                                                                    |
| Apptainer / Singularity | `apptainer --version` | Required - installed by `setup_platform.sh`                                                                                              |
| 16 GB RAM               | -                      | Recommended                                                                                                                                 |
| 50 GB free disk         | -                      | Containers + outputs                                                                                                                        |
| FreeSurfer licence      | -                      | Required for `process-all` / `process-recon`. Free from [surfer.nmr.mgh.harvard.edu](https://surfer.nmr.mgh.harvard.edu/registration.html) |

**Windows only - install WSL2 first:**

```powershell
# Run in PowerShell as Administrator, then reboot
wsl --install
```

After reboot, open the **Ubuntu** app and run all steps below inside it.

---

## Install

```bash
# 1. Get the code
git clone https://github.com/stylianosc/leukoquant.git
cd leukoquant

# 2. Install container runtime (once per machine)
bash setup_platform.sh

# 3. Install the Python package (once per machine)
bash install.sh
```

Open a **new terminal**, then verify:

```bash
leukoquant --version
```

---

## FreeSurfer licence

Place your FreeSurfer `license.txt` inside `leukoquant/licenses/` (this folder exists in the repo and is never committed):

```
leukoquant/licenses/license.txt
```

A FreeSurfer licence is free from [surfer.nmr.mgh.harvard.edu](https://surfer.nmr.mgh.harvard.edu/registration.html).

---

## Quick start

```bash
leukoquant process-all \
  --subject    sub-01 \
  --t1         /data/sub-01_T1.nii.gz \
  --flair      /data/sub-01_FLAIR.nii.gz \
  --dwi        /data/sub-01_DWI.nii.gz \
  --bvecs      /data/sub-01.bvec \
  --bvals      /data/sub-01.bval \
  --output-dir /data/output \
  --scheduler  sge
```

---

## Commands

| Command                                | Description                             |
| -------------------------------------- | --------------------------------------- |
| `process-all`                        | Full end-to-end pipeline                |
| `process-recon`                      | FreeSurfer recon-all                    |
| `process-gif`                        | GIF brain parcellation                  |
| `process-bamos`                      | BaMoS WMH segmentation                  |
| `process-dti`                        | DTI fitting (FA, MD, AD, RD)            |
| `process-noddi`                      | NODDI fitting (NDI, ODI, ISOVF)         |
| `process-tracula`                    | TRACULA tractography                    |
| `process-metrics`                    | WMH-informed lesion/tract metrics       |
| `process-tract-qc`                   | Tractography QC                         |
| `process-zscore`                     | Z-score normalisation vs healthy cohort |
| `process-atlas-conversion`           | Atlas format conversion                 |

```bash
leukoquant <command> --help   # full options for any command
```

---

## Documentation

Full documentation is available at **[https://stylianosc.github.io/leukoquant/](https://stylianosc.github.io/leukoquant/)**.

| Section | Description |
|---------|-------------|
| [User Guide](https://stylianosc.github.io/leukoquant/user_guide/) | Installation, configuration, and workflow concepts |
| [Examples](https://stylianosc.github.io/leukoquant/examples_guide/) | Worked examples with real command invocations |
| [Output Structure](https://stylianosc.github.io/leukoquant/output_structure/) | Directory layout produced by each command |
| [API Reference](https://stylianosc.github.io/leukoquant/api/) | Full parameter reference for all 11 CLI commands |
| [Citations](https://stylianosc.github.io/leukoquant/citations/) | How to cite LeukoQuant and its dependencies |

---

## Citing LeukoQuant

If you use LeukoQuant in your research, please cite (also available as [CITATION.bib](CITATION.bib)):

```bibtex
@software{leukoquant2026,
  author    = {Charalampous, Stylianos and Barkhof, Frederik and Cardoso, M. Jorge and Sudre, Carole},
  title     = {{leukoquant}: Lesion-Informed White Matter Damage Metrics Toolkit for Cerebral Small Vessel Disease},
  year      = {2026},
  version   = {1.5.0},
  url       = {https://github.com/stylianosc/leukoquant},
  license   = {Apache-2.0}
}
```

LeukoQuant depends on several published tools. Please also cite them as appropriate: see [CITATIONS.md](CITATIONS.md) for the full reference list.

---

## Contributing

Bug reports and feature requests are welcome via [GitHub Issues](https://github.com/stylianosc/leukoquant/issues). See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup and contribution guidelines.

---

## Licence

Apache 2.0 - see [LICENSE](LICENSE) for details.
