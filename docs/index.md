# LeukoQuant

**Lesion-informed white matter damage metrics toolkit for cerebral small vessel disease.**

LeukoQuant provides a complete, reproducible pipeline from raw MRI scans to quantitative white matter metrics. It wraps Snakemake workflows inside a clean CLI, supports local and SGE cluster execution, and handles all container-based software dependencies automatically.

---

## Pipeline overview

```
Raw MRI (T1, FLAIR, DWI)
        │
        ├─ process-gif          → GIF brain segmentation & parcellation
        ├─ process-recon        → FreeSurfer recon-all
        ├─ process-bamos        → BaMoS white matter hyperintensity segmentation
        ├─ process-dti          → DTI fitting (FA, MD, …)
        ├─ process-noddi        → NODDI fitting (ODI, ICVF, …)
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

## Quick start

```bash
pip install leukoquant

# Single subject, full pipeline
leukoquant process-all \
  --subject sub-001 \
  --t1    ./data/sub-001/T1/scan.nii.gz \
  --flair ./data/sub-001/FLAIR/scan.nii.gz \
  --dwi   ./data/sub-001/DWI/data.nii.gz \
  --output-dir ./outputs
```

## Documentation

| Section | Description |
|---------|-------------|
| [User Guide](user_guide.md) | Installation, configuration, and workflow concepts |
| [Examples](examples_guide.md) | Worked examples with real command invocations |
| [Output Structure](output_structure.md) | Directory layout produced by each command |
| [API Reference](api/index.md) | Full parameter reference for all 11 CLI commands |
