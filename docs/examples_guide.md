# leukoquant CLI Usage Guide

All commands accept either a **bare subject ID** or a **path to a text file** (one ID per line) via `--subject`.
File-based inputs use `{subject}` as a placeholder in path arguments so each subject's files are resolved automatically.

---

## Subject input formats

```bash
# Single subject - bare ID
--subject sub-001

# Multiple subjects - text file (one ID per line)
--subject subjects.txt
```

When using a subjects file, path arguments support `{subject}` glob patterns:

```bash
--t1 "./SCANS/{subject}/T1/I*.nii.gz"
```

If `{subject}` is absent the path is treated as a literal (single-subject back-compat).

---

## process-gif

Run GIF brain segmentation and parcellation.
*Mandatory:* `--subject`, and at least one of `--t1` or `--flair`, `--output-dir`.
*Optional:* `--mask-file`.

When only `--flair` is given, GIF runs on the FLAIR as the primary input (useful when no T1 is available or when a FLAIR-trained atlas is used). When both `--t1` and `--flair` are supplied, the T1 is used and FLAIR is ignored (the GIF script does not support multi-channel input).

**Single subject (T1 only):**
```bash
leukoquant process-gif \
  --subject sub-001/sub-001-ses-01 \
  --t1 ./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz \
  --output-dir ./examples/outputs/new/gif_single \
  --scheduler sge
```

**Multiple subjects (T1 only):**
```bash
leukoquant process-gif \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --t1 "./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz" \
  --output-dir ./examples/outputs/new/gif_multi \
  --scheduler sge
```

**Single subject - FLAIR only:**

When running GIF on FLAIR-only input, GIF uses a FLAIR-specific database (must contain `db.xml`). Pass `--gif-flair-db` to specify its path explicitly. The standard T1 database is used automatically when `--t1` is present.

```bash
leukoquant process-gif \
  --subject sub-001/sub-001-ses-01 \
  --flair ./examples/sample_files/SCANS/{subject}/FLAIR/I*.nii.gz \
  --output-dir ./examples/outputs/new/gif_flair_single \
  --scheduler sge
```

**Multiple subjects - FLAIR only:**
```bash
leukoquant process-gif \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --flair "./examples/sample_files/SCANS/{subject}/FLAIR/I*.nii.gz" \
  --output-dir ./examples/outputs/new/gif_flair_multi \
  --scheduler sge
```

**Multiple subjects - T1 primary, FLAIR also provided (T1 is used, FLAIR ignored by GIF):**
```bash
leukoquant process-gif \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --t1    "./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz" \
  --flair "./examples/sample_files/SCANS/{subject}/FLAIR/I*.nii.gz" \
  --output-dir ./examples/outputs/new/gif_t1_flair_multi \
  --scheduler sge
```


**With an optional brain mask:**
```bash
leukoquant process-gif \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --t1 "./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz" \
  --mask-file "./examples/sample_files/SCANS/{subject}/nodif_brain_mask.nii.gz" \
  --output-dir ./examples/outputs/new/gif_brain_mask \
  --scheduler sge
```

---

## process-bamos

Run BaMoS white matter hyperintensity detection (requires GIF output).
*Mandatory:* `--subject`, `--flair`, `--t1`, `--output-dir`.
*Optional:* `--gif-results-dir` (runs GIF automatically if omitted).

**Single subject:**
```bash
leukoquant process-bamos \
  --subject sub-001/sub-001-ses-01 \
  --flair ./examples/sample_files/SCANS/{subject}/FLAIR/I*.nii.gz \
  --t1    ./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz \
  --output-dir ./examples/outputs/new/bamos_single \
  --scheduler sge 
```

**Multiple subjects:**
```bash
leukoquant process-bamos \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --flair "./examples/sample_files/SCANS/{subject}/FLAIR/I*.nii.gz" \
  --t1    "./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz" \
  --output-dir ./examples/outputs/new/bamos_multi \
  --scheduler sge
```

**Multiple subjects - reusing GIF outputs from a previous `process-gif` run:**
```bash
leukoquant process-bamos \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --flair "./examples/sample_files/SCANS/{subject}/FLAIR/I*.nii.gz" \
  --t1    "./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz" \
  --gif-results-dir "./examples/outputs/new/gif_multi/{subject}/gif/outputs" \
  --output-dir ./examples/outputs/new/bamos_multi_existing_gif \
  --scheduler sge
```

---

## process-recon

Run FreeSurfer `recon-all`.
*Mandatory:* `--subject`, `--t1`, `--output-dir`.

**Single subject:**
```bash
leukoquant process-recon \
  --subject sub-001/sub-001-ses-01 \
  --t1 ./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz \
  --output-dir ./examples/outputs/new/recon_all_single \
  --scheduler sge
```

**Multiple subjects:**
```bash
leukoquant process-recon \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --t1 "./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz" \
  --output-dir ./examples/outputs/new/recon_all_multi \
  --scheduler sge
```

---

## process-tracula

Run TRACULA tractography.
*Mandatory:* `--subject`, `--dwi`, `--output-dir`.
*Optional:* `--freesurfer-recon-dir`, `--t1`, `--bvecs`, `--bvals`.

If `--freesurfer-recon-dir` is omitted, or if recon-all outputs are missing for any subject,
recon-all is run automatically before TRACULA. In that case `--t1` is required.


**Without existing recon-all (runs recon-all automatically):**
```bash
leukoquant process-tracula \
  --subject sub-001/sub-001-ses-01 \
  --dwi ./examples/sample_files/SCANS/{subject}/DWI/I*.zip \
  --t1 ./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz \
  --output-dir ./examples/outputs/new/tracula_single \
  --scheduler sge
```

**Multiple subjects:**
```bash
leukoquant process-tracula \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --dwi    "./examples/sample_files/SCANS/{subject}/DWI/I*.zip" \
  --t1 ./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz \
  --output-dir ./examples/outputs/new/tracula_multi \
  --scheduler sge
```

**With existing recon-all outputs:**
```bash
leukoquant process-tracula \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --dwi ./examples/sample_files/SCANS/{subject}/DWI/I*.zip \
  --freesurfer-recon-dir "./examples/outputs/new/recon_all_multi/{subject}/recon-all/outputs" \
  --output-dir ./examples/outputs/new/tracula_existing_recon_multi \
  --scheduler sge
```

---

## process-dti

Run DTI model fitting.
*Mandatory:* `--subject`, `--dwi`, `--output-dir`.
*Optional:* `--mask`, `--bvecs`, `--bvals`. If `--mask` is omitted, DTI is fit on the whole image.

**Single subject:**
```bash
leukoquant process-dti \
  --subject sub-001/sub-001-ses-01 \
  --dwi  ./examples/sample_files/SCANS/{subject}/DWI/I*.zip \
  --output-dir ./examples/outputs/new/dti_single \
  --scheduler sge
```

**Multiple subjects:**
```bash
leukoquant process-dti \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --dwi   "./examples/sample_files/SCANS/{subject}/DWI/I*.zip" \
  --output-dir ./examples/outputs/new/dti_multi \
  --scheduler sge 
```

---

## process-noddi

Run NODDI model fitting.
*Mandatory:* `--subject`, `--dwi`, `--output-dir`.
*Optional:* `--mask`, `--bvecs`, `--bvals`. If `--mask` is omitted, NODDI is fit on the whole image.

**Single subject:**
```bash
leukoquant process-noddi \
  --subject sub-001/sub-001-ses-01 \
  --dwi  ./examples/sample_files/SCANS/{subject}/DWI/I*.zip \
  --output-dir ./examples/outputs/new/noddi_single \
  --scheduler sge
```

**Multiple subjects:**
```bash
leukoquant process-noddi \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --dwi   "./examples/sample_files/SCANS/{subject}/DWI/I*.zip" \
  --output-dir ./examples/outputs/new/noddi_multi \
  --scheduler sge
```

---

## process-metrics

**Compute metrics along tracts, lesions, and ROIs.**
*Mandatory:* `--subject`, `--tractography-path`, `--output-dir`.
*Optional:* `--t1-path`, `--dwi-path`, `--lesion-path`, `--maps`, `--tract-mode`.

The `--tract-mode` argument controls the analysis approach:
- `tractography`: use tractography only
- `tractography-atlas`: use tractography, fall back to atlas if QC fails (default)
- `atlas`: atlas-based only

**Single subject (default tract-mode):**
```bash
leukoquant process-metrics \
  --subject sub-001/sub-001-ses-01 \
  --tractography-path ./results/tracula/{subject}:dpath/*/path.pd.nii.gz:dwi \
  --output-dir ./results/metrics
```

**Multiple subjects with lesion and DTI metrics:**
```bash
leukoquant process-metrics \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --tractography-path ./results/tracula:dpath/*/path.pd.nii.gz:dwi \
  --t1-path  "./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz:" \
  --dwi-path "./examples/sample_files/SCANS/{subject}/DWI/I*.zip:" \
  --lesion-path "wmh=./results/bamos:bamos/CorrectLesion_*.nii.gz:t1" \
  --maps dti_fa=./results/dti:fa.nii.gz:dwi \
  --maps dti_md=./results/dti:md.nii.gz:dwi \
  --tract-mode tractography-atlas \
  --output-dir ./results/metrics \
  --scheduler sge --cores 4
```

**Atlas-based only:**
```bash
leukoquant process-metrics \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --tractography-path ./results/tracula:dpath/*/path.pd.nii.gz:dwi \
  --tract-mode atlas \
  --output-dir ./results/metrics
```

---

## process-tract-qc

**Run tractography quality control.**
*Mandatory:* `--subject`, `--tractography-path`, `--output-dir`.

**Single subject:**
```bash
leukoquant process-tract-qc \
  --subject sub-001/sub-001-ses-01 \
  --tractography-path ./results/tracula/{subject}:dpath/*/path.pd.nii.gz:dwi \
  --output-dir ./results/tract_qc
```

**Multiple subjects:**
```bash
leukoquant process-tract-qc \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --tractography-path ./results/tracula:dpath/*/path.pd.nii.gz:dwi \
  --output-dir ./results/tract_qc \
  --scheduler sge --cores 4
```

---

## process-zscore

**Run a standalone Z-score analysis.**
*Mandatory:* `--healthy-list`, `--target-list`, `--t1-path`, `--metric`, `--output-dir`.
*Optional:* `--demographics-csv`, `--covariates`, `--polynomial-terms`.

```bash
leukoquant process-zscore \
  --healthy-list ./examples/sample_files/SCANS/healthy_subjects.txt \
  --target-list  ./examples/sample_files/SCANS/target_subjects.txt \
  --t1-path "./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz" \
  --demographics-csv ./examples/sample_files/SCANS/demographics.csv \
  --covariates age,sex \
  --polynomial-terms age:2 \
  --metric  "dti_fa=./results/dti/{subject}:fa.nii.gz:dwi" \
  --metric  "dti_md=./results/dti/{subject}:md.nii.gz:dwi" \
  --output-dir ./results/zscores \
  --scheduler sge --cores 4
```

---

## process-atlas-conversion

Convert a brain parcellation file from GIF label space to FreeSurfer aparc+aseg label space using a CSV mapping file.
Output is written to `{output_dir}/{subject}/outputs/converted_atlas.mgz`.

*Mandatory:* `--subject`, `--input-parcellation`, `--mapping-file`.
*Optional:* `--output-dir` (defaults to the parent directory of `--input-parcellation`), `--no-validate`, `--scheduler`, `--cores`.

**Single subject - convert GIF parcellation to FreeSurfer format:**
```bash
leukoquant process-atlas-conversion \
  --subject subj-001 \
  --input-parcellation ./results/subj-001/gif/outputs/subj-001_NeuroMorph_Parcellation.nii.gz \
  --mapping-file ./leukoquant/mappings/gif_to_freesurfer.csv \
  --output-dir ./results \
  --scheduler sge
```
The converted atlas is written to `./results/subj-001/atlas_conversion/outputs/converted_atlas.mgz`.

**Single subject - skip validation:**
```bash
leukoquant process-atlas-conversion \
  --subject subj-001 \
  --input-parcellation ./results/subj-001/gif/outputs/subj-001_NeuroMorph_Parcellation.nii.gz \
  --mapping-file ./leukoquant/mappings/gif_to_freesurfer.csv \
  --output-dir ./results \
  --no-validate \
  --scheduler local
```

---

## process-all

Run the complete end-to-end pipeline (recon-all → GIF → BaMoS → TRACULA → DTI → NODDI → metrics).
Multiple subjects are launched as concurrent Snakemake processes.

**Single subject (using real sample files):**
*Mandatory:* `--subject`, `--t1`, `--flair`, `--dwi`, `--output-dir`.

```bash
leukoquant process-all \
  --subject sub-001/sub-001-ses-01 \
  --t1    "./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz" \
  --flair "./examples/sample_files/SCANS/{subject}/FLAIR/I*.nii.gz" \
  --dwi   "./examples/sample_files/SCANS/{subject}/DWI/I*.zip" \
  --output-dir ./examples/outputs/new/process_all_single \
  --scheduler sge
```

**Multiple subjects:**
```bash
leukoquant process-all \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --t1    "./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz" \
  --flair "./examples/sample_files/SCANS/{subject}/FLAIR/I*.nii.gz" \
  --dwi   "./examples/sample_files/SCANS/{subject}/DWI/I*.zip" \
  --output-dir ./examples/outputs/new/process_all_multi \
  --scheduler sge
```

**Single subject - using GIF parcellation instead of FreeSurfer for TRACULA:**
GIF runs on the T1 first; its parcellation output is then converted to FreeSurfer label space and used by TRACULA in place of the standard `aparc+aseg`. Both single- and multi-subject runs are supported.
```bash
leukoquant process-all \
  --subject sub-001/sub-001-ses-01 \
  --t1    "./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz" \
  --flair "./examples/sample_files/SCANS/{subject}/FLAIR/I*.nii.gz" \
  --dwi   "./examples/sample_files/SCANS/{subject}/DWI/I*.zip" \
  --parcellation gif \
  --output-dir ./examples/outputs/new/process_all_gif_parcellation \
  --scheduler sge
```

**Multiple subjects - using GIF parcellation:**
```bash
leukoquant process-all \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --t1    "./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz" \
  --flair "./examples/sample_files/SCANS/{subject}/FLAIR/I*.nii.gz" \
  --dwi   "./examples/sample_files/SCANS/{subject}/DWI/I*.zip" \
  --parcellation gif \
  --output-dir ./examples/outputs/new/process_all_gif_parcellation_multi \
  --scheduler sge
```

**End-to-end with Z-score against a healthy cohort:**
Runs the full pipeline for all subjects (targets + healthy). Z-scores are computed only for target subjects against the healthy cohort. Metrics extraction for each target subject runs after its z-score maps are ready.
*Optional (for Z-score):* `--healthy-subjects`, `--demographics-csv`, `--covariates`, `--poly-terms`.

```bash
leukoquant process-all \
  --subject ./examples/sample_files/SCANS/target_subjects.txt \
  --healthy-subjects ./examples/sample_files/SCANS/healthy_subjects.txt \
  --t1    "./examples/sample_files/SCANS/{subject}/T1/I*.nii.gz" \
  --flair "./examples/sample_files/SCANS/{subject}/FLAIR/I*.nii.gz" \
  --dwi   "./examples/sample_files/SCANS/{subject}/DWI/I*.zip" \
  --demographics-csv ./examples/sample_files/SCANS/demographics.csv \
  --covariates age,sex \
  --poly-terms age:2 \
  --output-dir ./examples/outputs/new/process_all_z_score \
  --scheduler sge
```
