# Process All Workflow - Output Structure

This document defines the standardized output folder structure for all processing modules in the `process_all_workflow`. Each module follows the pattern: `outputs/`, `intermediate/`, and `logs/` folders.

## Overview Structure

```
process_all_results/
├── {SUBJ_ID}/
│   ├── recon-all/
│   │   ├── outputs/
│   │   ├── intermediate/
│   │   └── logs/
│   ├── tracula/
│   │   ├── outputs/
│   │   ├── intermediate/
│   │   └── logs/
│   ├── bamos/
│   │   ├── outputs/
│   │   └── logs/
│   ├── gif/
│   │   ├── outputs/
│   │   ├── intermediate/
│   │   └── logs/
│   ├── dti/
│   │   ├── outputs/
│   │   ├── intermediate/
│   │   └── logs/
│   ├── noddi/
│   │   ├── outputs/
│   │   ├── intermediate/
│   │   └── logs/
│   ├── metrics/
│   │   ├── outputs/
│   │   └── logs/
│   └── tract_qc/
│       ├── outputs/
│       └── logs/
```

---

## Module Details

### recon-all/

FreeSurfer anatomical segmentation and surface reconstruction.

```
recon-all/
├── outputs/
│   ├── mri/
│   │   ├── aparc+aseg.mgz                   [OUTPUT - aparc atlas with subcortical]
│   │   ├── aparc.DKTatlas+aseg.mgz          [OUTPUT - DKT atlas with subcortical]
│   │   ├── aseg.mgz                         [OUTPUT - automatic segmentation]
│   │   ├── wmparc.mgz                       [OUTPUT - white matter parcellation]
│   │   ├── ThalamicNuclei.v13.T1.FSvoxelSpace.mgz  [OUTPUT - thalamic nuclei segmentation]
│   │   ├── ThalamicNuclei.v13.T1.volumes.txt       [OUTPUT - thalamic volumes]
│   │   └── [...all FreeSurfer mri outputs...]
│   ├── surf/
│   │   ├── lh.pial
│   │   ├── rh.pial
│   │   ├── lh.white
│   │   ├── rh.white
│   │   └── [...all cortical surface files...]
│   ├── label/
│   │   └── [...all cortical label files...]
│   ├── stats/
│   │   └── [...all cortical statistics...]
│   ├── scripts/
│   │   └── [...FreeSurfer processing records...]
│   └── touch/
│       └── [...FreeSurfer status files...]
└── logs/
    ├── recon_all_subj-XXX.log               [LOG]
    └── thalamic_seg_subj-XXX.log            [LOG]
```

---

### tracula/

Tractography and fiber tract segmentation from FreeSurfer.

```
tracula/
├── outputs/
│   ├── dpath/
│   │   ├── merged_avg16_syn_bbr.mgz         [OUTPUT - merged tract template]
│   │   ├── cc.bodyc_avg16_syn_bbr/
│   │   │   ├── path.pd.nii.gz               [OUTPUT - streamline count map]
│   │   │   ├── path.pd.trk                  [OUTPUT - tractography file]
│   │   │   └── [...endpoint and statistics files...]
│   │   └── [... 42 tract directories with similar structure ...]
│   └── dmri/
│       └── nodif_brain_mask.nii.gz          [OUTPUT - brain mask used by metrics/dti/noddi]
├── intermediate/
│   └── dmrirc                               [INTERMEDIATE - TRACULA config file]
└── logs/
    └── tracula_all.log                      [LOG]
```

---

### bamos/

Lesion segmentation and correction from multimodal imaging.

```
bamos/
├── outputs/
│   └── CorrectLesion_{SUBJECT_ID}.nii.gz    [OUTPUT - corrected lesion mask]
├── intermediate/
│   ├── Correct_WS3WT3WC1Lesion_{SUBJECT_ID}_corr.nii.gz  [INTERMEDIATE - intermediate segmentation]
│   ├── Connect_WS3WT3WC1Lesion_{SUBJECT_ID}_corr.nii.gz  [INTERMEDIATE - connectivity map]
│   ├── Connect_WS3WT3WC1MergedLesion_{SUBJECT_ID}_Test_corr.nii.gz  [INTERMEDIATE]
│   ├── PrimaryLesions_{SUBJECT_ID}.nii.gz  [INTERMEDIATE]
│   ├── SecondaryLesions_{SUBJECT_ID}.nii.gz  [INTERMEDIATE]
│   ├── Mask_{SUBJECT_ID}.nii.gz             [INTERMEDIATE - brain mask]
│   ├── T1_{SUBJECT_ID}.nii.gz               [INTERMEDIATE - working copy]
│   ├── FLAIR_{SUBJECT_ID}.nii.gz            [INTERMEDIATE - working copy]
│   ├── DataCorrected_T1FLAIR_BiASM_{SUBJECT_ID}_TA.nii.gz  [INTERMEDIATE]
│   ├── T1FLAIR_BiASM_{SUBJECT_ID}_TA.nii.gz  [INTERMEDIATE]
│   ├── SecondaryCorrected_WS3WT3WC1WMI2Ventr1Parc1JC1SP1ST1CL2CIV1_TOT1FLAIR_BiASM_{SUBJECT_ID}_TA.nii.gz  [INTERMEDIATE]
│   ├── Secondary_WS3WT3WC1WMI2Ventr1Parc1JC1SP1ST1CL2CIV1_TOT1FLAIR_BiASM_{SUBJECT_ID}_TA.nii.gz  [INTERMEDIATE]
│   ├── LesionCorrected_WS3WT3WC1WMI2Ventr1Parc1JC1SP1ST1CL2CIV1_TOT1FLAIR_BiASM_{SUBJECT_ID}_TA.nii.gz  [INTERMEDIATE]
│   ├── [...many more intermediate processing files...]
│   ├── ScriptBaMoS_{SUBJECT_ID}_TA.sh      [INTERMEDIATE - generated script]
│   ├── Aff_FLAIRtoT1.txt                    [INTERMEDIATE - registration params]
│   └── T1FLAIR_BiASM_{SUBJECT_ID}_TA.txt   [INTERMEDIATE - registration info]
└── logs/
    ├── bamos.log                            [LOG]
    └── bamos_correction.log                 [LOG]
```

---

### gif/

GIF (Group-wise Image registration for Functional analysis) brain segmentation and parcellation.

```
gif/
├── outputs/
│   ├── sub-01_T1_NeuroMorph_Parcellation.nii.gz       [OUTPUT - brain parcellation]
│   ├── sub-01_T1_NeuroMorph_Segmentation.nii.gz       [OUTPUT - tissue segmentation]
│   ├── sub-01_T1_NeuroMorph_Brain.nii.gz              [OUTPUT - skull-stripped brain]
│   └── sub-01_T1_NeuroMorph_TIV.nii.gz                [OUTPUT - total intracranial volume]
├── intermediate/
│   ├── sub-01_T1_NeuroMorph_BiasCorrected.nii.gz      [INTERMEDIATE - bias corrected T1]
│   ├── sub-01_T1_NeuroMorph_prior.nii.gz              [INTERMEDIATE - probability maps]
│   ├── sub-01_T1_NeuroMorph.xml                       [INTERMEDIATE - GIF metadata]
│   └── sub-01_T1_gw_affine.txt                        [INTERMEDIATE - registration parameters]
└── logs/
    └── gif.log                                      [LOG]
```

---

### dti/

Diffusion Tensor Imaging metrics extraction from diffusion-weighted imaging.

```
dti/
├── outputs/
│   ├── fa.nii.gz                            [OUTPUT - Fractional Anisotropy]
│   ├── md.nii.gz                            [OUTPUT - Mean Diffusivity]
│   ├── ad.nii.gz                            [OUTPUT - Axial Diffusivity]
│   ├── rd.nii.gz                            [OUTPUT - Radial Diffusivity]
│   ├── rgb.nii.gz                           [OUTPUT - RGB directional encoding]
│   ├── fa_inverse.nii.gz                    [OUTPUT - Inverse FA (1-FA)]
│   ├── l1.nii.gz                            [OUTPUT - Primary eigenvalue]
│   ├── l2.nii.gz                            [OUTPUT - Secondary eigenvalue]
│   └── l3.nii.gz                            [OUTPUT - Tertiary eigenvalue]
├── intermediate/
│   ├── masked_dwi.nii.gz                    [INTERMEDIATE - masked diffusion volume]
│   └── metadata.json                        [INTERMEDIATE - processing metadata]
└── logs/
    └── dti.log                              [LOG]
```

---

### noddi/

Neurite Orientation Dispersion and Density Imaging microstructural modeling.

```
noddi/
├── outputs/
│   ├── odi.nii.gz                           [OUTPUT - Orientation Dispersion Index]
│   ├── ndi.nii.gz                           [OUTPUT - Neurite Density Index (multi-shell only)]
│   └── fwf.nii.gz                           [OUTPUT - Free Water Fraction (multi-shell only)]
├── intermediate/
│   ├── binarized_mask.nii.gz                [INTERMEDIATE - binary brain mask]
│   ├── dir.nii.gz                           [INTERMEDIATE - direction map]
│   ├── metadata.json                        [INTERMEDIATE - processing metadata]
│   └── scheme.txt                           [INTERMEDIATE - acquisition scheme]
└── logs/
    └── noddi.log                            [LOG]
```

---

### metrics/

Metric extraction along tracts, lesions, and ROIs.

```
metrics/
├── outputs/
│   ├── whole_brain_metrics.csv              [OUTPUT - summary statistics across whole brain]
│   ├── skeleton_metrics.csv                 [OUTPUT - FA skeleton statistics]
│   ├── tract_level_metrics.csv              [OUTPUT - per-tract metric values]
│   ├── wmh_lesion_metrics.csv               [OUTPUT - WMH lesion-specific metrics]
│   └── wmh_tract_aggregated_lesion_metrics.csv  [OUTPUT - tract-level aggregated lesion metrics]
└── logs/
    └── metrics_log.txt                      [LOG]
```

---

### tract_qc/

Tract quality control reporting and validation.

```
tract_qc/
├── outputs/
│   └── qc_report.csv                        [OUTPUT - tractography QC metrics per tract]
└── logs/
    └── tract_qc_log.txt                     [LOG]
```

