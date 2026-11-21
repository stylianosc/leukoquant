# Changelog

All notable changes to LeukoQuant are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [2.3.1] - 2026-11-21

Release 2.3.1 completes the skull-strip pin of 2.3.0. `MKL_CBWR=COMPATIBLE` alone did not give the same T1 mask on the AVX-only abner nodes. The metrics step now also exports `ONEDNN_MAX_CPU_ISA=SSE41` before any `mri_synthstrip` starts, and with both settings the b0 mask, the T1 mask and the `reg_aladin` affine were identical on every node family tested.

### Fixed
- `pin_mkl_kernel` (`bash_utils.sh`) now exports `ONEDNN_MAX_CPU_ISA=SSE41` as well as `MKL_CBWR=COMPATIBLE`, and logs one `[mkl]` line for each; a value already set by the caller wins for that variable. Evidence, `mri_synthstrip` masks and `reg_aladin` affines for 3 ADNI-3 sessions, 2 repeats each, on abner (AVX), a Haswell node (AVX2), a Xeon Gold node and an AMD EPYC node (AVX-512) and arbuckle (SSE4.2): 30 of 30 runs identical per session. Neither setting alone does it (`MKL_CBWR` differs on abner, `ONEDNN_MAX_CPU_ISA` differs between node families); one thread, `MKL_ENABLE_INSTRUCTIONS` and the OpenBLAS pin did not change the mask. No measurable change in run time.

### Known issues
- Masks and tables made without both settings, including those of 2.3.0 on AVX-only nodes, can differ from pinned ones, so they should not be mixed with 2.3.1 tables in one analysis.
- Equality of whole metrics tables across node types has still not been compared end to end, so byte-reproducibility across nodes is not claimed. The DWI-source rule described under 2.3.0 is unchanged.

---

## [2.3.0] - 2026-11-21

Release 2.3.0 pins the numerical kernels of the metrics step (OpenBLAS for numpy and MKL for the skull-strip network), which reduces the CPU as a source of differences between nodes in the parts we tested (one known exception, see Known issues), makes the lesion shape eigenvalues reproducible from run to run, and records where each number came from. End-to-end equality of the metrics tables across node types has not been verified (see Known issues), and tables generated before the pins differ slightly from pinned ones, so tables from the two should not be mixed in one analysis. The metrics tables of the accompanying data release are described with the paper, not here.

### Added
- `process-all` now prints a `[dwi]` line naming the TRACULA run that supplies the DWI to DTI, NODDI and the metrics step, and warns when it is not the FreeSurfer one (`utils/dwi_source.py`, the rule itself is unchanged).
- Per-job usage records with true peaks (`utils/usage_record.py`, `utils/usage_report.py`): every `run_with_peak_mem` stage of the TRACULA rule appends one line (exact peak RSS, sampled VmPeak/VmHWM/sum of address space, peak scratch, minimum free scratch on the node, elapsed, exit status, killed, host, job/task id, request) to the subject's `logs/usage_records.tsv`; `usage_report` merges them into `pipeline_status/usage/usage_table.tsv`, flags runs above 80 % of the request and refits the model from the true peaks. Sampling (default every 30 s) can miss short spikes; the scratch factor 2.0 applies only while the model's `scratch_basis` is not `peak`.
- TRACULA now also keeps `dmri/lowb_brain.nii.gz` (the skull-stripped b0) in the subject's `outputs/dmri/`, next to `dwi.nii.gz` and `nodif_brain_mask.nii.gz`. It is copied from the shell block only and is not a declared rule output, so subjects processed earlier are not rescheduled and simply lack the file.
- The memory and the node-local scratch (SGE `tscratch`) requested by a TRACULA job are now continuous functions of the DWI size, replacing the flat 8 GB and 3 GB. Size is megavoxel-volumes (voxels x volumes / 1e6), read from the zip (uncompressed bytes / 2) or the NIfTI header. Each quantity is a least-squares line plus a margin equal to the largest residual of the fit (so the request covers every measured run, in-sample), rounded up to whole GB and clamped only at the ends: memory 8-16 GB (the larger of the prep and path predictions), scratch 3-40 GB. The scratch fit is on the scratch in use at the END of a run, while the peak during the run is higher and unknown, so the fitted scratch is multiplied by a PROVISIONAL peak factor of 2.0 (a stated guess, `TRACULA_SCRATCH_PEAK_FACTOR`); refit the scratch from the per-stage `[node]` df lines of the TRACULA log once runs with them exist and then drop the factor. Coefficients are in `leukoquant/utils/tracula_model/tracula_resource_model.json`, refitted from `tracula_model/tracula_job_fit_input.tsv` (6560 session x parcellation runs of ADNI-3, EPAD, OASIS3 and DPUK, ids removed) by `python -m leukoquant.utils.tracula_resource_model [table.tsv]`. Fit: scratch n=6274, R2 0.79, y = 0.497 + 0.01106 x Mvv + 1.766 GB; prep RSS n=151, R2 0.55, y = 1.450 + 0.00632 x Mvv + 1.378 GB; path RSS n=151, R2 0.49, y = 0.883 + 0.01575 x Mvv + 5.295 GB. Resulting requests (memory / scratch): 116x116x80x55 gets 8 GB / 6 GB, 256x256x80x36 gets 10 GB / 9 GB, 256x256x80x54 gets 11 GB / 11 GB. Limits of the data: no peak virtual memory (maxvmem) exists for finished jobs while the cap acts on address space, only about 150 runs have memory lines and almost all are small images, and one 10.9 GB outlier (`trac-all -path`, 256x256x80x54) sets the path margin. `LEUKOQUANT_TRACULA_MEM_GB` now accepts 4-24 GB (the cap) and `LEUKOQUANT_TRACULA_SCRATCH_GB` 3-40 GB; both still override.
- Every stage run through `run_with_peak_mem` (TRACULA prep, bedpostx, path) now also logs a `[vm]` line with the largest peak address space (`VmPeak`) and peak resident memory (`VmHWM`) of any one process of the stage, sampled from `/proc` every `LEUKOQUANT_VM_SAMPLE_S` seconds (default 60), next to the address-space limit in force. The cap acts on address space, which GNU time does not measure. A process that starts and ends between two samples is missed. Data for the next refit of the memory model.
- The TRACULA log now carries `[node]` blocks (host name, `df -h` of the scratch and tmp locations, size of the job's scratch folder) at the start, after prep, before and after bedpostx, after path and at job end, so a kill for lack of node scratch can be told from a memory kill.
- When the bedpostx stage of TRACULA fails, its logs (`dmri.bedpostX/logs`, in particular `bedpostx_post.e*`) are copied from node scratch to `tracula-<parcellation>/logs/bedpostx_logs` before the scratch is deleted, so the reason of the failure is no longer lost.
- The metrics log now records provenance: a `[host]` line (host name, CPU model, widest SIMD extension) at the start, and before the Python step one `[sha] <label> <file> <md5>` line per intermediate: the DWI-to-T1 registration folder (affines, skull-stripped b0), the atlas registration folder (affine, control points, deformation) and every T1-space input. Two runs of the same job can then be diffed to the first file that differs. Added because some EPAD batch runs gave tract volumes 0.3-0.7% apart with identical code, inputs and logs (2026-11-18). Helpers `log_host_cpu` and `log_checksums` in `bash_utils.sh`.

### Fixed
- The skull-strip masks of the metrics step depended on the node. `mri_synthstrip` (a PyTorch network on the CPU) gave a different T1 mask on an AVX2 node than on AVX-512 and AMD EPYC nodes for the same input, which moved the DWI-to-T1 affine and with it every DTI and NODDI value sampled in the lesions (ADNI-3). MKL chooses its kernels from the CPU at run time. `metrics_calc.sh` now calls `pin_mkl_kernel` (`bash_utils.sh`), which exports `MKL_CBWR=COMPATIBLE` before any `mri_synthstrip` starts, and logs a `[mkl]` line; a value already set by the caller wins. With it, an AVX2 node, an AVX-512 Xeon Gold node and an AMD EPYC node gave identical b0 and T1 masks and identical `reg_aladin` affines for three sessions, two repeats each, with no change in run time. The AVX-only node (abner) is an exception: for one session (subj-130-s-6611) it still gave a different T1 mask, so the setting does not make the mask independent of every node. A stronger pin is being tested. Setting `ONEDNN_MAX_CPU_ISA` alone, one thread or the OpenBLAS pin did not change the mask on the AVX2 node. Tables produced before this change on an AVX2 node can differ from the pinned result.
- The peak-address-space sampler (`start_vm_sampler`) no longer prints `syntax error in expression` into the job log. From the second sample on, the stored line was read into four variables, so the last one swallowed the job-wide fields and the next line written carried extra fields. A process name with spaces or an empty name (read as `Name:` in `/proc/<pid>/status`) would have shifted the fields as well, so names are now always one word. Only the diagnostics were affected; the stage itself ran.
- A run started from another checkout than the production one failed at the start of every cluster task with `cannot import name 'tracula_scratch_mb' from 'leukoquant.utils.mem_override'` (EPAD and OASIS3 pinned metrics runs, 4,044 tasks, 2026-11-19). The tasks run without `PYTHONPATH` and the environment's editable install points at production, so the first `leukoquant` import decides which code is used. `process_all_workflow.smk` imported `leukoquant.utils.z_score_utils` (only when z-score covariates are on) before any workflow had put the configured `leukoquant_parent_dir` on `sys.path`, which loaded and cached the production package; the later `sys.path.insert` in the sub-workflows could not undo that. `process_all_workflow.smk` now inserts `leukoquant_parent_dir` first, before any import, and `tract_qc_workflow.smk` imports `leukoquant` only after its insert. `test_workflow_import_path.py` checks the order in every workflow file.
- A killed stage is no longer called out of memory unless the numbers say so. `run_with_peak_mem` printed `[oom] <stage>: out of memory` for any exit 137 or any "Killed" line in the stage's stderr. In the 021-s-7055, 021-s-7062 and 021-s-7092 FreeSurfer TRACULA runs (ADNI-3, 2026-11) that "Killed" came from `bedpostx_mgh`'s own monitor, which runs `kill -9` on bedpostx as soon as any of its error files is non-empty. It happened at a 2.1 GB peak under a 24 GB limit, and a standalone run of the bedpostx post stage at that size peaked at 2.4 GB of address space, so the label pointed the investigation at memory for nothing. `_finish_stage_report` now prints `[oom]` only for an allocation-failure message in the stage's own stderr (`std::bad_alloc`, `Cannot allocate memory`, `MemoryError`, ...), or for a kill whose peak (the larger of the exact RSS peak and the sampled `VmPeak`) reached 90% of the address-space limit (`_OOM_NEAR_CAP_PCT`). Any other kill prints `[killed] <stage>: killed, cause unknown (...)` with the peaks and the limit. The usage record's `killed` column is 1 for both. `report_step_failure` no longer counts a bare "Killed" line as out of memory; its exit-137 rule is unchanged, since it has no peak to compare with. Parsers that only want real memory failures (`repair_census.py` matches `[oom] metrics:`) now see fewer false positives.
- A SIGKILLed stage now reports status 137 on every host. GNU time 1.7 (the cluster hosts) returns the bare signal number (9) for a killed command and writes nothing about the signal, so it looked like a plain `exit 9`. Newer GNU time (the containers) returns 137. `run_with_peak_mem` now also records GNU time's `%x` (the command's own exit status, 0 after a signal) and turns a non-zero return with `%x` 0 into 128 + signal. A real `exit 9` keeps status 9. This also made `test_loud_failures.py` pass on the cluster hosts, where one test had always failed.
- The metrics step gave different tract masks on different nodes. numpy's OpenBLAS picks its matrix kernel from the CPU at run time, and a kernel with fused multiply-add rounds differently from one without, so the last bits of the transformed streamline coordinates differed and points near a voxel boundary landed in another voxel (abner nodes, AVX only, against AVX2 and AVX-512 nodes; EPAD 2026-11-18). Every registration output and every DTI/NODDI map was identical across nodes; the first differing files were the tractography density maps of `transform_trk.py`. `metrics_calc.sh` now calls `pin_blas_kernel` (`bash_utils.sh`) before any Python starts, which exports `OPENBLAS_CORETYPE=Core2`, a kernel that needs only SSSE3 and so runs on every node, including those without AVX, and logs a `[blas]` line; a value already set by the caller wins. With the pin, density maps of four tracts were byte-identical on a node without AVX, an AVX node and two AVX2/FMA nodes. The same session run twice on one node, once with `Core2` and once with `SandyBridge`, gave identical intermediate images and maps and five identical tables out of six. Only five lesion eigenvalue columns of `wmh_lesion_metrics.csv` differed, in the last digit (at most about 2e-15 relative), so the choice of kernel matters at that level for those columns. `Core2` tables from a node without AVX and from an AVX node were byte-identical. Tables computed before this change on AVX2 and AVX-512 nodes differ slightly from the pinned result (a few voxels per tract) and have to be regenerated to be comparable.
- Lesion shape eigenvalues (`eigenvalue_1..3` and the three `ratio_*_eigenvalue` columns of the lesion tables) are now identical bit for bit from run to run. The floating-point PCA summed the coordinates in an order- and threading-dependent way, and for a flat lesion the smallest eigenvalue, mathematically zero, came out as rounding noise (about 1e-15 to 1e-33 mm^2) that the ratios divided by: two reruns of the metrics step with the same code and inputs gave ratios up to 19% apart (EPAD, 2026-11-18). The second moments are now summed exactly in integer arithmetic, voxel sizes are applied to the exact covariance, and an eigenvalue at or below 1e-12 of the largest is reported as exactly 0, so a ratio with that denominator is NaN (as for any zero denominator). Values for non-degenerate lesions are unchanged up to rounding; flat and linear lesions now get eigenvalue 0 and NaN ratios where they used to get large, run-dependent numbers.

### Known issues
- Equality of the metrics tables across node types has not been verified end to end. With `MKL_CBWR=COMPATIBLE` the T1 skull-strip mask was identical on AVX2, AVX-512 and AMD nodes, but not on the AVX-only abner nodes for sensitive sessions (subj-130-s-6611 gave a different mask there), so the pin is not a full node-independence guarantee and a stronger one is under test. The b0 masks and the b0-to-T1 registration were identical across the nodes tested for the sessions tested, and repeated runs on one node gave identical outputs, but whole tables were not compared across node types after the pins, so byte-reproducibility across nodes is not claimed.
- DTI, NODDI and the metrics step read one DWI per session, from the TRACULA run of the first parcellation in the run's list, FreeSurfer if it is listed. A run that lists only another parcellation (for example `--parcellation gif`) reads that parcellation's TRACULA `dwi.nii.gz`, while the DTI/NODDI maps on disk were fitted on the FreeSurfer one. The two files hold the same data up to rounding (about 6e-6 of the b0 mean, bit-identical in no volume for most sessions), and the metrics step amplifies that: in 240 ADNI-3 sessions run both ways, most DTI/NODDI-derived columns, tract-level and skeleton metrics differed by a few tenths of a percent at the median, with some cells much larger. Run every parcellation of a cohort together (`--parcellation freesurfer,gif`) and the tables of one cohort share one DWI. The new `[dwi]` line shows which DWI a run used.

---

## [2.2.0] - 2026-11-07

### Changed
- Z-score maps are now written in DWI space by default (`output_space` of the z-score stage in `process-all`, previously T1). The metrics stage tags them `dwi` and projects them to the T1 grid with the same DWI-to-T1 transform and linear interpolation as the raw dti / noddi maps, so every map is resampled in one place. The standalone `process-zscore` command follows: `--metric-space` defaults to `dwi` and `--output-space` defaults to the metric space, so nothing is registered unless requested (the maps are DWI maps and stay in their native space); `--dwi-pattern` is therefore required unless the spaces are set to `t1` explicitly. Z-score maps from earlier releases (on the T1 grid) keep working in the metrics stage: see the next entry.

### Added
- `image_utils.py --op same_grid --ref REF` (exit status 0 when the image shares REF's voxel grid by shape and affine) and `same_voxel_grid()`. `metrics_calc.sh` uses it so that an input tagged `dwi` that already sits on the T1 grid is copied instead of being sent through the DWI-to-T1 transform.

---

## [2.1.4] - 2026-11-07

### Fixed
- The z-score stage now reorients its target T1 to RAS+ (after the first frame) as 2.1.3 did for the metrics stage: z-score maps are written on the target T1's grid, so for a T1 stored in another orientation they stayed on a grid the lesion map and the metrics-stage maps no longer shared.

### Changed
- `metrics_calc.py` now checks that the lesion maps, maps, tract atlas, tract files and skeleton sit on the reference T1's voxel grid by shape AND affine (`_assert_same_voxel_grid`), before any loading. The previous shape-only check could not see two grids that differ by an axis flip (a T1 stored LAS against a lesion map in RAS has the same shape), so such inputs were paired voxel by voxel with the left and right sides swapped. The run now stops with the offending inputs listed.

---

## [2.1.3] - 2026-11-07

### Fixed
- The metrics stage now reorients the reference T1 to RAS+ (after taking its first frame), as BaMoS and GIF already do. A raw T1 stored in another orientation (the EPAD sessions exported as a 4D P,S,R image) was left on its own voxel grid, so the lesion map, parcellation and DWI-space maps no longer matched it and `metrics_calc.py` stopped on its shape check. `image_utils.py` gains `--op canonical` (reorient only, voxel values untouched, so T1s that are already RAS+ are copied byte for byte); `reorient_to_ras` takes `clip_negatives` (default unchanged).

---

## [2.1.2] - 2026-11-06

### Added
- `LEUKOQUANT_TRACULA_MEM_GB` (whole GB, 4-64) sets the memory request of the TRACULA jobs for one run; unset or empty keeps the
  8 GB default, an invalid value is ignored with a printed warning. ADNI-3 sessions peaked at 4.7-8.1 GB and several were
  killed at the 8 GB limit (`utils/mem_override.py`, tests in `tests/test_mem_override.py`).

### Fixed
- BaMoS lesion maps no longer carry the input scan's NIfTI `scl_slope`/`scl_inter`. BaMoS never applies that slope to the voxels but copied it into every image it wrote, so a lesion map holding weights in [0, 1] read back multiplied by the slope of the scan it came from (EPAD `CorrectLesion.nii.gz` came out as 2.94-5.87 instead of 0.5-1). `correction_lesions_111125.py` made it permanent: `get_fdata()` applied the slope and the result was saved with a fresh header. The BaMoS header copy (`nifti_copy_nim_info`, `CreateNiiFromArray`) now resets the scaling, and the correction script reads the stored values (`dataobj.get_unscaled()`). Binary lesion metrics (volumes, counts, overlaps) are unchanged; the lesion-weighted `*lesion_pd*` metrics of an affected subject were multiplied by its slope and need `CorrectLesion.nii.gz` regenerated and the metrics rerun. `metrics_calc.py` now warns when a lesion map exceeds 1. The shipped `Seg_Analysis` and `Seg_BiASM` binaries are rebuilt from the patched source (`build_bamos.sh`, gcc 9.2.0).

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
