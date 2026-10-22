# LeukoQuant metrics: revision plan (version 2.1.0)

Status: **plan, nothing implemented yet.** Written 2026-10-17 after the audit of `leukoquant/utils/metrics_calc.py` and
`lesion_metrics.py`.

## Status (2026-10-17, evening)

| Phase | State |
|---|---|
| 0 Design | decisions recorded (section 7); open: ring distances re-confirmation |
| 1 Output isolation | **done and dry-run verified** (`metrics_output_dir`) |
| 2 Implementation | **done** on branch `metrics-v2.1` (local commits, not pushed). Not done: streamline-level disconnection (atlas streamlines are in diffusion space, tractography ones in T1 space; replaced by the density-fraction proxy) |
| 3 Tests | **done**: 244 tests in neuro-analyser (phantoms, brute-force AUC, end-to-end run of `main()`), 62 in LeukoGraph |
| 4 Pilot | waiting: needs the branch on the cluster (push), in a separate clone so running jobs do not change code |
| 5 Full recompute | on hold until the LeukoGraph run has finished (scheduling rule) |
| 6 Downstream | **done locally** (LeukoGraph commit 342ebbb, not pushed and not pulled on the cluster while the run is going) |
| 7 Second LeukoGraph run | after the first results |

Found on the way: the whole-brain skeleton statistics (PSMD, mean skeleton FA / MD) were statistics of map x fractional weight
(linear resampling of the thin skeleton); fixed in 2.1.0. Details in `CHANGELOG.md`.

## 0. Principles

0. **Scheduling rule (user, 2026-10-17): the LeukoGraph run has priority on the cluster.** Everything up to the pilot is
   local work and tiny test jobs. The *full* metrics recompute (Phase 5) is submitted only after the LeukoGraph run has finished,
   so it never queues in front of graph jobs.

1. **The running LeukoGraph experiment is not touched.** It finishes on the current (version 2.0.0) metrics and gives the first
   results. The revised metrics are developed and computed in the background, in separate folders, and used for a second
   LeukoGraph run only afterwards (section 7).
2. **Nothing is written over existing outputs.** The new tables go to their own location; the 2.0.0 tables stay as they are.
3. **No missing values.** Where a quantity is undefined it gets a neutral value (statistics 0, difference 0, AUC 0.5), and every
   region carries a voxel count so that a real 0 can be told apart from "nothing there".
4. **Every change has a test** that a reader can follow (a small synthetic example with a known answer).
5. **No push, no release** (GitHub, Zenodo) without your explicit go-ahead.

## 1. Why (evidence from the audit, ADNI-3 unless stated)

| Finding | Evidence |
|---|---|
| Lesion-vs-outside "change" ratio is unbounded on signed Z-score maps | 11% / 15% / 18% of values outside [-1, 1] (ADNI-3 / OASIS-3 / EPAD); maximum about 120,000 |
| `merged_after_dilation` carries no information | true for 23,995 of 23,995 lesions |
| `coalescence_score` can only be 0 or 1 | 22,513 zeros, 1,482 ones; size thresholds compared with a fraction, not a percentage |
| `fractal_dimension` undefined (0) for 14% of lesions; unreliable for most | median lesion 14 voxels; 73% under 30 voxels |
| Region statistics of empty regions are 0 and indistinguishable from real zeros | code, `_summary_nonzero` |
| Tract statistics use map x tract probability (not a weighted mean) | code, `_compute_tract_level_metrics_fast` |
| "Outside" reference includes the penumbra and other lesions | code, `_fast_summary_change` |
| 6- versus 26-connected lesion components in different tables | code, `_remove_small_or_flat_components` vs `lesion_metrics.py` |

## 2. The changes (agreed unless marked "to confirm")

**New**
- AUC per map (11 columns): chance that a random lesion voxel exceeds a random voxel of the reference tissue.
- Difference (inside minus reference) for each of the 10 statistics, per map (110 columns).
- `p5` and `p95` for every region; `n_voxels` for every region.
- NAWM region (normal-appearing white matter: the tract minus the dilated lesion) with all statistics; change measured against it.
- Penumbra rings 3-6 mm and 6-9 mm (reduced statistics: mean, median, std, p5, p95).
- Lesion thickness (largest distance from the lesion edge inward, mm); number of lesions merged after dilation.
- Tract-to-lesion distance: mean distance of tract voxels to the nearest lesion; fraction of the tract within 5 mm.
- Top-2 and top-3 tract matches with overlap fractions.
- Whole-brain lesion count and mean lesion size.
- To confirm: streamline disconnection (fraction of a tract's streamlines through lesions); periventricular / deep lesion volumes.

**Replaced**
- `merged_after_dilation`: voxel count -> number of distinct lesions in the dilated blob.
- Coalescence score: fraction -> fraction x 100 (same bands; supervisor's formula kept).
- Fractal dimension: scales limited to the lesion's extent.
- Surface area / sphericity / compactness: exposed voxel faces -> mesh surface.
- Region statistics weighted by tract probability: map x density -> weighted mean, std, percentiles.
- Penumbra: 3-voxel dilation -> the 0-3 mm ring.
- Component connectivity: 6 and 26 -> 26 everywhere.
- Centre-of-mass distance -> the two distances above. Tract match: top 1 -> top 3.

**Deleted**
- The `_change` and `_change_absolute` ratio features (220 columns per cohort), replaced by difference and AUC.
- `min` and `max` (replaced by `p5`, `p95`); centre-of-mass distance.
- To confirm: `volume_tiv_normalised_mm3` (exactly load % / 100).

Kept: IQR, 25th, 75th, sum, skew, kurtosis, peak width, the Z-score model, the skeleton (PSMD-type) statistics.

## 3. Phases

| Phase | What | Output / gate |
|---|---|---|
| **0. Freeze the design** | You decide the "to confirm" items, ring distances, version number. Check whether the paper uses any Z-score change feature. | Short decision list. **Gate: your approval.** |
| **1. Output isolation** | New config key `metrics_output_dir` (default: the usual `output_dir`). Only the metrics step's own folders (outputs, logs, work directory) move to it, in the same layout as today: `<metrics_output_dir>/<subject>/metrics-gif/outputs/metrics/*.csv`. Everything the step reads (T1, lesions, tracts, QC report, DTI / NODDI / Z-score maps) stays where it is. Edits: `metrics_workflow.smk` (about 10 path strings), `process_all_workflow.smk` (pass the key through; the two `makedirs` lines), and data collection is pointed at the new root. The new root holds only the new tables, so nothing can be mixed up with the old ones. | Workflow change + test on one subject |
| **2. Implementation** (branch `metrics-v2.1`, no push) | In order of risk: statistics set and `n_voxels`; difference + AUC; weighted statistics; NAWM and rings; per-lesion fixes; distances; whole-brain additions; tract matches; the optional ones. | One commit per change, tests green after each |
| **3. Tests** | Synthetic phantoms with known answers (sphere, two-tissue mixtures for AUC and difference); regression test that every unchanged metric is identical to 2.0.0; property checks (AUC in [0, 1], no NaN or Inf); timing test. | Test suite, benchmark |
| **4. Pilot** | 10 ADNI-3 subjects on the cluster, isolated output. Check ranges, counts, runtime, and the new column count. | Pilot report. **Gate: your approval.** |
| **5. Full recompute** | All four datasets on SGE, concurrency capped so the running experiment is not starved. | New tables in the isolated location |
| **6. Downstream adaptation** | Classify the new columns in `TractAutoencoder/model/metric_taxonomy.py`, update `shared/feature_rules.py` (empty-region and lesion-free rules), the report metadata, and the data collection pointer. | LeukoGraph can read either version |
| **7. Second LeukoGraph run** | Same pipeline and analysis plan, new tables. Compared with the first run as a sensitivity analysis. Start only after the first results are in. | Second set of results |

## 4. Size and cost (estimates, to be measured in the pilot)

- Tract features: 1,512 now. Expect roughly 2,100 to 2,400 afterwards (rings, NAWM, `n_voxels`, difference and AUC, minus the 220
  ratio columns). The exact number comes from the pilot. This matters for autoencoder memory.
- Compute: the AUC adds a ranking step for about 11 maps x 42 tracts per subject; expected modest next to the existing work.
- Calendar: implementation and tests are a few working sessions; the pilot about a day; the full recompute depends on cluster load
  (measured in the pilot); the second LeukoGraph run takes about as long as the first (sweeps about 12 to 14 hours in this run).

## 5. Risks and how they are handled

| Risk | Handling |
|---|---|
| New run overwrites tables the experiment reads | Separate output root (Phase 1), tested before anything else |
| Two process-all runs on the same output folder can collide | Run the metrics step alone, after the other processing has settled |
| Silent value changes in "unchanged" metrics | Regression test against 2.0.0 outputs |
| Column explosion exhausts memory | Pilot reports the count first; reduced statistics for outer rings |
| Cluster contention with the running experiment | Concurrency cap; no priority changes |
| Fast path and legacy path disagree | New metrics get brute-force reference implementations in the tests |
| New columns unclassified downstream | Phase 6 checks that every column is classified before any run |

## 6. What does not change

The running LeukoGraph experiment and its analysis plan; the published release 2.0.0 and its DOI; the finished paper (apart from
the check in Phase 0).

## 7. Decisions

| # | Question | Decision (2026-10-17) |
|---|---|---|
| 1 | Streamline disconnection | include |
| 2 | Periventricular / deep lesion volumes | include |
| 3 | Delete `volume_tiv_normalised_mm3` | yes |
| 4 | Ring distances | include rings; distances 0-3, 3-6, 6-9 mm (recommended; to be re-confirmed) |
| 5 | Version number | 2.1.0 (the sentence continues "but ..."; awaiting the rest) |
| 6 | LeukoGraph run in progress | uses the existing metrics with the clip at +-5 SD (already in the running pipeline) as the stopgap for the Z-score change features |
