"""End-to-end run of metrics_calc.main() on a small synthetic subject: columns, ranges and the absence of missing values."""

import sys

import nibabel as nib
import numpy as np
import pandas as pd
import pytest

from leukoquant.utils import metrics_calc

SHAPE = (40, 24, 24)


def _save(arr, path, dtype=np.float32):
    nib.save(nib.Nifti1Image(np.asarray(arr).astype(dtype), np.eye(4)), str(path))
    return path


@pytest.fixture
def subject(tmp_path):
    rng = np.random.default_rng(0)
    t1 = _save(rng.random(SHAPE), tmp_path / "t1.nii.gz")
    brain = _save(np.ones(SHAPE), tmp_path / "brain.nii.gz")

    lesion = np.zeros(SHAPE)
    lesion[5:9, 8:12, 8:12] = 1                                   # near the ventricle slab (periventricular)
    lesion[28:32, 8:12, 8:12] = 1                                 # deep
    (tmp_path / "lesion").mkdir()
    _save(lesion, tmp_path / "lesion" / "wmh.nii.gz")

    tract = np.zeros(SHAPE)
    tract[2:38, 6:14, 6:14] = rng.uniform(0.2, 1.0, size=(36, 8, 8))     # a probabilistic bar through both lesions
    (tmp_path / "tracts" / "cst_avg16_syn").mkdir(parents=True)
    _save(tract, tmp_path / "tracts" / "cst_avg16_syn" / "path.pd.nii.gz")

    # a map that is lower inside the lesions, and a signed (Z-score-like) one
    fa = 0.5 + 0.05 * rng.standard_normal(SHAPE)
    fa[lesion > 0] -= 0.2
    z = rng.standard_normal(SHAPE)
    z[lesion > 0] -= 1.0
    (tmp_path / "maps").mkdir()
    _save(fa, tmp_path / "maps" / "fa.nii.gz")
    _save(z, tmp_path / "maps" / "fa_z.nii.gz")

    skeleton = np.zeros(SHAPE, np.uint8)
    skeleton[3:37, 10, 10] = 1
    (tmp_path / "skeleton").mkdir()
    _save(skeleton, tmp_path / "skeleton" / "skeleton.nii.gz", np.uint8)

    labels = np.zeros(SHAPE)
    labels[0:3, :, :] = 52                                        # "ventricle" slab at x = 0..2
    _save(labels, tmp_path / "parcellation.nii.gz")

    out = tmp_path / "out"
    argv = ["metrics_calc.py", "--subject", "sub-01", "--t1-path", str(t1), "--t1-brain-path", str(brain),
            "--lesion-path", f"wmh={tmp_path / 'lesion'}", "--tracts-path", str(tmp_path / "tracts"),
            "--tract-mode", "tractography", "--maps", f"fa={tmp_path / 'maps' / 'fa.nii.gz'},fa_z={tmp_path / 'maps' / 'fa_z.nii.gz'}",
            "--skeleton-path", str(tmp_path / "skeleton"), "--parcellation-path", str(tmp_path / "parcellation.nii.gz"),
            "--ventricle-labels", "50,51,52,53", "--output", str(out)]
    return argv, out / "metrics"


def _run(monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", argv)
    assert metrics_calc.main() == 0


def test_all_tables_are_written_finite_and_in_range(subject, monkeypatch):
    argv, metrics_dir = subject
    _run(monkeypatch, argv)

    tract = pd.read_csv(metrics_dir / "tract_level_metrics.csv")
    assert len(tract) == 1 and tract.loc[0, "tract_name"] == "cst"
    numeric = tract.select_dtypes("number")
    assert not numeric.isna().any().any() and np.isfinite(numeric.to_numpy()).all()       # no missing values, no Inf

    # the replaced change features
    assert not [c for c in tract.columns if c.endswith("_change") or c.endswith("_change_absolute")]
    aucs = tract[[c for c in tract.columns if c.endswith("lesion_vs_nawm_auc")]]
    assert len(aucs.columns) == 3                                  # the tract density, fa and fa_z
    assert ((aucs >= 0) & (aucs <= 1)).all().all()
    assert tract.loc[0, "fa_wmh_map_tract_lesion_vs_nawm_auc"] < 0.1                     # FA is clearly lower inside lesions
    assert tract.loc[0, "fa_z_wmh_map_tract_lesion_vs_nawm_mean_diff"] == pytest.approx(-1.0, abs=0.4)   # about 1 SD lower, finite

    # new statistics and regions
    for c in ("fa_wmh_map_tract_binary_penumbra_p5", "fa_wmh_map_tract_binary_penumbra_n_voxels",
              "fa_wmh_map_tract_binary_penumbra_3_6_mean", "fa_wmh_map_tract_binary_nawm_p95",
              "wmh_tract_mean_lesion_distance", "wmh_tract_fraction_within_5mm_of_lesion", "wmh_tract_density_fraction_in_lesion"):
        assert c in tract.columns, c
    assert not [c for c in tract.columns if c.endswith("_min") or c.endswith("_max")]
    assert "wmh_tract_centre_of_mass_nearest_lesion_distance" not in tract.columns

    whole = pd.read_csv(metrics_dir / "whole_brain_metrics.csv").iloc[0]
    assert whole["wmh_n_lesions"] == 2
    assert whole["wmh_pvwmh_volume_mm3"] == pytest.approx(4 * 4 * 4)                     # the lesion 2 voxels from the slab
    assert whole["wmh_dwmh_volume_mm3"] == pytest.approx(4 * 4 * 4)
    assert whole["wmh_pvwmh_fraction"] == pytest.approx(0.5)
    assert "wmh_volume_tiv_normalised_mm3" not in whole.index

    skeleton = pd.read_csv(metrics_dir / "skeleton_metrics.csv").set_index("map_name")
    assert skeleton.loc["fa", "skeleton_n_voxels"] == 34
    assert 0.35 < skeleton.loc["fa", "skeleton_mean"] < 0.65                              # a real FA, not FA x a fraction

    lesions = pd.read_csv(metrics_dir / "wmh_lesion_metrics.csv")
    assert len(lesions) == 2
    assert {"n_lesions_merged_after_dilation", "max_inscribed_radius_mm", "tract_match_2", "tract_match_overlap_fraction"} <= set(lesions.columns)
    assert (lesions["tract_match"] == "cst").all() and (lesions["tract_match_overlap_fraction"] == 1.0).all()


def test_runs_without_a_parcellation(subject, monkeypatch):
    argv, metrics_dir = subject
    i = argv.index("--parcellation-path")
    _run(monkeypatch, argv[:i] + argv[i + 4:])                     # drop --parcellation-path and --ventricle-labels (+ values)
    whole = pd.read_csv(metrics_dir / "whole_brain_metrics.csv").iloc[0]
    assert "wmh_pvwmh_volume_mm3" not in whole.index and whole["wmh_n_lesions"] == 2
