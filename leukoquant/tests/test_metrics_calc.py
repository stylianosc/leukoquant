"""Unit tests for metrics_calc module."""

import numpy as np
import pytest
from leukoquant.utils.metrics_calc import (
    _compute_load,
    _compute_tiv_from_array,
    _compute_whole_brain_lesion_microstructure_metrics_fast as _compute_whole_brain_lesion_microstructure_metrics,
    DISTANCE_CAP_MM,
    _distance_to_lesion,
    _summary_nonzero,
    _auc,
    _lesion_vs_reference,
    _parse_lesion_spec,
    _parse_named_spec,
)


class TestComputeLoad:
    """Test _compute_load function."""

    def test_empty_mask(self):
        """Test with empty mask (no overlap)."""
        mask = np.zeros((10, 10, 10), dtype=np.uint8)
        region = np.ones((10, 10, 10), dtype=np.uint8)
        voxel_vol = 1.0
        result = _compute_load(mask, region, voxel_vol)
        assert result == 0.0

    def test_full_overlap(self):
        """Test with complete overlap."""
        mask = np.ones((10, 10, 10), dtype=np.uint8)
        region = np.ones((10, 10, 10), dtype=np.uint8)
        voxel_vol = 2.0
        result = _compute_load(mask, region, voxel_vol)
        assert result == 1000 * 2.0  # 10*10*10 * 2.0

    def test_partial_overlap(self):
        """Test with partial overlap."""
        mask = np.zeros((10, 10, 10), dtype=np.uint8)
        mask[:5, :5, :5] = 1
        region = np.zeros((10, 10, 10), dtype=np.uint8)
        region[2:7, 2:7, 2:7] = 1
        voxel_vol = 1.0
        # Intersection is [2:5, 2:5, 2:5] = 3*3*3 = 27
        result = _compute_load(mask, region, voxel_vol)
        assert result == 27.0

    def test_different_voxel_volume(self):
        """Test with non-unit voxel volume."""
        mask = np.ones((5, 5, 5), dtype=np.uint8)
        region = np.ones((5, 5, 5), dtype=np.uint8)
        voxel_vol = 8.0  # 2x2x2 mm³
        result = _compute_load(mask, region, voxel_vol)
        assert result == 125 * 8.0  # 5*5*5 * 8.0


class TestComputeTiv:
    """Test _compute_tiv_from_array function."""

    def test_empty_brain(self):
        """Test with all-zero brain array."""
        brain = np.zeros((10, 10, 10), dtype=np.uint8)
        voxel_vol = 1.0
        result = _compute_tiv_from_array(brain, voxel_vol)
        assert result == 0.0

    def test_full_brain(self):
        """Test with fully filled brain."""
        brain = np.ones((10, 10, 10), dtype=np.uint8)
        voxel_vol = 1.0
        result = _compute_tiv_from_array(brain, voxel_vol)
        assert result == 1000.0

    def test_with_realistic_voxel_volume(self):
        """Test with realistic voxel volume (1mm³)."""
        brain = np.ones((100, 100, 100), dtype=np.uint8)
        voxel_vol = 1.0
        result = _compute_tiv_from_array(brain, voxel_vol)
        assert result == 1_000_000.0

    def test_partial_brain(self):
        """Test with partially filled brain."""
        brain = np.zeros((10, 10, 10), dtype=np.uint8)
        brain[:5, :5, :5] = 1
        voxel_vol = 8.0
        result = _compute_tiv_from_array(brain, voxel_vol)
        assert result == 125 * 8.0


class TestSummaryNonzero:
    """Test _summary_nonzero function."""

    def test_empty_array(self):
        """Test with all-zero array."""
        arr = np.array([0.0, 0.0, 0.0])
        result = _summary_nonzero(arr)
        assert result["n_voxels"] == 0.0
        assert result["p5"] == 0.0
        assert result["p95"] == 0.0
        assert result["mean"] == 0.0
        assert result["sum"] == 0.0
        assert not any(np.isnan(v) for v in result.values())

    def test_single_nonzero(self):
        """Test with one non-zero value."""
        arr = np.array([0.0, 5.0, 0.0])
        result = _summary_nonzero(arr)
        assert result["n_voxels"] == 1.0
        assert result["p5"] == 5.0
        assert result["p95"] == 5.0
        assert result["mean"] == 5.0
        assert result["sum"] == 5.0

    def test_multiple_nonzero(self):
        """Test with multiple non-zero values."""
        arr = np.array([0.0, 1.0, 2.0, 3.0, 0.0])
        result = _summary_nonzero(arr)
        assert result["n_voxels"] == 3.0
        assert result["p5"] == pytest.approx(1.1)       # linear interpolation of [1, 2, 3]
        assert result["p95"] == pytest.approx(2.9)
        assert result["mean"] == 2.0
        assert result["sum"] == 6.0

    def test_all_nonzero(self):
        """Test with all non-zero values."""
        arr = np.array([1.0, 2.0, 3.0, 4.0])
        result = _summary_nonzero(arr)
        assert result["n_voxels"] == 4.0
        assert result["p5"] == pytest.approx(1.15)
        assert result["p95"] == pytest.approx(3.85)
        assert result["peak_width"] == pytest.approx(2.7)
        assert result["IQR"] == pytest.approx(result["75th"] - result["25th"])
        assert result["mean"] == 2.5
        assert result["sum"] == 10.0

    def test_shape_stats_match_scipy(self):
        """Skewness and excess kurtosis equal scipy's defaults (bias=True, fisher=True) on the non-zero values."""
        scipy_stats = pytest.importorskip("scipy.stats")
        rng = np.random.default_rng(0)
        values = rng.gamma(2.0, 0.3, size=500)          # right-skewed, heavy-tailed
        arr = np.concatenate([values, np.zeros(100)])   # zeros are excluded, as for every statistic
        result = _summary_nonzero(arr)
        assert result["skew"] == pytest.approx(scipy_stats.skew(values), rel=1e-9)
        assert result["kurtosis"] == pytest.approx(scipy_stats.kurtosis(values), rel=1e-9)
        assert result["skew"] > 0

    def test_shape_stats_float32_input(self):
        """A float32 map gives the float64 moments (no precision loss in the third and fourth powers)."""
        scipy_stats = pytest.importorskip("scipy.stats")
        values = np.random.default_rng(1).normal(0.7, 0.05, size=2000).astype(np.float32)
        result = _summary_nonzero(values)
        assert result["skew"] == pytest.approx(scipy_stats.skew(values.astype(np.float64)), rel=1e-9)
        assert result["kurtosis"] == pytest.approx(scipy_stats.kurtosis(values.astype(np.float64)), rel=1e-9)

    def test_shape_stats_without_spread_are_zero(self):
        """Empty, single-value and constant regions have no shape: 0.0, never NaN."""
        for arr in (np.zeros(5), np.array([0.0, 5.0]), np.full(10, 0.3)):
            result = _summary_nonzero(arr)
            assert result["skew"] == 0.0
            assert result["kurtosis"] == 0.0


class TestLesionVersusReference:
    """Difference and AUC between the lesion and the reference tissue replace the (a-b)/(a+b) ratio."""

    def test_auc_worked_example(self):
        """Inside 1..5, outside 3..7: 25 pairs, 3 won (4>3, 5>3, 5>4), 3 tied (3,4,5), 19 lost -> (3 + 1.5) / 25 = 0.18."""
        assert _auc(np.arange(1.0, 6.0), np.arange(3.0, 8.0)) == pytest.approx(0.18)

    def test_auc_matches_brute_force_with_ties(self):
        rng = np.random.default_rng(5)
        x = rng.integers(0, 6, size=40).astype(float)
        y = rng.integers(0, 6, size=55).astype(float)
        brute = np.mean([(a > b) + 0.5 * (a == b) for a in x for b in y])
        assert _auc(x, y) == pytest.approx(brute)

    def test_auc_is_bounded_and_neutral_when_a_group_is_empty(self):
        assert _auc(np.array([5.0, 6.0]), np.array([1.0, 2.0])) == 1.0
        assert _auc(np.array([1.0, 2.0]), np.array([5.0, 6.0])) == 0.0
        assert _auc(np.array([]), np.array([1.0])) == 0.5
        assert _auc(np.array([1.0]), np.array([])) == 0.5

    def test_signed_values_that_cancel_do_not_explode(self):
        """Inside +0.51, outside -0.50: the old ratio was (a-b)/(a+b) = 101; the difference is simply 1.01."""
        out = _lesion_vs_reference(np.array([0.51]), np.array([-0.50]))
        assert out["mean_diff"] == pytest.approx(1.01)
        assert out["median_diff"] == pytest.approx(1.01)
        assert out["auc"] == 1.0
        assert all(np.isfinite(v) for v in out.values())

    def test_lesion_lower_than_reference_has_the_right_sign(self):
        """Z-scores -1.0 inside and -0.2 outside: the old ratio gave +0.67 (wrong sign); the difference is -0.8."""
        out = _lesion_vs_reference(np.array([-1.0, -1.0]), np.array([-0.2, -0.2, -0.2]))
        assert out["mean_diff"] == pytest.approx(-0.8)
        assert out["auc"] == 0.0

    def test_empty_groups_give_neutral_values(self):
        out = _lesion_vs_reference(np.zeros(3), np.array([1.0, 2.0]))
        assert out["auc"] == 0.5
        assert all(v == 0.0 for k, v in out.items() if k != "auc")
        assert set(out) == {f"{s}_diff" for s in ("p5", "25th", "median", "75th", "p95", "mean", "sum", "std", "IQR",
                                                   "peak_width")} | {"auc"}


class TestParseLesionSpec:
    """Test _parse_lesion_spec function."""

    def test_empty_spec(self):
        """Test with empty specification."""
        result = _parse_lesion_spec("")
        assert result == []

    def test_single_path_no_name(self):
        """Test single path without name uses default."""
        result = _parse_lesion_spec("/path/to/lesion.nii.gz")
        assert result == [("lesion", "/path/to/lesion.nii.gz")]

    def test_single_path_with_equals(self):
        """Test single path with name using equals."""
        result = _parse_lesion_spec("wmh=/path/to/wmh.nii.gz")
        assert result == [("wmh", "/path/to/wmh.nii.gz")]

    def test_single_path_with_colon(self):
        """Test single path with name using colon."""
        result = _parse_lesion_spec("la:/path/to/la.nii.gz")
        assert result == [("la", "/path/to/la.nii.gz")]

    def test_multiple_lesions(self):
        """Test multiple lesion specifications."""
        spec = "wmh=/data/wmh.nii.gz,la=/data/la.nii.gz"
        result = _parse_lesion_spec(spec)
        assert result == [
            ("wmh", "/data/wmh.nii.gz"),
            ("la", "/data/la.nii.gz"),
        ]

    def test_empty_path_raises(self):
        """Test that empty path raises ValueError."""
        with pytest.raises(ValueError, match="path must be non-empty"):
            _parse_lesion_spec("wmh=")

    def test_custom_default_name(self):
        """Test with custom default name."""
        result = _parse_lesion_spec("/path/to/lesion.nii.gz", default_name="custom")
        assert result == [("custom", "/path/to/lesion.nii.gz")]


class TestParseNamedSpec:
    """Test _parse_named_spec function."""

    def test_empty_spec(self):
        """Test with empty specification."""
        result = _parse_named_spec("")
        assert result == []

    def test_single_named_spec(self):
        """Test single named specification."""
        result = _parse_named_spec("fa=/path/to/fa.nii.gz")
        assert result == [("fa", "/path/to/fa.nii.gz")]

    def test_multiple_specs(self):
        """Test multiple specifications."""
        spec = "fa=/path/fa.nii.gz,md=/path/md.nii.gz"
        result = _parse_named_spec(spec)
        assert result == [
            ("fa", "/path/fa.nii.gz"),
            ("md", "/path/md.nii.gz"),
        ]

    def test_path_without_name_raises(self):
        """Test that path without name raises when no default provided."""
        with pytest.raises(ValueError, match="name must be non-empty"):
            _parse_named_spec("/path/to/file.nii.gz")

    def test_path_without_name_uses_default(self):
        """Test that path without name uses default."""
        result = _parse_named_spec("/path/to/file.nii.gz", default_name="map")
        assert result == [("map", "/path/to/file.nii.gz")]

    def test_empty_name_raises(self):
        """Test that empty name raises."""
        with pytest.raises(ValueError, match="name must be non-empty"):
            _parse_named_spec("=/path/to/file.nii.gz")

    def test_empty_path_raises(self):
        """Test that empty path raises."""
        with pytest.raises(ValueError, match="path must be non-empty"):
            _parse_named_spec("fa=")


class TestComputeWholeBrainLesionMicrostructureMetrics:
    """Test _compute_whole_brain_lesion_microstructure_metrics function."""

    def test_column_naming(self):
        """Column names follow {map}_{lesion}_map_wb_{region}_{stat}."""
        shape = (4, 4, 4)
        map_arrays = {"fa": np.ones(shape)}
        lesion_arrays = {"wmh": np.zeros(shape)}
        lesion_arrays["wmh"][0, 0, 0] = 0.8
        binary_lesion_arrays = {"wmh": (lesion_arrays["wmh"] > 0).astype(np.uint8)}
        penumbras_arrays = {"wmh": np.zeros(shape, dtype=np.uint8)}
        penumbras_arrays["wmh"][0, 0, 1] = 1
        dilated_lesion_arrays = {
            "wmh": binary_lesion_arrays["wmh"] | penumbras_arrays["wmh"]
        }

        df = _compute_whole_brain_lesion_microstructure_metrics(
            subject="sub-01",
            map_arrays=map_arrays,
            lesion_arrays=lesion_arrays,
            binary_lesion_arrays=binary_lesion_arrays,
            penumbras_arrays=penumbras_arrays,
            dilated_lesion_arrays=dilated_lesion_arrays,
        )

        assert len(df) == 1
        assert df.iloc[0]["subject"] == "sub-01"
        for region in ("lesion_binary", "lesion_pd", "penumbra", "dilated_lesion"):
            for stat in ("n_voxels", "p5", "p95", "mean", "sum", "median", "std", "25th", "75th",
                        "IQR", "peak_width"):
                assert f"fa_wmh_map_wb_{region}_{stat}" in df.columns

    def test_empty_lesion_gives_zero_stats(self):
        """An all-zero lesion mask yields the _summary_nonzero zero-fill, not NaN/crash."""
        shape = (4, 4, 4)
        map_arrays = {"fa": np.ones(shape)}
        lesion_arrays = {"wmh": np.zeros(shape)}
        binary_lesion_arrays = {"wmh": np.zeros(shape, dtype=np.uint8)}
        penumbras_arrays = {"wmh": np.zeros(shape, dtype=np.uint8)}
        dilated_lesion_arrays = {"wmh": np.zeros(shape, dtype=np.uint8)}

        df = _compute_whole_brain_lesion_microstructure_metrics(
            subject="sub-01",
            map_arrays=map_arrays,
            lesion_arrays=lesion_arrays,
            binary_lesion_arrays=binary_lesion_arrays,
            penumbras_arrays=penumbras_arrays,
            dilated_lesion_arrays=dilated_lesion_arrays,
        )
        assert df.iloc[0]["fa_wmh_map_wb_lesion_binary_mean"] == 0.0
        assert df.iloc[0]["fa_wmh_map_wb_lesion_binary_sum"] == 0.0

    def test_partial_overlap_values(self):
        """Values are map x region, summarised over the whole brain (no tract mask)."""
        shape = (2, 2, 2)
        fa = np.full(shape, 2.0)
        binary = np.zeros(shape, dtype=np.uint8)
        binary[0, 0, 0] = 1
        binary[0, 0, 1] = 1  # 2 lesion voxels, both map value 2.0

        df = _compute_whole_brain_lesion_microstructure_metrics(
            subject="sub-01",
            map_arrays={"fa": fa},
            lesion_arrays={"wmh": binary.astype(float)},
            binary_lesion_arrays={"wmh": binary},
            penumbras_arrays={"wmh": np.zeros(shape, dtype=np.uint8)},
            dilated_lesion_arrays={"wmh": binary},
        )
        row = df.iloc[0]
        assert row["fa_wmh_map_wb_lesion_binary_sum"] == 4.0  # 2 voxels x 2.0
        assert row["fa_wmh_map_wb_lesion_binary_mean"] == 2.0
        assert row["fa_wmh_map_wb_lesion_binary_p95"] == 2.0
        assert row["fa_wmh_map_wb_lesion_binary_p5"] == 2.0
        assert row["fa_wmh_map_wb_lesion_binary_n_voxels"] == 2.0

    def test_multiple_lesions_and_maps(self):
        """Every (map, lesion) pair gets its own column set -- no cross-contamination."""
        shape = (3, 3, 3)
        map_arrays = {"fa": np.ones(shape), "md": np.full(shape, 3.0)}
        binary_lesion_arrays = {
            "lesion_a": np.zeros(shape, dtype=np.uint8),
            "lesion_b": np.zeros(shape, dtype=np.uint8),
        }
        binary_lesion_arrays["lesion_a"][0, 0, 0] = 1
        binary_lesion_arrays["lesion_b"][1, 1, 1] = 1
        lesion_arrays = {k: v.astype(float) for k, v in binary_lesion_arrays.items()}
        penumbras_arrays = {k: np.zeros(shape, dtype=np.uint8) for k in binary_lesion_arrays}
        dilated_lesion_arrays = binary_lesion_arrays

        df = _compute_whole_brain_lesion_microstructure_metrics(
            subject="sub-01",
            map_arrays=map_arrays,
            lesion_arrays=lesion_arrays,
            binary_lesion_arrays=binary_lesion_arrays,
            penumbras_arrays=penumbras_arrays,
            dilated_lesion_arrays=dilated_lesion_arrays,
        )
        row = df.iloc[0]
        assert row["fa_lesion_a_map_wb_lesion_binary_sum"] == 1.0
        assert row["md_lesion_a_map_wb_lesion_binary_sum"] == 3.0
        assert row["fa_lesion_b_map_wb_lesion_binary_sum"] == 1.0
        assert row["md_lesion_b_map_wb_lesion_binary_sum"] == 3.0

    def test_missing_optional_region_arrays_no_crash(self):
        """A lesion with no penumbra/dilated entry just omits those columns."""
        shape = (2, 2, 2)
        binary = np.zeros(shape, dtype=np.uint8)
        binary[0, 0, 0] = 1

        df = _compute_whole_brain_lesion_microstructure_metrics(
            subject="sub-01",
            map_arrays={"fa": np.ones(shape)},
            lesion_arrays={},  # no pd array for this lesion
            binary_lesion_arrays={"wmh": binary},
            penumbras_arrays={},  # no penumbra array for this lesion
            dilated_lesion_arrays={},  # no dilated array for this lesion
        )
        assert "fa_wmh_map_wb_lesion_binary_sum" in df.columns
        assert "fa_wmh_map_wb_lesion_pd_sum" not in df.columns
        assert "fa_wmh_map_wb_penumbra_sum" not in df.columns
        assert "fa_wmh_map_wb_dilated_lesion_sum" not in df.columns

    def test_no_lesions_returns_subject_only_row(self):
        """No lesion masks at all still returns a valid single-row DataFrame."""
        df = _compute_whole_brain_lesion_microstructure_metrics(
            subject="sub-01",
            map_arrays={"fa": np.ones((2, 2, 2))},
            lesion_arrays={},
            binary_lesion_arrays={},
            penumbras_arrays={},
            dilated_lesion_arrays={},
        )
        assert len(df) == 1
        assert df.iloc[0]["subject"] == "sub-01"
        assert list(df.columns) == ["subject"]


class TestTractBinaryDilatedLesionStatistics:
    """The binary-tract x dilated-lesion statistics summarise the map, not the map squared.

    Until this fix both the fast and the legacy implementations multiplied the map in twice, so
    {map}_{lesion}_map_tract_binary_dilated_lesion_* held statistics of map**2. With a constant map of 2.0
    the mean over the region must be 2.0 (it was 4.0), the same value as the other lesion regions.
    """

    @staticmethod
    def _inputs():
        shape = (6, 6, 6)
        tract = np.zeros(shape)
        tract[1:5, 1:5, 1:5] = 1.0
        lesion = np.zeros(shape)
        lesion[2:4, 2:4, 2:4] = 1.0
        binary = (lesion > 0).astype(np.uint8)
        dilated = np.zeros(shape, dtype=np.uint8)
        dilated[1:5, 1:5, 1:5] = 1
        penumbra = (dilated - binary).astype(np.uint8)
        return dict(
            subject="sub-01",
            tract_arrays={"lh.cst": tract},
            tract_binary_arrays={"lh.cst": (tract > 0).astype(np.uint8)},
            tract_skeletons_arrays={"lh.cst": np.zeros(shape, dtype=np.uint8)},
            tract_source={"lh.cst": "tractography"},
            lesion_arrays={"wmh": lesion},
            binary_lesion_arrays={"wmh": binary},
            dilated_lesion_arrays={"wmh": dilated},
            penumbras_arrays={"wmh": penumbra},
            map_arrays={"fa": np.full(shape, 2.0)},
            voxel_vol=1.0,
            voxel_sizes=(1.0, 1.0, 1.0),
        )

    def test_fast_path_summarises_the_map(self):
        from leukoquant.utils.metrics_calc import _compute_tract_level_metrics_fast
        df = _compute_tract_level_metrics_fast(**self._inputs())
        row = df.iloc[0]
        assert row["fa_wmh_map_tract_binary_dilated_lesion_mean"] == pytest.approx(2.0)
        assert row["fa_wmh_map_tract_binary_dilated_lesion_p95"] == pytest.approx(2.0)
        assert row["fa_wmh_map_tract_binary_dilated_lesion_mean"] == pytest.approx(
            row["fa_wmh_map_tract_binary_penumbra_mean"])


class TestNonFiniteVoxels:
    """NaN / Inf voxels are missing data: set to 0 (excluded by every summary), never propagated."""

    def test_non_finite_voxels_become_zero_and_summaries_stay_finite(self):
        from leukoquant.utils.metrics_calc import _treat_non_finite_as_missing
        arr = np.array([1.0, np.nan, 3.0, np.inf, -np.inf, 5.0], dtype=np.float32)
        out = _treat_non_finite_as_missing(arr, "test map")
        np.testing.assert_array_equal(out, [1.0, 0.0, 3.0, 0.0, 0.0, 5.0])
        stats = _summary_nonzero(out)
        assert stats["n_voxels"] == 3.0
        assert all(np.isfinite(v) for v in stats.values())

    def test_integer_arrays_pass_through_untouched(self):
        from leukoquant.utils.metrics_calc import _treat_non_finite_as_missing
        arr = np.arange(5, dtype=np.int16)
        assert _treat_non_finite_as_missing(arr) is arr


class TestWeightedSummary:
    """_weighted_summary: values untouched, weights decide how much each voxel counts."""

    def test_equal_weights_equal_the_unweighted_summary(self):
        from leukoquant.utils.metrics_calc import _weighted_summary
        rng = np.random.default_rng(3)
        v = rng.gamma(2.0, 0.3, size=400)
        plain = _summary_nonzero(v)
        weighted = _weighted_summary(v, np.ones_like(v))
        for key in plain:
            assert weighted[key] == pytest.approx(plain[key], rel=1e-9), key

    def test_edge_voxel_is_not_shrunk(self):
        """FA 0.6 at tract probability 1.0 and FA 0.4 at probability 0.2: weighted mean 0.5667, product mean 0.34."""
        from leukoquant.utils.metrics_calc import _weighted_summary
        v = np.array([0.6, 0.4])
        w = np.array([1.0, 0.2])
        stats = _weighted_summary(v, w)
        assert stats["mean"] == pytest.approx((0.6 * 1.0 + 0.4 * 0.2) / 1.2)
        assert _summary_nonzero(v * w)["mean"] == pytest.approx(0.34)
        assert stats["sum"] == pytest.approx(0.68)              # sum(w v)
        assert stats["n_voxels"] == 2.0

    def test_integer_weights_equal_repeating_the_values(self):
        """Weighted mean, std, skew and kurtosis equal those of the values repeated weight times; percentiles are close."""
        from leukoquant.utils.metrics_calc import _weighted_summary
        rng = np.random.default_rng(4)
        v = rng.normal(1.0, 0.3, size=300)
        w = rng.integers(1, 6, size=300).astype(float)
        repeated = np.repeat(v, w.astype(int))
        stats = _weighted_summary(v, w)
        ref = _summary_nonzero(repeated)
        for key in ("mean", "std", "skew", "kurtosis"):
            assert stats[key] == pytest.approx(ref[key], rel=1e-9), key
        for key in ("p5", "25th", "median", "75th", "p95"):
            assert stats[key] == pytest.approx(ref[key], abs=0.06), key

    def test_zero_values_and_zero_weights_are_excluded_and_empty_is_zero(self):
        from leukoquant.utils.metrics_calc import _weighted_summary
        stats = _weighted_summary(np.array([0.0, 2.0, 3.0]), np.array([1.0, 0.0, 1.0]))
        assert stats["n_voxels"] == 1.0 and stats["mean"] == 3.0
        empty = _weighted_summary(np.zeros(4), np.ones(4))
        assert empty["n_voxels"] == 0.0 and all(v == 0.0 for v in empty.values())


class TestTissueAroundLesions:
    """Penumbra rings and NAWM are defined by distance in mm; the lesion-free case has no spurious distances."""

    @staticmethod
    def _inputs(with_lesion=True):
        shape = (24, 8, 8)
        tract = np.zeros(shape)
        tract[:, 2:6, 2:6] = 1.0                       # a bar along x, 24 voxels long
        lesion = np.zeros(shape)
        if with_lesion:
            lesion[0:2, 3:5, 3:5] = 1.0                # lesion at one end of the bar
        binary = (lesion > 0).astype(np.uint8)
        dist = _distance_to_lesion(binary, (1.0, 1.0, 1.0))
        dilated = (dist <= 3.0).astype(np.uint8)
        penumbra = ((dist > 0) & (dist <= 3.0)).astype(np.uint8)
        fa = np.zeros(shape)
        fa[...] = 0.5
        fa[dist <= 3.0] = 0.2                          # damaged tissue within 3 mm
        fa[(dist > 3.0) & (dist <= 6.0)] = 0.4         # milder change up to 6 mm
        return dict(
            subject="s", tract_arrays={"t": tract}, tract_binary_arrays={"t": (tract > 0).astype(np.uint8)},
            tract_skeletons_arrays={"t": np.zeros(shape, dtype=np.uint8)}, tract_source={"t": "tractography"},
            lesion_arrays={"wmh": lesion}, binary_lesion_arrays={"wmh": binary}, dilated_lesion_arrays={"wmh": dilated},
            penumbras_arrays={"wmh": penumbra}, map_arrays={"fa": fa}, voxel_vol=1.0, voxel_sizes=(1.0, 1.0, 1.0)), dist

    def test_distance_cap_and_lesion_free_subject(self):
        d = _distance_to_lesion(np.zeros((5, 5, 5)), (1.0, 1.0, 1.0))
        assert np.all(d == DISTANCE_CAP_MM)             # not the meaningless distance to the array border

    def test_rings_and_nawm_see_the_graded_change(self):
        from leukoquant.utils.metrics_calc import _compute_tract_level_metrics_fast
        kw, dist = self._inputs()
        row = _compute_tract_level_metrics_fast(**kw).iloc[0]
        assert row["fa_wmh_map_tract_binary_penumbra_mean"] == pytest.approx(0.2)           # within 3 mm
        assert row["fa_wmh_map_tract_binary_penumbra_3_6_mean"] == pytest.approx(0.4)       # 3 to 6 mm
        assert row["fa_wmh_map_tract_binary_penumbra_6_9_mean"] == pytest.approx(0.5)       # unchanged tissue
        # NAWM is everything of the tract farther than 3 mm from the lesion
        tract_mask = kw["tract_arrays"]["t"] > 0
        n_nawm = int((tract_mask & (dist > 3.0)).sum())
        assert row["wmh_tract_nawm_n_voxels"] == n_nawm
        assert row["fa_wmh_map_tract_binary_nawm_n_voxels"] == n_nawm
        assert row["fa_wmh_map_tract_binary_nawm_p5"] >= 0.4 - 1e-9                         # never the damaged 0.2

    def test_lesion_versus_nawm_columns_in_the_tract_table(self):
        from leukoquant.utils.metrics_calc import _compute_tract_level_metrics_fast
        kw, _ = self._inputs()
        row = _compute_tract_level_metrics_fast(**kw).iloc[0]
        # FA is 0.2 inside the lesion; the NAWM holds 0.4 (3 to 6 mm) and 0.5 (beyond)
        assert row["fa_wmh_map_tract_lesion_vs_nawm_auc"] == 0.0
        assert row["fa_wmh_map_tract_lesion_vs_nawm_mean_diff"] == pytest.approx(0.2 - row["fa_wmh_map_tract_binary_nawm_mean"])
        assert row["fa_wmh_map_tract_lesion_vs_nawm_median_diff"] < 0
        assert row["wmh_tract_lesion_vs_nawm_auc"] == pytest.approx(0.5)                # tract probability is 1 everywhere
        assert not [c for c in row.index if c.endswith("_change") or c.endswith("_change_absolute")]

    def test_share_of_tract_density_inside_lesions(self):
        from leukoquant.utils.metrics_calc import _compute_tract_level_metrics_fast
        kw, _ = self._inputs()
        # tract: 24 x 4 x 4 voxels of density 1; lesion: 2 x 2 x 2 voxels inside it
        row = _compute_tract_level_metrics_fast(**kw).iloc[0]
        assert row["wmh_tract_density_fraction_in_lesion"] == pytest.approx(8 / 384)
        kw["tract_arrays"]["t"][0:2, 3:5, 3:5] = 3.0        # density three times higher exactly where the lesion is
        row = _compute_tract_level_metrics_fast(**kw).iloc[0]
        assert row["wmh_tract_density_fraction_in_lesion"] == pytest.approx(24 / (24 + 376))
        kw_free, _ = self._inputs(with_lesion=False)
        assert _compute_tract_level_metrics_fast(**kw_free).iloc[0]["wmh_tract_density_fraction_in_lesion"] == 0.0

    def test_lesion_free_subject_has_all_tract_as_nawm_and_empty_rings(self):
        from leukoquant.utils.metrics_calc import _compute_tract_level_metrics_fast
        kw, _ = self._inputs(with_lesion=False)
        row = _compute_tract_level_metrics_fast(**kw).iloc[0]
        assert row["fa_wmh_map_tract_binary_nawm_n_voxels"] == row["tract_volume_mm3"]
        assert row["fa_wmh_map_tract_binary_penumbra_n_voxels"] == 0.0
        assert row["fa_wmh_map_tract_binary_penumbra_3_6_n_voxels"] == 0.0
        assert row["wmh_tract_nearest_lesion_distance"] == DISTANCE_CAP_MM


class TestTractLesionDistance:
    """Distances of a tract to the lesions: nearest, mean, and the share within 5 mm."""

    def test_distances_of_a_bar_next_to_a_lesion(self):
        from leukoquant.utils.metrics_calc import _compute_tract_level_metrics_fast
        shape = (30, 8, 8)
        tract = np.zeros(shape); tract[10:20, 2:6, 2:6] = 1.0           # tract voxels at x = 10..19
        lesion = np.zeros(shape); lesion[0:6, 2:6, 2:6] = 1.0           # lesion at x = 0..5, so the tract starts 5 mm away
        binary = (lesion > 0).astype(np.uint8)
        dist = _distance_to_lesion(binary, (1.0, 1.0, 1.0))
        kw = dict(subject="s", tract_arrays={"t": tract}, tract_binary_arrays={"t": (tract > 0).astype(np.uint8)},
                  tract_skeletons_arrays={"t": np.zeros(shape, dtype=np.uint8)}, tract_source={"t": "tractography"},
                  lesion_arrays={"wmh": lesion}, binary_lesion_arrays={"wmh": binary},
                  dilated_lesion_arrays={"wmh": (dist <= 3).astype(np.uint8)},
                  penumbras_arrays={"wmh": ((dist > 0) & (dist <= 3)).astype(np.uint8)},
                  map_arrays={"fa": np.full(shape, 0.5)}, voxel_vol=1.0, voxel_sizes=(1.0, 1.0, 1.0))
        row = _compute_tract_level_metrics_fast(**kw).iloc[0]
        assert row["wmh_tract_nearest_lesion_distance"] == pytest.approx(5.0)       # x = 10 is 5 voxels from x = 5
        assert row["wmh_tract_mean_lesion_distance"] == pytest.approx(9.5)          # mean of 5..14
        assert row["wmh_tract_fraction_within_5mm_of_lesion"] == pytest.approx(0.1) # only the x = 10 slab (distance 5.0)
        assert "wmh_tract_centre_of_mass_nearest_lesion_distance" not in row.index

    def test_a_tract_containing_lesion_has_nearest_distance_zero(self):
        from leukoquant.utils.metrics_calc import _compute_tract_level_metrics_fast
        kw, _ = TestTissueAroundLesions._inputs()
        row = _compute_tract_level_metrics_fast(**kw).iloc[0]
        assert row["wmh_tract_nearest_lesion_distance"] == 0.0
        assert 0.0 < row["wmh_tract_fraction_within_5mm_of_lesion"] < 1.0


class TestWholeBrainLesionCounts:
    def test_count_and_mean_size_of_separate_lesions(self):
        from leukoquant.utils.metrics_calc import _compute_whole_brain_metrics
        shape = (20, 10, 10)
        brain = np.ones(shape)
        lesion = np.zeros(shape)
        lesion[1:3, 1:3, 1:3] = 1          # 8 voxels
        lesion[8:12, 1:4, 1:3] = 1         # 24 voxels, not touching the first (also not diagonally)
        row = _compute_whole_brain_metrics("s", brain, {"wmh": lesion}, 1.0).iloc[0]
        assert row["wmh_n_lesions"] == 2
        assert row["wmh_mean_lesion_volume_mm3"] == pytest.approx(16.0)
        assert row["wmh_load_pct"] == pytest.approx(32 / 2000 * 100)
        assert "wmh_volume_tiv_normalised_mm3" not in row.index

    def test_diagonal_neighbours_are_one_lesion_and_no_lesion_is_zero(self):
        from leukoquant.utils.metrics_calc import _compute_whole_brain_metrics
        shape = (6, 6, 6)
        lesion = np.zeros(shape); lesion[1, 1, 1] = 1; lesion[2, 2, 2] = 1       # touch at a corner: one 26-connected lesion
        row = _compute_whole_brain_metrics("s", np.ones(shape), {"wmh": lesion}, 1.0).iloc[0]
        assert row["wmh_n_lesions"] == 1
        empty = _compute_whole_brain_metrics("s", np.ones(shape), {"wmh": np.zeros(shape)}, 1.0).iloc[0]
        assert empty["wmh_n_lesions"] == 0 and empty["wmh_mean_lesion_volume_mm3"] == 0.0


class TestPeriventricularDeep:
    def test_split_at_ten_millimetres(self):
        from leukoquant.utils.metrics_calc import _periventricular_deep_metrics
        shape = (50, 10, 10)
        ventricles = np.zeros(shape); ventricles[0:5, :, :] = 1               # ventricle slab at x = 0..4
        lesion = np.zeros(shape)
        lesion[8:10, 2:4, 2:4] = 1          # x = 8,9: 4-5 mm from the ventricle -> periventricular (8 voxels)
        lesion[30:33, 2:4, 2:4] = 1         # x = 30..32: 26-28 mm away -> deep (12 voxels)
        out = _periventricular_deep_metrics("wmh", lesion, ventricles, (1.0, 1.0, 1.0), 2.0)
        assert out["wmh_pvwmh_volume_mm3"] == pytest.approx(8 * 2.0)
        assert out["wmh_dwmh_volume_mm3"] == pytest.approx(12 * 2.0)
        assert out["wmh_pvwmh_fraction"] == pytest.approx(0.4)

    def test_no_lesion_and_no_ventricle_are_zero_not_nan(self):
        from leukoquant.utils.metrics_calc import _periventricular_deep_metrics
        shape = (10, 10, 10)
        out = _periventricular_deep_metrics("wmh", np.zeros(shape), np.ones(shape), (1.0, 1.0, 1.0), 1.0)
        assert out == {"wmh_pvwmh_volume_mm3": 0.0, "wmh_dwmh_volume_mm3": 0.0, "wmh_pvwmh_fraction": 0.0}
        lesion = np.zeros(shape); lesion[2:4, 2:4, 2:4] = 1
        out = _periventricular_deep_metrics("wmh", lesion, np.zeros(shape), (1.0, 1.0, 1.0), 1.0)
        assert out["wmh_pvwmh_volume_mm3"] == 0.0 and out["wmh_dwmh_volume_mm3"] == 8.0   # nothing is near a missing ventricle

    def test_shape_mismatch_is_an_error(self):
        from leukoquant.utils.metrics_calc import _periventricular_deep_metrics
        with pytest.raises(ValueError):
            _periventricular_deep_metrics("wmh", np.zeros((4, 4, 4)), np.zeros((5, 4, 4)), (1.0, 1.0, 1.0), 1.0)


class TestVentricleMask:
    @staticmethod
    def _save(arr, affine, path):
        import nibabel as nib
        nib.save(nib.Nifti1Image(arr.astype(np.int16), affine), str(path))

    def test_labels_become_a_mask_on_the_same_grid(self, tmp_path):
        import nibabel as nib
        from leukoquant.utils.metrics_calc import _ventricle_mask
        labels = np.zeros((10, 10, 10)); labels[2:4, 2:4, 2:4] = 52; labels[6:8, 6:8, 6:8] = 99
        self._save(labels, np.eye(4), tmp_path / "p.nii.gz")
        ref = nib.Nifti1Image(np.zeros((10, 10, 10)), np.eye(4))
        mask = _ventricle_mask(str(tmp_path / "p.nii.gz"), "50,51,52,53", ref)
        assert mask.sum() == 8 and mask[2:4, 2:4, 2:4].all()              # label 99 is not a ventricle

    def test_a_parcellation_on_a_finer_grid_is_resampled_to_the_reference(self, tmp_path):
        import nibabel as nib
        from leukoquant.utils.metrics_calc import _ventricle_mask
        labels = np.zeros((20, 20, 20)); labels[4:8, 4:8, 4:8] = 4         # 0.5 mm voxels: a 2 mm cube at 2..4 mm
        self._save(labels, np.diag([0.5, 0.5, 0.5, 1.0]), tmp_path / "p.nii.gz")
        ref = nib.Nifti1Image(np.zeros((10, 10, 10)), np.eye(4))
        mask = _ventricle_mask(str(tmp_path / "p.nii.gz"), "4,43", ref)
        assert mask.shape == (10, 10, 10) and 0 < mask.sum() <= 27

    def test_missing_labels_are_an_error_not_a_silent_all_deep(self, tmp_path):
        import nibabel as nib
        from leukoquant.utils.metrics_calc import _ventricle_mask
        self._save(np.ones((6, 6, 6)), np.eye(4), tmp_path / "p.nii.gz")
        with pytest.raises(ValueError):
            _ventricle_mask(str(tmp_path / "p.nii.gz"), "50,51", nib.Nifti1Image(np.zeros((6, 6, 6)), np.eye(4)))


class TestSkeletonIsBinary:
    """Skeleton statistics (PSMD is the peak width of MD on the skeleton) use the map values on the skeleton, not map x weight."""

    @staticmethod
    def _md_map():
        md = np.full((10, 10, 10), 0.0008)
        md[5, 5, :] += np.linspace(-0.0002, 0.0002, 10)                  # the skeleton line carries the spread
        return md

    def test_fractional_and_255_valued_skeletons_give_the_same_statistics_as_a_binary_one(self):
        from leukoquant.utils.metrics_calc import _skeleton_metrics_fast
        md = self._md_map()
        binary = np.zeros((10, 10, 10), np.uint8); binary[5, 5, :] = 1
        reference = _skeleton_metrics_fast("s", {"md": md}, binary).iloc[0]
        blurred = binary.astype(float)
        blurred[5, 6, :] = 0.3; blurred[5, 4, :] = 0.3                   # a linearly resampled thin mask: 1 on the line, a faint halo
        for weighted in (binary * 255, blurred):                          # scikit-image's uint8 255; the blurred mask
            out = _skeleton_metrics_fast("s", {"md": md}, weighted).iloc[0]
            assert out["skeleton_mean"] == pytest.approx(reference["skeleton_mean"])
            assert out["skeleton_peak_width"] == pytest.approx(reference["skeleton_peak_width"])
        assert reference["skeleton_mean"] == pytest.approx(md[5, 5, :].mean())

    def test_voxels_with_a_faint_blur_are_not_on_the_skeleton(self):
        from leukoquant.utils.metrics_calc import _skeleton_metrics_fast
        md = self._md_map()
        sk = np.zeros((10, 10, 10)); sk[5, 5, :] = 1.0; sk[5, 6, :] = 0.2    # blur next to the line
        out = _skeleton_metrics_fast("s", {"md": md}, sk).iloc[0]
        assert out["skeleton_n_voxels"] == 10
