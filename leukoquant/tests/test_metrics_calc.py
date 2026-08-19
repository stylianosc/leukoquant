"""Unit tests for metrics_calc module."""

import numpy as np
import pytest
from leukoquant.utils.metrics_calc import (
    _compute_load,
    _compute_tiv_from_array,
    _summary_nonzero,
    _summary_change_nonzero,
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
        assert result["min"] == 0.0
        assert result["max"] == 0.0
        assert result["mean"] == 0.0
        assert result["sum"] == 0.0

    def test_single_nonzero(self):
        """Test with one non-zero value."""
        arr = np.array([0.0, 5.0, 0.0])
        result = _summary_nonzero(arr)
        assert result["min"] == 5.0
        assert result["max"] == 5.0
        assert result["mean"] == 5.0
        assert result["sum"] == 5.0

    def test_multiple_nonzero(self):
        """Test with multiple non-zero values."""
        arr = np.array([0.0, 1.0, 2.0, 3.0, 0.0])
        result = _summary_nonzero(arr)
        assert result["min"] == 1.0
        assert result["max"] == 3.0
        assert result["mean"] == 2.0
        assert result["sum"] == 6.0

    def test_all_nonzero(self):
        """Test with all non-zero values."""
        arr = np.array([1.0, 2.0, 3.0, 4.0])
        result = _summary_nonzero(arr)
        assert result["min"] == 1.0
        assert result["max"] == 4.0
        assert result["mean"] == 2.5
        assert result["sum"] == 10.0


class TestSummaryChangeNonzero:
    """Test _summary_change_nonzero function."""

    def test_returns_dict_with_change_keys(self):
        """Test that function returns dict with change metrics."""
        lesion_map = np.array([1.0, 1.0, 0.0])
        values = np.array([1.0, 2.0, 3.0])
        result = _summary_change_nonzero(lesion_map, values)

        # Function returns keys like "min_change", "max_change", etc.
        assert isinstance(result, dict)
        # Verify it has change-related keys
        change_keys = [k for k in result.keys() if "change" in k]
        assert len(change_keys) > 0

    def test_binary_lesion_map(self):
        """Test with binary lesion map."""
        lesion_map = np.array([1.0, 1.0, 1.0, 0.0, 0.0])
        values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result = _summary_change_nonzero(lesion_map, values)

        # Should contain change statistics
        assert isinstance(result, dict)
        assert len(result) > 0


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
