"""Smoke tests for metrics_calc module.

These tests verify basic functionality without requiring full NIfTI fixtures.
"""

import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import nibabel as nib
import pytest
from leukoquant.utils.metrics_calc import (
    _load,
    _voxel_vol_mm3,
    _assert_shape_compatibility,
)


@pytest.fixture
def temp_dir():
    """Create temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def simple_nifti(temp_dir):
    """Create a simple test NIfTI file."""
    data = np.random.rand(10, 10, 10)
    affine = np.diag([2.0, 2.0, 2.0, 1.0])  # 2mm isotropic
    img = nib.Nifti1Image(data, affine=affine)
    path = temp_dir / "test.nii.gz"
    nib.save(img, path)
    return path


class TestLoadFunction:
    """Test _load function."""

    def test_load_valid_nifti(self, simple_nifti):
        """Test loading a valid NIfTI file."""
        img = _load(str(simple_nifti))
        assert img is not None
        assert img.shape == (10, 10, 10)

    def test_load_nonexistent_raises(self):
        """Test that loading nonexistent file raises error."""
        with pytest.raises(Exception):  # FileNotFoundError or nibabel error
            _load("/nonexistent/file.nii.gz")


class TestVoxelVolFunction:
    """Test _voxel_vol_mm3 function."""

    def test_isotropic_voxel(self, simple_nifti):
        """Test voxel volume calculation for isotropic voxels."""
        img = _load(str(simple_nifti))
        voxel_vol = _voxel_vol_mm3(img)
        # 2mm x 2mm x 2mm = 8 mm³
        assert abs(voxel_vol - 8.0) < 0.01

    def test_anisotropic_voxel(self, temp_dir):
        """Test voxel volume calculation for anisotropic voxels."""
        data = np.zeros((5, 5, 5))
        affine = np.diag([1.0, 2.0, 4.0, 1.0])  # 1mm x 2mm x 4mm
        img = nib.Nifti1Image(data, affine=affine)
        path = temp_dir / "aniso.nii.gz"
        nib.save(img, path)

        loaded = _load(str(path))
        voxel_vol = _voxel_vol_mm3(loaded)
        # 1 x 2 x 4 = 8 mm³
        assert abs(voxel_vol - 8.0) < 0.01


class TestAssertShapeCompatibility:
    """Test _assert_shape_compatibility function."""

    def test_compatible_shapes(self):
        """Test arrays with matching shapes."""
        arrays = {
            "arr1": np.zeros((5, 5, 5)),
            "arr2": np.ones((5, 5, 5)),
        }
        # Should not raise
        _assert_shape_compatibility((5, 5, 5), arrays, "test_group")

    def test_incompatible_shapes_raises(self):
        """Test arrays with mismatched shapes raise ValueError."""
        arrays = {
            "arr1": np.zeros((5, 5, 5)),
            "arr2": np.ones((6, 6, 6)),
        }
        with pytest.raises(ValueError, match="Shape mismatch detected"):
            _assert_shape_compatibility((5, 5, 5), arrays, "test_group")

    def test_empty_arrays_dict(self):
        """Test with empty arrays dictionary."""
        arrays = {}
        # Should not raise
        _assert_shape_compatibility((5, 5, 5), arrays, "test_group")

    def test_single_mismatched_raises(self):
        """Test that single mismatched array raises."""
        arrays = {
            "correct": np.zeros((5, 5, 5)),
            "wrong": np.ones((4, 4, 4)),
        }
        with pytest.raises(ValueError, match="Shape mismatch detected"):
            _assert_shape_compatibility((5, 5, 5), arrays, "test_group")
