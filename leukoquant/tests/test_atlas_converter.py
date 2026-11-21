"""Unit tests for atlas_converter module."""

import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import nibabel as nib
import pytest
from leukoquant.utils.atlas_converter import convert_atlas_to_freesurfer


@pytest.fixture
def temp_dir():
    """Create temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def simple_parcellation(temp_dir):
    """Create a simple test parcellation NIfTI file."""
    data = np.zeros((10, 10, 10), dtype=np.int32)
    # Create regions with labels 1 and 2
    data[0:5, 0:5, 0:5] = 1
    data[5:10, 5:10, 5:10] = 2

    img = nib.Nifti1Image(data, affine=np.eye(4))
    path = temp_dir / "input.nii.gz"
    nib.save(img, path)
    return path


@pytest.fixture
def one_to_one_mapping(temp_dir):
    """Create a 1:1 mapping CSV."""
    mapping_data = {
        "Label": [1, 2],
        "Source": [1, 2],
        "Target": [10, 20],
        "Nearest": [10, 20],
    }
    df = pd.DataFrame(mapping_data)
    path = temp_dir / "mapping.csv"
    df.to_csv(path, index=False)
    return path


class TestConvertAtlas:
    """Test atlas conversion function."""

    def test_one_to_one_mapping(self, simple_parcellation, one_to_one_mapping, temp_dir):
        """Test basic 1:1 mapping conversion."""
        output = temp_dir / "output.nii.gz"
        result_path = convert_atlas_to_freesurfer(
            simple_parcellation,
            one_to_one_mapping,
            output_path=output,
            validate=False,
            verbose=False,
        )

        assert Path(result_path).exists()

        # Verify output
        output_img = nib.load(result_path)
        output_data = output_img.get_fdata().astype(np.int32)

        # First region should be remapped to 10
        assert np.all(output_data[0:5, 0:5, 0:5] == 10)
        # Second region should be remapped to 20
        assert np.all(output_data[5:10, 5:10, 5:10] == 20)

    def test_output_generated_if_none(self, simple_parcellation, one_to_one_mapping):
        """Test that output path is auto-generated if None."""
        result_path = convert_atlas_to_freesurfer(
            simple_parcellation,
            one_to_one_mapping,
            output_path=None,
            validate=False,
            verbose=False,
        )

        assert Path(result_path).exists()
        # Should have generated filename
        assert "converted" in result_path or "nii.gz" in result_path

    def test_missing_input_raises(self, one_to_one_mapping):
        """Test that missing input file raises error."""
        with pytest.raises(Exception):  # FileNotFoundError or nibabel error
            convert_atlas_to_freesurfer(
                "/nonexistent/file.nii.gz",
                one_to_one_mapping,
                validate=False,
            )

    def test_missing_mapping_raises(self, simple_parcellation):
        """Test that missing mapping file raises error."""
        with pytest.raises(Exception):  # FileNotFoundError or pandas error
            convert_atlas_to_freesurfer(
                simple_parcellation,
                "/nonexistent/mapping.csv",
                validate=False,
            )
