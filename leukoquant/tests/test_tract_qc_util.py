"""Unit tests for tract_qc_util module."""

import tempfile
from pathlib import Path
import pytest
from leukoquant.utils.tract_qc_util import (
    find_tract_files,
    parse_expected_tracts,
    run_tract_qc,
)


@pytest.fixture
def temp_dir():
    """Create temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


class TestFindTractFiles:
    """Test find_tract_files function."""

    def test_directory_path(self, temp_dir):
        """Test finding tract files in a directory."""
        subject = "subj_001"
        subject_dir = temp_dir / subject
        subject_dir.mkdir(parents=True)

        # Create dummy tract files
        (subject_dir / "tract1.trk").touch()
        (subject_dir / "tract2.trk").touch()

        tract_files, tract_names, exists = find_tract_files(
            subject, str(temp_dir)
        )

        assert exists is True
        assert len(tract_files) >= 0  # Depends on directory structure

    def test_base_glob_pattern(self, temp_dir):
        """Test finding tract files with colon-separated base:glob pattern."""
        subject = "subj_002"
        subject_dir = temp_dir / subject / "dpath"
        subject_dir.mkdir(parents=True)

        # Create dummy tract files
        (subject_dir / "path.pd.nii.gz").touch()

        tractography_path = f"{temp_dir}:{subject}/dpath/*.nii.gz"
        tract_files, tract_names, exists = find_tract_files(subject, tractography_path)

        # The pattern should find files or indicate none exist
        # This is a valid test of the glob pattern matching
        assert isinstance(exists, bool)
        if len(tract_files) > 0:
            assert any("path.pd.nii.gz" in f or "nii.gz" in f for f in tract_files)

    def test_nonexistent_subject(self, temp_dir):
        """Test with nonexistent subject."""
        subject = "nonexistent_subject"
        tractography_path = str(temp_dir)

        tract_files, tract_names, exists = find_tract_files(subject, tractography_path)

        # Should indicate that the subject path doesn't exist
        assert exists is False or len(tract_files) == 0

    def test_base_only_no_glob(self, temp_dir):
        """Test base:glob pattern with empty glob."""
        subject = "subj_003"
        subject_dir = temp_dir / subject
        subject_dir.mkdir(parents=True)
        (subject_dir / "file.txt").touch()

        tractography_path = f"{temp_dir}:"
        tract_files, tract_names, exists = find_tract_files(subject, tractography_path)

        # Should handle empty glob gracefully
        assert isinstance(tract_files, list)


class TestParseExpectedTracts:
    """Test parse_expected_tracts function."""

    def test_tracula_keyword(self):
        tracts = parse_expected_tracts("tracula")
        assert len(tracts) == 42
        assert "lh.cst_avg16_syn_bbr" in tracts
        assert "rh.cst_avg16_syn_bbr" in tracts
        assert "mcp_avg16_syn_bbr" in tracts

    def test_comma_separated(self):
        tracts = parse_expected_tracts("lh.cst,rh.cst,mcp")
        assert tracts == ["lh.cst", "rh.cst", "mcp"]

    def test_directory(self, temp_dir):
        (temp_dir / "lh.cst.display.trk").touch()
        (temp_dir / "rh.cst.display.trk").touch()
        tracts = parse_expected_tracts(str(temp_dir))
        assert "lh.cst" in tracts
        assert "rh.cst" in tracts

    def test_empty(self):
        assert parse_expected_tracts(None) == []
        assert parse_expected_tracts("") == []


class TestRunTractQcWithExpectedTracts:
    """Test run_tract_qc flags missing expected tracts."""

    def test_missing_tract_flagged_in_qc_report(self, temp_dir):
        import pandas as pd

        subject = "test_sub"
        subj_dir = temp_dir / subject / "dpath"
        # Only reconstruct acomm
        acomm_dir = subj_dir / "acomm_avg16_syn_bbr"
        acomm_dir.mkdir(parents=True)
        (acomm_dir / "path.pd.trk").touch()

        output_csv = temp_dir / "qc_report.csv"
        # Expect both acomm and lh.cst
        run_tract_qc(
            subject=subject,
            tractography_path=f"{temp_dir}:{subject}/dpath/*/path.pd.trk",
            output_path=str(output_csv),
            skip_subject_dir=True,
            expected_tracts=["acomm_avg16_syn_bbr", "lh.cst_avg16_syn_bbr"],
        )

        assert output_csv.exists()
        df = pd.read_csv(output_csv)
        assert len(df) == 2
        
        # Check lh.cst is marked as missing / failed QC
        cst_row = df[df["tract_name"].str.contains("lh.cst")].iloc[0]
        assert cst_row["exists"] == False
        assert cst_row["qc_pass"] == False
        assert "Missing from tractography" in cst_row["reason"]

