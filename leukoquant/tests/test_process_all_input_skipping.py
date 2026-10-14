"""Tests for process-all input failure marker skipping and --retry-input-failures."""

import os
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from leukoquant.core.process_all_processor import ProcessAllProcessor
from leukoquant.cli.main import main
from leukoquant.utils.snakemake_utils import write_run_config
from leukoquant.utils.dwi_utils import check_input_failure_marker


@pytest.fixture
def test_env(tmp_path):
    out_dir = tmp_path / "output"
    out_dir.mkdir()

    subj1 = "subj-001"
    subj2 = "subj-002"

    subj_file = tmp_path / "subjects.txt"
    subj_file.write_text(f"{subj1}\n{subj2}\n")

    # Mock inputs
    data_dir = tmp_path / "data"
    data_dir.mkdir()

    for s in (subj1, subj2):
        s_dir = data_dir / s
        s_dir.mkdir()
        (s_dir / "t1.nii.gz").write_text("t1")
        (s_dir / "flair.nii.gz").write_text("flair")
        (s_dir / "dwi.nii.gz").write_text("dwi")
        (s_dir / "dwi.bvec").write_text("0 0 0")
        (s_dir / "dwi.bval").write_text("0 1000 1000")

    return {
        "out_dir": out_dir,
        "subj_file": subj_file,
        "data_dir": data_dir,
        "subj1": subj1,
        "subj2": subj2,
    }


def test_skips_subject_with_dwi_failed_marker(test_env):
    out_dir = test_env["out_dir"]
    subj1 = test_env["subj1"]
    data_dir = test_env["data_dir"]

    # Mark subj1 as failed
    checks_dir = out_dir / subj1 / "input_checks"
    checks_dir.mkdir(parents=True)
    marker = checks_dir / "dwi_failed.txt"
    marker.write_text("Insufficient directions: 3 < 6")

    proc = ProcessAllProcessor()

    with patch("leukoquant.core.process_all_processor.SnakemakeApi") as mock_api, \
         patch("leukoquant.core.process_all_processor._resolve_fs_license", return_value=Path("/fake/lic")), \
         patch("leukoquant.core.process_all_processor._resolve_gif_home", return_value=Path("/fake/gif")), \
         patch("leukoquant.core.process_all_processor.write_run_config", wraps=write_run_config) as spy_write_config:

        mock_api_inst = MagicMock()
        mock_api.return_value.__enter__.return_value = mock_api_inst
        mock_wf = MagicMock()
        mock_api_inst.workflow.return_value = mock_wf
        mock_dag = MagicMock()
        mock_wf.dag.return_value = mock_dag
        mock_dag.execute_workflow.return_value = True

        proc.run_process_all(
            subject_input=str(test_env["subj_file"]),
            t1_pattern=str(data_dir / "{subject}" / "t1.nii.gz"),
            flair_pattern=str(data_dir / "{subject}" / "flair.nii.gz"),
            dwi_pattern=str(data_dir / "{subject}" / "dwi.nii.gz"),
            bvecs_pattern="",
            bvals_pattern="",
            mask_pattern=None,
            output_dir=str(out_dir),
            retry_input_failures=False,
        )

        call_args = spy_write_config.call_args
        config_passed = call_args[0][1]
        assert config_passed["subjects"] == ["subj-002"]
        assert marker.exists()


def test_retries_subject_when_flag_set(test_env):
    out_dir = test_env["out_dir"]
    subj1 = test_env["subj1"]
    data_dir = test_env["data_dir"]

    checks_dir = out_dir / subj1 / "input_checks"
    checks_dir.mkdir(parents=True)
    marker = checks_dir / "dwi_failed.txt"
    marker.write_text("Insufficient directions: 3 < 6")

    proc = ProcessAllProcessor()

    with patch("leukoquant.core.process_all_processor.SnakemakeApi") as mock_api, \
         patch("leukoquant.core.process_all_processor._resolve_fs_license", return_value=Path("/fake/lic")), \
         patch("leukoquant.core.process_all_processor._resolve_gif_home", return_value=Path("/fake/gif")), \
         patch("leukoquant.core.process_all_processor.write_run_config", wraps=write_run_config) as spy_write_config:

        mock_api_inst = MagicMock()
        mock_api.return_value.__enter__.return_value = mock_api_inst
        mock_wf = MagicMock()
        mock_api_inst.workflow.return_value = mock_wf
        mock_dag = MagicMock()
        mock_wf.dag.return_value = mock_dag
        mock_dag.execute_workflow.return_value = True

        proc.run_process_all(
            subject_input=str(test_env["subj_file"]),
            t1_pattern=str(data_dir / "{subject}" / "t1.nii.gz"),
            flair_pattern=str(data_dir / "{subject}" / "flair.nii.gz"),
            dwi_pattern=str(data_dir / "{subject}" / "dwi.nii.gz"),
            bvecs_pattern="",
            bvals_pattern="",
            mask_pattern=None,
            output_dir=str(out_dir),
            retry_input_failures=True,
        )

        call_args = spy_write_config.call_args
        config_passed = call_args[0][1]
        assert "subj-001" in config_passed["subjects"]
        assert "subj-002" in config_passed["subjects"]
        assert not marker.exists()


def test_newer_input_file_clears_marker_and_reprocesses(test_env):
    out_dir = test_env["out_dir"]
    subj1 = test_env["subj1"]
    data_dir = test_env["data_dir"]

    checks_dir = out_dir / subj1 / "input_checks"
    checks_dir.mkdir(parents=True)
    marker = checks_dir / "dwi_failed.txt"
    marker.write_text("Insufficient directions")

    # Set marker mtime in the past
    past = time.time() - 100
    os.utime(marker, (past, past))

    # DWI file is newer (now)
    dwi_file = data_dir / subj1 / "dwi.nii.gz"
    now = time.time()
    os.utime(dwi_file, (now, now))

    proc = ProcessAllProcessor()

    with patch("leukoquant.core.process_all_processor.SnakemakeApi") as mock_api, \
         patch("leukoquant.core.process_all_processor._resolve_fs_license", return_value=Path("/fake/lic")), \
         patch("leukoquant.core.process_all_processor._resolve_gif_home", return_value=Path("/fake/gif")), \
         patch("leukoquant.core.process_all_processor.write_run_config", wraps=write_run_config) as spy_write_config:

        mock_api_inst = MagicMock()
        mock_api.return_value.__enter__.return_value = mock_api_inst
        mock_wf = MagicMock()
        mock_api_inst.workflow.return_value = mock_wf
        mock_dag = MagicMock()
        mock_wf.dag.return_value = mock_dag
        mock_dag.execute_workflow.return_value = True

        proc.run_process_all(
            subject_input=str(test_env["subj_file"]),
            t1_pattern=str(data_dir / "{subject}" / "t1.nii.gz"),
            flair_pattern=str(data_dir / "{subject}" / "flair.nii.gz"),
            dwi_pattern=str(data_dir / "{subject}" / "dwi.nii.gz"),
            bvecs_pattern="",
            bvals_pattern="",
            mask_pattern=None,
            output_dir=str(out_dir),
            retry_input_failures=False,
        )

        call_args = spy_write_config.call_args
        config_passed = call_args[0][1]
        assert "subj-001" in config_passed["subjects"]
        assert not marker.exists()


def test_cli_accepts_retry_input_failures_flag():
    runner = CliRunner()
    result = runner.invoke(main, ["process-all", "--help"])
    assert result.exit_code == 0
    assert "--retry-input-failures" in result.output


def test_check_input_failure_marker_helper(tmp_path):
    session_dir = tmp_path / "sub-01"
    checks_dir = session_dir / "input_checks"
    checks_dir.mkdir(parents=True)
    marker = checks_dir / "dwi_failed.txt"
    marker.write_text("Error: insufficient gradient directions")

    raw_file = tmp_path / "dwi.nii.gz"
    raw_file.write_text("dummy")

    # When marker is newer than raw file: should skip
    past = time.time() - 100
    os.utime(raw_file, (past, past))
    assert check_input_failure_marker(session_dir, raw_input_paths=[raw_file], retry=False) is True
    assert marker.exists()

    # When retry=True: should not skip and should unlink marker
    assert check_input_failure_marker(session_dir, raw_input_paths=[raw_file], retry=True) is False
    assert not marker.exists()

    # When marker doesn't exist: should not skip
    assert check_input_failure_marker(session_dir, raw_input_paths=[raw_file], retry=False) is False

    # When raw file is newer than marker: should not skip and should unlink marker
    marker.write_text("Error again")
    os.utime(marker, (past, past))
    now = time.time()
    os.utime(raw_file, (now, now))
    assert check_input_failure_marker(session_dir, raw_input_paths=[raw_file], retry=False) is False
    assert not marker.exists()

