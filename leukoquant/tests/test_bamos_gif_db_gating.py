"""Tests for BaMoS's gated GIF-atlas-database auto-download.

BaMoS only needs the T1 GIF database (db_mideface) when at least one subject
in the batch lacks externally-supplied --gif-results-dir output and must
therefore fall back to running GIF internally (BaMoS always requires T1, so
that internal fallback never touches the FLAIR db). These tests verify the
download is attempted exactly when that condition holds, and never
otherwise -- e.g. a batch where every subject supplies external GIF results
must never trigger a ~900MB download just because a GIF install happens to
be configured.

Network calls (`ensure_gif_db`) and the Snakemake/container invocation are
mocked out, so these tests are fast, hermetic, and need no HF_TOKEN, network
access, or real GIF/container install.
"""
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from leukoquant.core.bamos_processor import BaMoSProcessor


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


def _make_files(temp_dir, modality, subjects):
    d = temp_dir / f"{modality}_input"
    d.mkdir(exist_ok=True)
    for s in subjects:
        (d / f"{s}_{modality}.nii.gz").touch()
    return str(d / ("{subject}_" + modality + ".nii.gz"))


def _make_gif_results_dirs(temp_dir, subjects_with_results):
    d = temp_dir / "gif_results"
    d.mkdir(exist_ok=True)
    for s in subjects_with_results:
        (d / s).mkdir(exist_ok=True)
    return str(d / "{subject}")


def _write_subjects_file(temp_dir, subjects):
    p = temp_dir / "subjects.txt"
    p.write_text("\n".join(subjects) + "\n")
    return str(p)


def _fake_download(db_dir, filename):
    Path(db_dir).mkdir(parents=True, exist_ok=True)
    (Path(db_dir) / "db.xml").write_text("<document/>")


class TestBamosGifDbGating:
    def _run(self, temp_dir, subjects, subjects_with_external_results, fake_gif_home):
        subjects_file = _write_subjects_file(temp_dir, subjects)
        t1_pattern = _make_files(temp_dir, "t1", subjects)
        flair_pattern = _make_files(temp_dir, "flair", subjects)
        gif_results_pattern = (
            _make_gif_results_dirs(temp_dir, subjects_with_external_results)
        )
        out_dir = temp_dir / "out"

        processor = BaMoSProcessor()

        with patch(
            "leukoquant.core.bamos_processor._find_gif_home_dir",
            return_value=fake_gif_home,
        ), patch(
            "leukoquant.utils.container_utils.ensure_gif_db",
            side_effect=_fake_download,
        ) as mock_ensure, patch(
            "leukoquant.core.bamos_processor.ensure_container"
        ), patch(
            "leukoquant.core.bamos_processor.subprocess.run",
            return_value=MagicMock(returncode=0, stdout="", stderr=""),
        ):
            processor.run_bamos(
                subject_input=subjects_file,
                flair_pattern=flair_pattern,
                t1_pattern=t1_pattern,
                output_dir=str(out_dir),
                gif_results_pattern=gif_results_pattern,
            )
        return mock_ensure

    def test_all_subjects_have_external_results_never_downloads(self, temp_dir):
        """Every subject supplies external GIF results -> BaMoS never needs
        the GIF database at all, even though a GIF install (missing its db)
        is configured."""
        subjects = ["sub-01", "sub-02"]
        fake_gif_home = temp_dir / "GIF"
        fake_gif_home.mkdir()  # exists, but db_mideface/db.xml is absent

        mock_ensure = self._run(
            temp_dir, subjects, subjects_with_external_results=subjects, fake_gif_home=fake_gif_home
        )
        mock_ensure.assert_not_called()

    def test_some_subjects_need_internal_gif_downloads_t1_db(self, temp_dir):
        """One subject lacks external results -> falls back to internal GIF
        -> T1 db must be auto-downloaded."""
        subjects = ["sub-01", "sub-02"]
        fake_gif_home = temp_dir / "GIF"
        fake_gif_home.mkdir()

        mock_ensure = self._run(
            temp_dir, subjects, subjects_with_external_results=["sub-01"], fake_gif_home=fake_gif_home
        )
        mock_ensure.assert_called_once()
        db_dir_arg = Path(mock_ensure.call_args.args[0])
        assert db_dir_arg.name == "db_mideface"
        assert (db_dir_arg / "db.xml").exists()

    def test_no_gif_home_available_skips_download_and_excludes_subjects(self, temp_dir, capsys):
        """No bundled/external GIF install at all (_find_gif_home_dir
        returns None) -> no download attempted (nowhere to place it), and
        the existing graceful skip-with-warning behavior still applies."""
        subjects = ["sub-01", "sub-02"]

        with patch(
            "leukoquant.core.bamos_processor._find_gif_home_dir", return_value=None
        ), patch(
            "leukoquant.utils.container_utils.ensure_gif_db", side_effect=_fake_download
        ) as mock_ensure:
            subjects_file = _write_subjects_file(temp_dir, subjects)
            t1_pattern = _make_files(temp_dir, "t1", subjects)
            flair_pattern = _make_files(temp_dir, "flair", subjects)
            out_dir = temp_dir / "out"
            processor = BaMoSProcessor()

            with pytest.raises(ValueError, match="No subjects remain"):
                processor.run_bamos(
                    subject_input=subjects_file,
                    flair_pattern=flair_pattern,
                    t1_pattern=t1_pattern,
                    output_dir=str(out_dir),
                    gif_results_pattern=None,
                )
        mock_ensure.assert_not_called()
