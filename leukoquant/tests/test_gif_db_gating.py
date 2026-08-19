"""Tests for GIF-atlas-database auto-download and its per-modality gating.

Covers two things:
  1. `_resolve_gif_home(ensure_t1_db=...)` actually honors the flag (skips
     the T1/db_mideface check-and-download entirely when False).
  2. `run_gif()` derives that flag (and whether to check/download the FLAIR
     db) correctly from the actual per-subject T1/FLAIR files in a batch,
     matching gif_workflow.smk's per-subject db selection: a subject uses
     the T1 db whenever it has a T1 file, and only falls back to the FLAIR
     db when it doesn't. A batch can therefore need T1 only, FLAIR only,
     both (mixed batch), or -- when both modalities are given for every
     subject -- T1 only, since T1 always wins per-subject.

Network calls (`ensure_gif_db`) and the actual Snakemake/container
invocation are mocked out, so these tests are fast, hermetic, and need no
HF_TOKEN, network access, or real GIF/container install.
"""
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from leukoquant.core.gif_processor import GIFProcessor
from leukoquant.utils.container_utils import GIF_DB_FLAIR_FILENAME
from leukoquant.utils.external_utils import _resolve_gif_home


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def gif_processor(temp_dir):
    """GIFProcessor pointed at a throwaway external/gif/ dir with just enough
    structure to pass __init__'s validation (a stub GIF_111125.sh) -- never a
    real GIF install."""
    external_dir = temp_dir / "external"
    gif_dir = external_dir / "gif"
    gif_dir.mkdir(parents=True)
    (gif_dir / "GIF_111125.sh").write_text("#!/bin/bash\necho stub\n")
    return GIFProcessor(external_dir=str(external_dir))


def _make_pattern(temp_dir, modality, present_for):
    """Touch empty placeholder files for `present_for` subjects and return
    the {subject}-glob pattern for them."""
    d = temp_dir / f"{modality}_input"
    d.mkdir(exist_ok=True)
    for s in present_for:
        (d / f"{s}_{modality}.nii.gz").touch()
    return str(d / ("{subject}_" + modality + ".nii.gz"))


def _write_subjects_file(temp_dir, subjects):
    p = temp_dir / "subjects.txt"
    p.write_text("\n".join(subjects) + "\n")
    return str(p)


class TestResolveGifHomeGating:
    """_resolve_gif_home()'s ensure_t1_db flag in isolation."""

    def test_ensure_t1_db_false_skips_check_and_download(self, temp_dir):
        gif_dir = temp_dir / "gif_install"
        gif_dir.mkdir()
        with patch("leukoquant.utils.container_utils.ensure_gif_db") as mock_ensure:
            result = _resolve_gif_home(gif_dir, ensure_t1_db=False)
        mock_ensure.assert_not_called()
        assert result == gif_dir / "GIF"
        assert not (result / "db_mideface").exists()

    def test_ensure_t1_db_true_checks_and_downloads(self, temp_dir):
        gif_dir = temp_dir / "gif_install"
        gif_dir.mkdir()

        def _fake_download(db_dir, filename):
            # Simulate a successful download: materialize db.xml.
            Path(db_dir).mkdir(parents=True, exist_ok=True)
            (Path(db_dir) / "db.xml").write_text("<document/>")

        with patch(
            "leukoquant.utils.container_utils.ensure_gif_db", side_effect=_fake_download
        ) as mock_ensure:
            result = _resolve_gif_home(gif_dir, ensure_t1_db=True)
        mock_ensure.assert_called_once()
        assert (result / "db_mideface" / "db.xml").exists()

    def test_ensure_t1_db_true_raises_if_download_fails(self, temp_dir):
        gif_dir = temp_dir / "gif_install"
        gif_dir.mkdir()
        with patch("leukoquant.utils.container_utils.ensure_gif_db"):  # no-op, leaves db.xml missing
            with pytest.raises(FileNotFoundError, match="db_mideface/db.xml"):
                _resolve_gif_home(gif_dir, ensure_t1_db=True)


class TestRunGifModalityGating:
    """run_gif()'s end-to-end decision of which db(s) to require/download."""

    @pytest.mark.parametrize(
        "case,t1_present,flair_present,expect_ensure_t1_db,expect_flair_download",
        [
            ("t1_only",    ["sub-01", "sub-02"], [],                    True,  False),
            ("flair_only", [],                    ["sub-01", "sub-02"], False, True),
            ("mixed",      ["sub-01"],            ["sub-02"],            True,  True),
            ("both_given", ["sub-01", "sub-02"], ["sub-01", "sub-02"], True,  False),
        ],
    )
    def test_gating(
        self,
        temp_dir,
        gif_processor,
        case,
        t1_present,
        flair_present,
        expect_ensure_t1_db,
        expect_flair_download,
    ):
        subjects = ["sub-01", "sub-02"]
        subjects_file = _write_subjects_file(temp_dir, subjects)
        t1_pattern = _make_pattern(temp_dir, "t1", t1_present)
        flair_pattern = _make_pattern(temp_dir, "flair", flair_present)
        out_dir = temp_dir / "out"

        fake_gif_home = temp_dir / "GIF"
        fake_gif_home.mkdir(exist_ok=True)

        def _fake_download(db_dir, filename):
            # Simulate a successful download: materialize db.xml, same as a
            # real ensure_gif_db() call would after extracting the tarball.
            Path(db_dir).mkdir(parents=True, exist_ok=True)
            (Path(db_dir) / "db.xml").write_text("<document/>")

        with patch(
            "leukoquant.core.gif_processor._resolve_gif_home", return_value=fake_gif_home
        ) as mock_resolve, patch(
            "leukoquant.utils.container_utils.ensure_gif_db", side_effect=_fake_download
        ) as mock_ensure_flair_db, patch(
            "leukoquant.core.gif_processor.ensure_container"
        ), patch(
            "leukoquant.core.gif_processor.subprocess.run",
            return_value=MagicMock(returncode=0, stdout="", stderr=""),
        ):
            gif_processor.run_gif(
                subject_input=subjects_file,
                output_dir=str(out_dir),
                t1_pattern=t1_pattern,
                flair_pattern=flair_pattern,
            )

            # _resolve_gif_home must be told whether the batch needs the T1 db.
            assert mock_resolve.call_args.kwargs["ensure_t1_db"] == expect_ensure_t1_db, (
                f"[{case}] expected ensure_t1_db={expect_ensure_t1_db}, "
                f"got call: {mock_resolve.call_args}"
            )

            # ensure_gif_db (FLAIR path) must only be invoked when at least
            # one subject in the batch actually falls back to the FLAIR db.
            flair_calls = [
                c for c in mock_ensure_flair_db.call_args_list
                if c.args[1:2] == (GIF_DB_FLAIR_FILENAME,)
                or c.kwargs.get("filename") == GIF_DB_FLAIR_FILENAME
            ]
            assert bool(flair_calls) == expect_flair_download, (
                f"[{case}] expected FLAIR db download={expect_flair_download}, "
                f"got calls: {mock_ensure_flair_db.call_args_list}"
            )

    def test_pure_t1_batch_never_requires_flair_db_directory_to_exist(
        self, temp_dir, gif_processor
    ):
        """A T1-only batch must succeed even when db_FLAIR_mideface doesn't
        exist at all under GIF_HOME -- the FileNotFoundError guard for a
        missing FLAIR db must be skipped entirely, not just its download."""
        subjects = ["sub-01"]
        subjects_file = _write_subjects_file(temp_dir, subjects)
        t1_pattern = _make_pattern(temp_dir, "t1", subjects)
        out_dir = temp_dir / "out"

        fake_gif_home = temp_dir / "GIF"
        fake_gif_home.mkdir(exist_ok=True)
        # Deliberately do NOT create db_FLAIR_mideface/ under fake_gif_home.

        with patch(
            "leukoquant.core.gif_processor._resolve_gif_home", return_value=fake_gif_home
        ), patch("leukoquant.utils.container_utils.ensure_gif_db") as mock_ensure, patch(
            "leukoquant.core.gif_processor.ensure_container"
        ), patch(
            "leukoquant.core.gif_processor.subprocess.run",
            return_value=MagicMock(returncode=0, stdout="", stderr=""),
        ):
            # Should not raise, despite db_FLAIR_mideface being entirely absent.
            gif_processor.run_gif(
                subject_input=subjects_file,
                output_dir=str(out_dir),
                t1_pattern=t1_pattern,
                flair_pattern=None,
            )
        mock_ensure.assert_not_called()
