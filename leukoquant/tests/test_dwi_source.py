"""Which TRACULA run supplies the DWI to DTI, NODDI and the metrics step, and the line that says so.

Found 2026-11-21 (ADNI-3): the DWI comes from the first parcellation of the run's list, FreeSurfer if it is listed. A
single-parcellation GIF run therefore reads tracula-gif's dwi.nii.gz while the DTI/NODDI maps on disk were fitted on
tracula-freesurfer's. The two files hold the same data up to rounding (about 6e-6 of the b0 mean), and the metrics step
amplifies that into differences of a few tenths of a percent, with some larger. The rule is kept; the run now states it."""
import re
from pathlib import Path

import pytest

from leukoquant.utils.dwi_source import describe_dwi_source, primary_parcellation

SMK = Path(__file__).resolve().parents[1] / "workflow" / "workflows" / "process_all_workflow.smk"


@pytest.mark.parametrize("parcellations, expected", [
    (["freesurfer"], "freesurfer"),
    (["gif", "freesurfer"], "freesurfer"),      # freesurfer wins wherever it stands in the list
    (["freesurfer", "gif"], "freesurfer"),
    (["gif"], "gif"),
    (["custom", "gif"], "custom"),              # otherwise the first one
])
def test_primary_parcellation_rule(parcellations, expected):
    assert primary_parcellation(parcellations) == expected


def test_describe_names_the_tracula_run_the_dwi_comes_from():
    msg = describe_dwi_source(["gif", "freesurfer"])
    assert msg.startswith("[dwi] ")
    assert "tracula-freesurfer" in msg
    assert "WARNING" not in msg


def test_describe_warns_when_the_dwi_is_not_the_freesurfer_one():
    msg = describe_dwi_source(["gif"])
    assert "tracula-gif" in msg
    assert "WARNING" in msg and "freesurfer" in msg      # tells how to get the full-run behaviour


def test_workflow_uses_the_helper_and_prints_the_line():
    text = SMK.read_text()
    assert "primary_parcellation(PARCELLATIONS)" in text
    assert "describe_dwi_source(PARCELLATIONS)" in text
    # The helper is imported only after the configured checkout is on sys.path.
    assert text.index("sys.path.insert(0, leukoquant_parent_dir)") < text.index("from leukoquant.utils.dwi_source import")
