"""Every Snakemake workflow must put the run's own LeukoQuant tree first on sys.path before it imports leukoquant.

Background (EPAD and OASIS3 pinned metrics runs, 2026-11-19): all 4,044 array tasks failed at start with
"cannot import name 'tracula_scratch_mb' from 'leukoquant.utils.mem_override'". The cluster tasks are started
without PYTHONPATH and the virtual environment's editable install points at the production checkout. The tasks
were meant to run a different checkout (leukoquant_parent_dir of the run config). process_all_workflow.smk
imported leukoquant.utils.z_score_utils (line ~482, only when z-score covariates are on) before any workflow
had inserted that checkout into sys.path, so the production package was loaded and cached, and the later
sys.path.insert in tracula_workflow.smk could not change it. The rule is therefore: in every .smk file the first
import of the leukoquant package comes after the sys.path.insert of the configured checkout.
"""
import re
from pathlib import Path

import pytest

WORKFLOWS = sorted((Path(__file__).resolve().parents[1] / "workflow" / "workflows").glob("*.smk"))

# Any import of the package, at any indentation (the failing one sat inside an ``if`` block).
_IMPORT = re.compile(r"^\s*(?:from\s+leukoquant[\w.]*\s+import|import\s+leukoquant\b)", re.M)
# The insertion of the configured checkout (parent dir, or the parent of leukoquant_dir) at the front of sys.path.
_INSERT = re.compile(r"^\s*sys\.path\.insert\(\s*0\s*,\s*[^)]*(?:PARENT_DIR|parent_dir|dirname\(\s*LEUKOQUANT_DIR\s*\))", re.M | re.I)


def test_workflows_found():
    assert len(WORKFLOWS) >= 10


@pytest.mark.parametrize("path", WORKFLOWS, ids=lambda p: p.name)
def test_sys_path_is_set_before_the_first_leukoquant_import(path):
    text = path.read_text()
    imp = _IMPORT.search(text)
    if imp is None:
        return  # nothing imported, nothing to order
    ins = _INSERT.search(text)
    assert ins is not None, f"{path.name} imports leukoquant but never puts the configured checkout on sys.path"
    assert ins.start() < imp.start(), (
        f"{path.name}: leukoquant is imported at line {text[:imp.start()].count(chr(10)) + 1} before sys.path is set "
        f"at line {text[:ins.start()].count(chr(10)) + 1}; an editable install elsewhere would win"
    )
