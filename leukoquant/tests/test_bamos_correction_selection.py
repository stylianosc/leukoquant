"""Which lesion corrections each -corr option selects (BaMoS correction script).

Every correction is independent and selected purely by its -corr option names."""
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "external" / "bamos" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import correction_lesions_101126 as cl  # noqa: E402


def selected(requested):
    return [c.log_name for c in cl.CORRECTIONS if c.is_selected(set(requested))]


@pytest.mark.parametrize("requested, expected", [
    (["all"], ["choroid", "cortical", "sheet"]),
    (["choroid"], ["choroid"]),
    (["choroid", "cortex", "sheet"], ["choroid", "cortical", "sheet"]),
    (["choroid-cortex"], ["choroid", "cortical"]),
    (["choroid-sheet"], ["choroid", "sheet"]),
    (["cortical-sheet"], ["cortical", "sheet"]),
    (["cortex"], ["cortical"]),
    (["sheet"], ["sheet"]),
    (["hypot1"], ["hypo corr"]),
    (["ventr"], ["ventr"]),
    (["hypot1", "ventr", "cortex"], ["cortical", "hypo corr", "ventr"]),
])
def test_selection(requested, expected):
    assert selected(requested) == expected


def test_all_does_not_select_hypot1_or_ventr():
    assert "hypo corr" not in selected(["all"])
    assert "ventr" not in selected(["all"])


def test_only_choroid_options_select_choroid():
    for opt in ("cortex", "sheet", "hypot1", "ventr", "cortical-sheet"):
        assert "choroid" not in selected([opt])
