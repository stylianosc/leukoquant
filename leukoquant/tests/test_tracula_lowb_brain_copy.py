"""The skull-stripped b0 (dmri/lowb_brain.nii.gz) is copied to the outputs directory but is not a declared rule output.

Declaring it would add a dependency: every subject processed before the file was kept would then look incomplete and be
rescheduled. Copying it from the shell block only keeps it for new runs and changes nothing for finished subjects."""
import re
from pathlib import Path

SMK = Path(__file__).resolve().parents[1] / "workflow" / "workflows" / "tracula_workflow.smk"


def _copy_list():
    text = SMK.read_text()
    block = re.search(r"folder_paths=\((.*?)\n\s*\)", text, re.S)
    assert block, "folder_paths list not found"
    return re.findall(r'"([^"]+)"', block.group(1))


def test_lowb_brain_is_copied_to_outputs():
    assert "dmri/lowb_brain.nii.gz" in _copy_list()


def test_lowb_brain_is_not_a_declared_rule_output():
    # Only the copy list may name it; no input/output/target line may require it.
    lines = [l for l in SMK.read_text().splitlines() if "lowb_brain" in l]
    assert lines == ['                "dmri/lowb_brain.nii.gz"'], lines
