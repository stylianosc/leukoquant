"""The scratch folder of a metrics task must be unique per task.

The GIF and FreeSurfer metrics tasks of one session start together on the same node and share /scratch0. A folder
named from the session alone was shared, and the task that finished first deleted it under the other one.
"""
import re
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "utils" / "metrics_calc.sh"


def _scratch_name(subject, output_dir, job="1", task="2"):
    text = SCRIPT.read_text()
    start = text.index("_scratch_subject=")
    end = text.index('mkdir -p "$SCRATCH_DIR"')
    snippet = text[start:end]
    script = f'SUBJECT="{subject}"\nOUTPUT_DIR="{output_dir}"\n_scratch_base=/scratch0\n{snippet}\necho "$SCRATCH_DIR"\n'
    env = {"JOB_ID": job, "SGE_TASK_ID": task, "PATH": "/usr/bin:/bin"}
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env, check=True).stdout.strip()


def test_gif_and_freesurfer_tasks_of_one_session_get_different_scratch_folders():
    gif = _scratch_name("subj-021-s-6890", "/output/subj-021-s-6890/metrics-gif/outputs/metrics")
    fs = _scratch_name("subj-021-s-6890", "/output/subj-021-s-6890/metrics-freesurfer/outputs/metrics")
    assert gif != fs
    assert "metrics-gif" in gif and "metrics-freesurfer" in fs


def test_session_ids_with_a_slash_do_not_create_nested_folders():
    name = _scratch_name("sub-OAS30573/ses-d0515", "/output/x/metrics-freesurfer/outputs/metrics")
    assert re.fullmatch(r"/scratch0/[^/]+", name), name


def test_two_runs_of_the_same_task_still_differ():
    a = _scratch_name("s", "/o/metrics-gif/outputs", job="7", task="3")
    b = _scratch_name("s", "/o/metrics-gif/outputs", job="7", task="3")
    assert a != b  # process id and $RANDOM
