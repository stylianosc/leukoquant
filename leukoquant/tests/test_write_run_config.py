"""write_run_config: each run's Snakemake config is its own file.

Jobs re-read their config by path when they start, possibly long after
submission. If every run wrote the same fixed file, a later run would change an
earlier run's still-queued jobs (EPAD, 2026-09-27). Each run must therefore get
a config file no later run touches, while the fixed-name file keeps recording
the latest run's settings.
"""
import re
from pathlib import Path

import yaml

from leukoquant.utils.snakemake_utils import RUN_CONFIG_DIR, write_run_config


def test_each_run_gets_its_own_file_and_earlier_runs_stay_intact(tmp_path):
    record = tmp_path / "process_all_config.yaml"
    first = write_run_config(record, {"subjects": ["a", "b", "c"]})
    second = write_run_config(record, {"subjects": ["b"]})

    assert first != second
    for p in (first, second):
        assert Path(p).parent == tmp_path / RUN_CONFIG_DIR
        assert re.fullmatch(r"process_all_config_\d{8}-\d{6}_[0-9a-f]{8}\.yaml", Path(p).name)
    # The earlier run's config is unchanged by the later run ...
    assert yaml.safe_load(open(first)) == {"subjects": ["a", "b", "c"]}
    assert yaml.safe_load(open(second)) == {"subjects": ["b"]}
    # ... while the fixed-name record holds the latest run's settings.
    assert yaml.safe_load(open(record)) == {"subjects": ["b"]}


def test_no_temp_files_left_behind(tmp_path):
    record = tmp_path / "z_score_config.yaml"
    write_run_config(str(record), {"covariates": ["sex_binary", "age"]})
    leftovers = [p.name for p in tmp_path.iterdir() if p.name.startswith(".") and p.is_file()]
    assert leftovers == []
