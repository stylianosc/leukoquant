"""usage_report: merge per-job usage tables, flag runs near their request, refit the request model from TRUE peaks."""
from leukoquant.utils import usage_record as ur
from leukoquant.utils import usage_report as rep


def _rec(**kw):
    base = {"dataset": "D", "subject": "s", "session": "", "rule": "trac_all_freesurfer", "stage": "trac-all path",
            "job_id": "1", "task_id": "1", "dwi_mvv": "100", "req_mem_mb": "8192", "req_scratch_mb": "6144",
            "peak_rss_mb": "2000", "peak_scratch_mb": "3000", "killed": "0", "timestamp_utc": "t"}
    base.update(kw)
    return base


def test_near_miss_flags_runs_above_80_percent_or_killed():
    rows = [_rec(), _rec(subject="hot", peak_scratch_mb="5500"), _rec(subject="dead", killed="1"),
            _rec(subject="mem", peak_rss_mb="7000")]
    flagged = {r["subject"]: why for r, why in rep.near_misses(rows)}
    assert set(flagged) == {"hot", "dead", "mem"}
    assert "scratch" in flagged["hot"] and "killed" in flagged["dead"] and "memory" in flagged["mem"]


def test_refit_uses_true_peaks_and_reports_fit_quality(tmp_path):
    rows = [_rec(subject=f"s{i}", dwi_mvv=str(50 + 10 * i), peak_scratch_mb=str(1000 * (1 + i)),
                 peak_rss_mb=str(500 * (2 + i))) for i in range(6)]
    fit = rep.refit(rows)
    assert fit["scratch_gb"]["n"] == 6 and fit["scratch_gb"]["r2"] > 0.99
    assert fit["memory_gb"]["slope"] > 0


def test_refit_takes_the_largest_peak_over_the_stages_of_one_job():
    rows = [_rec(stage="a", peak_scratch_mb="1000"), _rec(stage="b", peak_scratch_mb="4000")]
    assert rep.per_job_peaks(rows)[0]["peak_scratch_mb"] == 4000


def test_main_collects_then_reports(tmp_path, capsys):
    src = tmp_path / "logs" / "usage_records.tsv"
    src.parent.mkdir()
    for i in range(4):
        ur.append_record(src, _rec(subject=f"s{i}", dwi_mvv=str(50 + 10 * i), peak_scratch_mb=str(1000 + 500 * i),
                                   peak_rss_mb=str(2000 + 300 * i), timestamp_utc=f"t{i}"))
    central = tmp_path / "usage_table.tsv"
    assert rep.main(["--central", str(central), "--source", str(src)]) == 0
    out = capsys.readouterr().out
    assert "4 new" in out and "scratch" in out and central.exists()
