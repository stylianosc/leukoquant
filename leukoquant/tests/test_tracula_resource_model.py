"""The TRACULA memory and scratch requests are fitted lines of the DWI size, not size classes.

x = megavoxel-volumes (voxels x volumes / 1e6). Each quantity y is fitted as y = intercept + slope * x by least squares;
the margin is an upper quantile of the residuals (so the request covers that share of the measured runs), not a factor."""
import csv
import json
import math

import numpy as np
import pytest

from leukoquant.utils import tracula_resource_model as rm


def test_fit_line_recovers_a_known_line_and_margin_is_the_residual_quantile():
    rng = np.random.default_rng(0)
    x = np.linspace(10, 300, 400)
    noise = rng.normal(0, 0.3, x.size)
    y = 1.5 + 0.02 * x + noise
    fit = rm.fit_line(x, y, quantile=0.99)
    assert fit["intercept"] == pytest.approx(1.5, abs=0.1)
    assert fit["slope"] == pytest.approx(0.02, abs=0.001)
    resid = y - (fit["intercept"] + fit["slope"] * x)
    assert fit["margin"] == pytest.approx(np.quantile(resid, 0.99))
    assert fit["n"] == 400
    assert 0.0 < fit["r2"] <= 1.0
    assert fit["quantile"] == 0.99


def test_fit_line_margin_covers_the_chosen_share_of_the_data():
    rng = np.random.default_rng(1)
    x = rng.uniform(10, 300, 2000)
    y = 1.0 + 0.01 * x + rng.exponential(0.2, x.size)
    fit = rm.fit_line(x, y, quantile=0.99)
    covered = np.mean(y <= fit["intercept"] + fit["slope"] * x + fit["margin"])
    assert covered >= 0.99


def test_fit_needs_enough_points():
    with pytest.raises(ValueError):
        rm.fit_line([1.0, 2.0], [1.0, 2.0], quantile=0.99)


def _row(**kw):
    base = dict(dataset="X", parc="gif", dim_x=100, dim_y=100, dim_z=50, n_vol=40, prep_rss_mb=2000, bedp_rss_mb=2000,
                path_rss_mb=2000, scratch_gb_at_cleanup=1.0, completed=1)
    base.update(kw)
    return base


def test_rows_are_selected_for_each_quantity():
    rows = [_row(),                                                    # usable for all three
            _row(completed=0, path_rss_mb="", scratch_gb_at_cleanup=2.0),   # failed run: prep RSS valid, path/scratch not
            _row(n_vol=72),                                            # header of an older run than the log: dropped
            _row(dim_x=""),                                            # no DWI size: dropped
            _row(prep_rss_mb="", bedp_rss_mb="", path_rss_mb="")]                      # no memory lines (log before 2026-09-30)
    sel = rm.select_rows(rows)
    assert len(sel["scratch_gb"]) == 2      # the complete run and the one without memory lines
    assert len(sel["prep_rss_gb"]) == 2     # the complete run and the failed one (prep finished before it failed)
    assert len(sel["path_rss_gb"]) == 1     # only the complete run
    assert sel["path_rss_gb"][0]["mvv"] == pytest.approx(100 * 100 * 50 * 40 / 1e6)
    # a failed run must not feed the path or scratch fits
    assert all(r["completed"] == 1 for r in sel["path_rss_gb"])
    assert all(r["completed"] == 1 for r in sel["scratch_gb"])


def test_build_model_has_the_three_fits_and_provenance(tmp_path):
    rng = np.random.default_rng(2)
    rows = []
    for _ in range(60):
        nvol = int(rng.integers(20, 130))
        dim = int(rng.choice([96, 112, 128, 256]))
        mvv = dim * dim * 80 * nvol / 1e6
        rows.append(_row(dim_x=dim, dim_y=dim, dim_z=80, n_vol=nvol, prep_rss_mb=(1 + 0.005 * mvv) * 1024,
                         path_rss_mb=(1.2 + 0.01 * mvv) * 1024, scratch_gb_at_cleanup=0.5 + 0.011 * mvv))
    model = rm.build_model(rows, quantile=0.99)
    for key in ("scratch_gb", "prep_rss_gb", "path_rss_gb"):
        for field in ("intercept", "slope", "margin", "n", "r2", "quantile"):
            assert field in model[key], (key, field)
    assert model["scratch_gb"]["slope"] == pytest.approx(0.011, rel=0.02)
    assert model["path_rss_gb"]["slope"] == pytest.approx(0.01, rel=0.02)
    assert "limits" in model and "rows" in model


def test_committed_model_is_what_the_committed_data_gives():
    """The JSON the pipeline reads must be reproducible from the committed fit input with the documented script."""
    with open(rm.FIT_INPUT_PATH, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    fresh = rm.build_model(rows, quantile=rm.MARGIN_QUANTILE)
    committed = json.load(open(rm.MODEL_PATH))
    for key in ("scratch_gb", "prep_rss_gb", "path_rss_gb"):
        for field in ("intercept", "slope", "margin"):
            assert committed[key][field] == pytest.approx(fresh[key][field], rel=1e-6), (key, field)


def test_committed_model_covers_the_measured_extremes():
    """Sanity of the fit against the data it was fitted on: scratch covers the large ADNI-3 runs (4.35 and 4.88 GB)."""
    model = rm.load_model()
    s = model["scratch_gb"]
    for vols, used in ((36, 4.35), (54, 4.88)):
        mvv = 256 * 256 * 80 * vols / 1e6
        assert s["intercept"] + s["slope"] * mvv + s["margin"] >= used - 0.1


def test_predict_and_request_are_clamped_only_at_the_ends():
    fit = {"intercept": 1.0, "slope": 0.01, "margin": 0.5}
    assert rm.predict_gb(fit, 100.0) == pytest.approx(2.5)
    assert rm.request_gb(fit, 100.0, low=8, high=16) == 8                 # floor binds
    assert rm.request_gb(fit, 900.0, low=8, high=16) == math.ceil(10.5)   # inside: ceil of the prediction
    assert rm.request_gb(fit, 5000.0, low=8, high=16) == 16               # cap binds
