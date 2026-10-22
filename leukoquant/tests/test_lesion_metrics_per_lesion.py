"""Per-lesion metrics on synthetic lesions with known answers (merge count, coalescence band, thickness, surface, fractal, tracts)."""

import numpy as np
import pytest

from leukoquant.utils.lesion_metrics import _compute_lesion_metrics, _calculate_fractal_dimension, _mesh_surface_area

SHAPE = (60, 40, 40)


def _cube(arr, x0, size=5, y0=5, z0=5):
    arr[x0:x0 + size, y0:y0 + size, z0:z0 + size] = 1


def _sphere(arr, centre, radius):
    xx, yy, zz = np.indices(arr.shape)
    arr[(xx - centre[0]) ** 2 + (yy - centre[1]) ** 2 + (zz - centre[2]) ** 2 <= radius ** 2] = 1


def _run(lesion, tiv=1e6, tracts=None):
    tracts = tracts if tracts is not None else {"far": np.zeros(SHAPE, dtype=bool)}
    lesion_df, _ = _compute_lesion_metrics("s", lesion.astype(np.uint8), tiv, tracts, {"fa": np.full(SHAPE, 0.5)}, (1.0, 1.0, 1.0))
    return lesion_df


@pytest.fixture(scope="module")
def lesions():
    lesion = np.zeros(SHAPE)
    _cube(lesion, 2)            # cube A: x 2..6
    _cube(lesion, 8)            # cube B: x 8..12, one empty voxel layer (x=7) between A and B
    _cube(lesion, 30, y0=20, z0=20)   # cube C, far from the others
    return lesion


def _row_at(df, x_centre):
    return df.iloc[(df["centre_of_mass_x_mm"] - x_centre).abs().argmin()]


def test_merged_after_dilation_counts_distinct_lesions(lesions):
    df = _run(lesions)
    a, b, c = _row_at(df, 4), _row_at(df, 10), _row_at(df, 32)
    assert a["n_lesions_merged_after_dilation"] == 2 and bool(a["merged_after_dilation"])
    assert b["n_lesions_merged_after_dilation"] == 2 and bool(b["merged_after_dilation"])
    assert c["n_lesions_merged_after_dilation"] == 1 and not bool(c["merged_after_dilation"])     # was True for every lesion


def test_thickness_of_a_cube(lesions):
    df = _run(lesions)
    # the centre voxel of a 5-voxel cube is 3 voxels from the nearest voxel outside it
    assert _row_at(df, 32)["max_inscribed_radius_mm"] == pytest.approx(3.0)


def test_coalescence_reaches_the_upper_bands_for_a_large_lesion():
    lesion = np.zeros(SHAPE)
    _sphere(lesion, (30, 20, 20), 9)
    volume = lesion.sum()
    df = _run(lesion, tiv=volume / 0.02)          # the lesion is 2% of the TIV (>= 1.35%), and it survives erosion
    assert df.iloc[0]["coalescence_score"] == 4
    small = np.zeros(SHAPE)
    _cube(small, 10)
    assert _run(small).iloc[0]["coalescence_score"] == 1                 # a small lesion that survives erosion


def test_mesh_surface_area_is_close_to_the_true_area_of_a_sphere():
    lesion = np.zeros(SHAPE)
    _sphere(lesion, (30, 20, 20), 9)
    volume = float(lesion.sum())
    radius = (3 * volume / (4 * np.pi)) ** (1 / 3)
    true_area = 4 * np.pi * radius ** 2
    padded = np.pad(lesion > 0, 1)
    # marching cubes on a binary mask keeps a small residual bias (+8% here), against +50% for exposed voxel faces
    assert _mesh_surface_area(padded, (1.0, 1.0, 1.0)) == pytest.approx(true_area, rel=0.12)
    faces = (np.diff(padded.astype(int), axis=0) != 0).sum() + (np.diff(padded.astype(int), axis=1) != 0).sum() \
        + (np.diff(padded.astype(int), axis=2) != 0).sum()
    assert faces > 1.3 * true_area                                       # the old exposed-face count overestimated
    assert abs(_mesh_surface_area(padded, (1.0, 1.0, 1.0)) - true_area) < abs(faces - true_area) / 3


def test_fractal_dimension_of_a_solid_body_is_about_three_and_tiny_lesions_get_zero():
    lesion = np.zeros(SHAPE)
    _sphere(lesion, (30, 20, 20), 9)
    coords = np.argwhere(lesion > 0)
    assert 2.4 < _calculate_fractal_dimension(coords, SHAPE) < 3.3
    assert _calculate_fractal_dimension(coords[:4], SHAPE) == 0.0         # under 6 voxels
    assert _calculate_fractal_dimension(np.argwhere(np.ones((2, 2, 2))), SHAPE) == 0.0   # extent under 4 voxels


def test_the_three_best_matching_tracts_with_their_overlap_fractions():
    lesion = np.zeros(SHAPE)
    _cube(lesion, 10)                              # x 10..14 (5 slabs of 25 voxels), y,z 5..9
    tract_a = np.zeros(SHAPE, dtype=bool); tract_a[10:13, 5:10, 5:10] = True      # 3 of 5 slabs: 60%
    tract_b = np.zeros(SHAPE, dtype=bool); tract_b[13:14, 5:10, 5:10] = True      # 1 of 5 slabs: 20%
    row = _run(lesion, tracts={"a": tract_a, "b": tract_b, "far": np.zeros(SHAPE, dtype=bool)}).iloc[0]
    assert row["tract_match"] == "a" and row["tract_match_overlap_fraction"] == pytest.approx(0.6)
    assert row["tract_match_2"] == "b" and row["tract_match_2_overlap_fraction"] == pytest.approx(0.2)
    assert row["tract_match_3"] == "No Tract" and row["tract_match_3_overlap_fraction"] == 0.0


def test_a_lesion_touching_no_tract_is_no_tract():
    lesion = np.zeros(SHAPE)
    _cube(lesion, 10)
    row = _run(lesion).iloc[0]
    assert row["tract_match"] == "No Tract" and row["tract_match_overlap_fraction"] == 0.0
