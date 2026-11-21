"""Lesion shape eigenvalues must be reproducible bit for bit, and a flat lesion's ratios must not be rounding noise.

Background (EPAD determinism check, 2026-11-18): two reruns of the metrics step with identical code and inputs gave
eigenvalue ratios up to 19% apart for flat lesions. The smallest PCA eigenvalue of a lesion whose voxels lie in one
plane is mathematically zero, but the floating-point PCA returned rounding noise (about 1e-15 mm^2) that depended on
summation order, and the ratios divided by it. These tests pin the replacement: exact integer second moments, an
eigenvalue at rounding level reported as exactly 0, and NaN for a ratio whose denominator is 0.
"""

import numpy as np
import pytest

from leukoquant.utils.lesion_metrics import _compute_pca_eigenvalues


def _ellipsoid_coords(a=6, b=4, c=2):
    """Integer voxel coordinates of a solid axis-aligned ellipsoid (semi-axes a > b > c voxels)."""
    xx, yy, zz = np.mgrid[-a:a + 1, -b:b + 1, -c:c + 1]
    inside = (xx / a) ** 2 + (yy / b) ** 2 + (zz / c) ** 2 <= 1.0
    return np.stack([xx[inside], yy[inside], zz[inside]], axis=1) + 20


def test_eigenvalues_match_sample_covariance():
    coords = _ellipsoid_coords()
    expected = np.sort(np.linalg.eigvalsh(np.cov(coords.T.astype(float))))[::-1]
    got = np.asarray(_compute_pca_eigenvalues(coords))
    np.testing.assert_allclose(got, expected, rtol=1e-12)
    assert got[0] >= got[1] >= got[2] > 0


def test_voxel_sizes_scale_the_axes():
    coords = _ellipsoid_coords()
    vs = (1.0, 1.0, 2.0)
    expected = np.sort(np.linalg.eigvalsh(np.cov((coords * np.asarray(vs)).T)))[::-1]
    np.testing.assert_allclose(_compute_pca_eigenvalues(coords, vs), expected, rtol=1e-12)


def test_order_of_voxels_does_not_change_a_single_bit():
    coords = _ellipsoid_coords()
    rng = np.random.default_rng(0)
    reference = np.asarray(_compute_pca_eigenvalues(coords, (0.9375, 0.9375, 1.2)))
    for _ in range(5):
        shuffled = coords[rng.permutation(len(coords))]
        assert np.array_equal(np.asarray(_compute_pca_eigenvalues(shuffled, (0.9375, 0.9375, 1.2))), reference)


@pytest.mark.parametrize("plane", ["axial", "diagonal"])
def test_flat_lesion_has_exactly_zero_third_eigenvalue(plane):
    xx, yy = np.mgrid[0:9, 0:5]
    if plane == "axial":                       # all voxels at one z
        coords = np.stack([xx.ravel(), yy.ravel(), np.full(xx.size, 7)], axis=1)
    else:                                      # oblique plane z = x
        coords = np.stack([xx.ravel(), yy.ravel(), xx.ravel()], axis=1)
    eig = _compute_pca_eigenvalues(coords + 10, (1.0, 1.0, 1.0))
    assert eig[2] == 0.0
    assert eig[0] > 0 and eig[1] > 0


def test_line_lesion_has_two_zero_eigenvalues():
    coords = np.stack([np.arange(8), np.full(8, 3), np.full(8, 4)], axis=1)
    eig = _compute_pca_eigenvalues(coords)
    assert eig[0] > 0 and eig[1] == 0.0 and eig[2] == 0.0


def test_fewer_than_three_voxels_gives_nan():
    assert all(np.isnan(_compute_pca_eigenvalues(np.array([[1, 2, 3], [1, 2, 4]]))))


def test_large_coordinates_and_many_voxels_stay_exact():
    # A big lesion far from the origin: the integer moments must not overflow or lose precision.
    coords = _ellipsoid_coords(40, 30, 20) + 400
    expected = np.sort(np.linalg.eigvalsh(np.cov((coords - coords.mean(axis=0)).T)))[::-1]
    np.testing.assert_allclose(_compute_pca_eigenvalues(coords), expected, rtol=1e-10)
