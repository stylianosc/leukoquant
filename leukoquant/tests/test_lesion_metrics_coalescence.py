"""Unit tests for lesion_metrics._calculate_coalescence_score.

The score grades a lesion from 0 to 4 by its TIV-normalised volume and by the fraction of its volume that survives
erosion. A lesion that vanishes on erosion scores 0; before the fix the next size test was a separate `if`, so that 0
was always overwritten.
"""

import pytest

from leukoquant.utils.lesion_metrics import _calculate_coalescence_score


def score(tiv_normalised, surviving_fraction, vanished=False):
    """Score of a 100 mm3 lesion whose eroded volume is `surviving_fraction` of the original."""
    volume = 100.0
    return _calculate_coalescence_score(
        tiv_normalised,
        volume * surviving_fraction,
        volume,
        "Vanished" if vanished else "Survived",
    )


def test_vanished_lesion_scores_zero_whatever_its_size():
    assert score(0.01, 0.0, vanished=True) == 0
    assert score(2.0, 0.0, vanished=True) == 0


@pytest.mark.parametrize(
    "tiv_normalised, surviving_fraction, expected",
    [
        (0.01, 0.5, 1),   # small lesion
        (0.10, 0.1, 1),   # 0.05-0.4, mostly eroded
        (0.10, 0.5, 2),   # 0.05-0.4, survives erosion
        (0.50, 0.1, 2),   # 0.4-1.35, mostly eroded
        (0.50, 0.5, 3),   # 0.4-1.35, survives erosion
        (2.00, 0.1, 3),   # >= 1.35, mostly eroded
        (2.00, 0.5, 4),   # >= 1.35, survives erosion
    ],
)
def test_score_bands(tiv_normalised, surviving_fraction, expected):
    assert score(tiv_normalised, surviving_fraction) == expected
