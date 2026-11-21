"""Tests for resolving pipeline subject IDs against a demographics table.

The pipeline identifies a subject-session by its output-directory path, so
multi-session cohorts use composite IDs while single-session cohorts use a plain
subject ID. The harmonised clinical_data_nn.csv schema splits these across
`subject`, `session` and `subject_session`, so a composite ID never matches
`subject` on its own.

That mismatch silently prevented covariate-adjusted pseudo-healthy z-scores from
ever running for OASIS-3 and EPAD: their healthy-reference lists use composite
IDs, so every healthy subject looked like it was missing age/sex and the run
aborted before submitting anything. ADNI-3 was unaffected purely because its IDs
are subject-only, which is why the bug went unnoticed -- hence a test covering
all three real conventions rather than just one.
"""
import pandas as pd
import pytest

from leukoquant.utils.z_score_utils import demographics_rows_for_subject


def _demo(rows):
    return pd.DataFrame(rows, columns=["subject", "session", "subject_session", "age", "sex_binary"])


OASIS3 = _demo([
    ["sub-OAS30005", "ses-d1274", "sub-OAS30005_ses-d1274", 71.0, 1],
    ["sub-OAS30005", "ses-d9999", "sub-OAS30005_ses-d9999", 74.0, 1],
    ["sub-OAS30052", "ses-d2737", "sub-OAS30052_ses-d2737", 68.0, 0],
])

EPAD = _demo([
    ["sub-011EPAD23687", "ses-01", "sub-011EPAD23687_ses-01", 64.0, 0],
    ["sub-011EPAD23687", "ses-03", "sub-011EPAD23687_ses-03", 66.0, 0],
])

ADNI3 = _demo([
    ["subj-002-s-6404", "sess-2017-03-31", "subj-002-s-6404-sess-2017-03-31", 70.0, 1],
    ["subj-003-s-6307", "sess-2017-05-02", "subj-003-s-6307-sess-2017-05-02", 73.0, 0],
])


def test_oasis3_composite_id_matches_on_session():
    """OASIS-3 IDs are "<subject>/<session>"; the tail is the bare session."""
    rows = demographics_rows_for_subject(OASIS3, "sub-OAS30005/ses-d1274")
    assert len(rows) == 1
    assert rows.iloc[0]["age"] == 71.0


def test_epad_composite_id_matches_on_subject_session():
    """EPAD IDs are "<subject>/<subject>_<session>"; the tail is subject_session."""
    rows = demographics_rows_for_subject(EPAD, "sub-011EPAD23687/sub-011EPAD23687_ses-01")
    assert len(rows) == 1
    assert rows.iloc[0]["age"] == 64.0


def test_adni3_plain_subject_id_unchanged():
    """Single-session cohorts keep the original exact-match behaviour."""
    rows = demographics_rows_for_subject(ADNI3, "subj-002-s-6404")
    assert len(rows) == 1
    assert rows.iloc[0]["age"] == 70.0


def test_composite_id_selects_the_right_session():
    """The whole point of resolving the tail: two sessions of one subject must
    not collapse into the same row."""
    first = demographics_rows_for_subject(OASIS3, "sub-OAS30005/ses-d1274")
    second = demographics_rows_for_subject(OASIS3, "sub-OAS30005/ses-d9999")
    assert first.iloc[0]["age"] == 71.0
    assert second.iloc[0]["age"] == 74.0


def test_unknown_subject_returns_empty():
    assert demographics_rows_for_subject(OASIS3, "sub-NOPE/ses-d0001").empty
    assert demographics_rows_for_subject(ADNI3, "subj-999-s-0000").empty


def test_known_subject_unknown_session_returns_empty():
    """A real subject with a session the table does not cover is still a miss --
    silently falling back to another session would attach the wrong age."""
    assert demographics_rows_for_subject(OASIS3, "sub-OAS30005/ses-dXXXX").empty


def test_missing_id_column_returns_empty_not_keyerror():
    df = pd.DataFrame({"not_subject": ["x"], "age": [70.0]})
    assert demographics_rows_for_subject(df, "sub-OAS30005/ses-d1274").empty


def test_schema_without_session_columns_falls_back_to_subject():
    """A demographics file carrying no session information cannot disambiguate,
    so the subject match is accepted rather than failing outright."""
    df = pd.DataFrame({"subject": ["sub-OAS30005"], "age": [71.0]})
    rows = demographics_rows_for_subject(df, "sub-OAS30005/ses-d1274")
    assert len(rows) == 1


def test_whitespace_in_id_is_tolerated():
    rows = demographics_rows_for_subject(OASIS3, "  sub-OAS30005/ses-d1274  ")
    assert len(rows) == 1


@pytest.mark.parametrize("cohort,sid", [
    (OASIS3, "sub-OAS30005/ses-d1274"),
    (EPAD, "sub-011EPAD23687/sub-011EPAD23687_ses-01"),
    (ADNI3, "subj-002-s-6404"),
])
def test_covariates_are_reachable_for_every_cohort_convention(cohort, sid):
    """The regression this guards: the covariate-completeness gate asks whether
    age and sex_binary are present for each healthy subject. Before the fix this
    was False for every OASIS-3 and EPAD subject regardless of the data."""
    rows = demographics_rows_for_subject(cohort, sid)
    assert not rows[["age", "sex_binary"]].isna().any(axis=None)
