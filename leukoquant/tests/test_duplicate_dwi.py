"""Tests for choosing which ADNI DWI archive of one session is used (no DICOM converter needed)."""
import logging
import zipfile

import pytest

from leukoquant.utils.duplicate_dwi import Scan, _acquisition_seconds, read_zip_info, select_scans

BVALS = (0.0, 0.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0)
BVECS = tuple([0.0] * 6 + [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0])
SHORT_BVALS = BVALS[:5]  # start of the same gradient table: 2 b0 + 3 diffusion volumes


def scan(shape=(116, 116, 81, 9), zooms=(2.0, 2.0, 2.0), origin=(116.0, 140.0, -100.0), bvals=BVALS, bvecs=BVECS,
         pe="j", acq=1000.0):
    return Scan(shape=shape, zooms=zooms, origin=origin, bvals=bvals, bvecs=bvecs, pe_dir=pe, acq_time=acq)


def short_scan(**kw):
    """A truncated run: the first volumes of the main scan."""
    return scan(shape=(116, 116, 81, 5), bvals=SHORT_BVALS, bvecs=BVECS[:15], **kw)


def b0_scan(pe, **kw):
    return scan(shape=(116, 116, 81), bvals=(0.0,), bvecs=(0.0, 0.0, 0.0), pe=pe, **kw)


def make_zip(folder, name, series, series_no):
    """Zip holding one empty file named like an ADNI DICOM."""
    path = folder / name
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr(f"ADNI_021_S_7045_MR_{series}_br_raw_20220223120056490_2300_S{series_no}_I{name[1:8]}.dcm", b"")
    return str(path)


def choose(paths, scans, **kw):
    return select_scans(paths, scans.__getitem__, **kw)


def test_read_zip_info(tmp_path):
    info = read_zip_info(make_zip(tmp_path, "I1547797.zip", "Axial_DTI_", 1107498))
    assert info.series == "Axial_DTI_" and info.series_no == 1107498


def test_non_adni_or_non_zip_input_is_never_touched(tmp_path):
    other = tmp_path / "dwi.nii.gz"
    other.write_bytes(b"")
    assert read_zip_info(str(other)) is None
    assert select_scans([str(other)]).paths == [str(other)]


@pytest.mark.parametrize("text, expected", [
    ("16:30:7.117500", 16 * 3600 + 30 * 60 + 7.1175),       # not zero padded, as written by some scanners
    ("08:51:04.555", 8 * 3600 + 51 * 60 + 4.555),
    ("", 0.0),
    ("garbage", 0.0),
])
def test_acquisition_time_parsing(text, expected):
    assert _acquisition_seconds({"AcquisitionTime": text}) == pytest.approx(expected)


def test_acquisition_datetime_orders_scans_on_different_days():
    day1 = _acquisition_seconds({"AcquisitionDateTime": "2021-07-16T23:00:00.0"})
    day2 = _acquisition_seconds({"AcquisitionDateTime": "2021-07-17T01:00:00.0"})
    assert day2 > day1


def test_orig_copy_of_the_same_acquisition_loses_even_when_it_is_the_later_series(tmp_path):
    # 021-s-7045: identical acquisition (same real time), the plain series has the LOWER series number.
    plain = make_zip(tmp_path, "I1.zip", "Axial_DTI_", 1107498)
    orig = make_zip(tmp_path, "I2.zip", "ORIG_Axial_DTI", 1107499)
    assert choose([plain, orig], {plain: scan(), orig: scan()}).paths == [plain]


def test_orig_is_matched_in_any_case(tmp_path):
    a = make_zip(tmp_path, "I1.zip", "orig_Axial_DTI", 5)
    b = make_zip(tmp_path, "I2.zip", "Axial_DTI_", 4)
    assert choose([a, b], {a: scan(), b: scan()}).paths == [b]


def test_identical_scans_keep_the_latest_acquisition(tmp_path):
    a = make_zip(tmp_path, "I1.zip", "Axial_DTI_", 11)   # higher series number but acquired first
    b = make_zip(tmp_path, "I2.zip", "Axial_DTI_", 10)
    assert choose([a, b], {a: scan(acq=100.0), b: scan(acq=2000.0)}).paths == [b]


def test_truncated_first_attempt_is_dropped_and_the_full_scan_kept(tmp_path):
    # 036-s-6878 / 052-s-6832 / 032-s-6602: the small scan is the start of the big one's gradient table.
    full = make_zip(tmp_path, "I1.zip", "Axial_MB_DTI", 949778)
    part = make_zip(tmp_path, "I2.zip", "Axial_MB_DTI", 949779)
    scans = {full: scan(acq=2000.0), part: short_scan(acq=1000.0)}
    assert choose([full, part], scans).paths == [full]
    assert choose([part, full], scans).paths == [full]  # independent of the order given


def test_different_voxel_size_keeps_the_isotropic_scan(tmp_path):
    # 024-s-6005: slab tops within 0.3 mm of each other, 2.0 mm cubic voxels versus 2.52 mm in-plane.
    iso = make_zip(tmp_path, "I1.zip", "Axial_DTI_", 544997)
    aniso = make_zip(tmp_path, "I2.zip", "Axial_DTI_", 544999)
    scans = {iso: scan(origin=(114.6, 141.8, -74.9)), aniso: scan(zooms=(2.52, 2.52, 2.0), origin=(145.2, 167.9, -74.6), acq=2000.0)}
    assert choose([iso, aniso], scans).paths == [iso]
    assert choose([aniso, iso], scans).paths == [iso]


def test_different_slab_keeps_the_brain_slab(tmp_path):
    # 052-s-4944: the later series is the LOWER slab (top at -37 mm versus +55 mm); the brain slab must win.
    brain = make_zip(tmp_path, "I1.zip", "Axial_MB_DTI", 678957)
    lower = make_zip(tmp_path, "I2.zip", "Axial_MB_DTI", 678960)
    scans = {brain: scan(origin=(116.0, 145.1, -105.4)), lower: scan(origin=(116.0, 146.4, -197.3), acq=2000.0)}
    assert choose([brain, lower], scans).paths == [brain]


def test_most_diffusion_volumes_outranks_slab_position(tmp_path):
    big = make_zip(tmp_path, "I1.zip", "Axial_DTI_", 10)
    high = make_zip(tmp_path, "I2.zip", "Axial_DTI_", 11)
    scans = {big: scan(origin=(116.0, 140.0, -200.0)), high: short_scan(origin=(116.0, 140.0, 0.0))}
    assert choose([big, high], scans).paths == [big]


def test_different_phase_encoding_direction_keeps_one_scan(tmp_path):
    a = make_zip(tmp_path, "I1.zip", "Axial_DTI_", 10)
    b = make_zip(tmp_path, "I2.zip", "Axial_DTI_", 11)
    assert choose([a, b], {a: scan(pe="j"), b: scan(pe="i", acq=2000.0)}).paths == [b]  # tie -> latest acquisition


def _073_session(tmp_path):
    dwi = make_zip(tmp_path, "I1.zip", "Axial_MB_DTI_low_sar", 1044829)
    pa = make_zip(tmp_path, "I2.zip", "Axial_MB_DTI_b0_PA", 1044824)
    ap = make_zip(tmp_path, "I3.zip", "Axial_MB_DTI_b0_AP", 1044836)
    return dwi, pa, ap, {dwi: scan(), pa: b0_scan("j"), ap: b0_scan("j-")}


def test_single_b0_scans_are_dropped_by_default(tmp_path):
    dwi, pa, ap, scans = _073_session(tmp_path)
    result = choose([pa, dwi, ap], scans)
    assert result.paths == [dwi] and not result.reverse_polarity


def test_topup_auto_keeps_one_b0_of_each_direction(tmp_path):
    dwi, pa, ap, scans = _073_session(tmp_path)
    result = choose([pa, dwi, ap], scans, topup="auto")
    assert result.paths == [pa, dwi, ap] and result.reverse_polarity


def test_topup_auto_needs_both_directions(tmp_path):
    dwi, pa, ap, scans = _073_session(tmp_path)
    result = choose([pa, dwi], scans, topup="auto")
    assert result.paths == [dwi] and not result.reverse_polarity


def test_topup_auto_ignores_b0_on_another_grid(tmp_path):
    dwi, pa, ap, scans = _073_session(tmp_path)
    scans[pa] = scan(shape=(128, 128, 81), bvals=(0.0,), bvecs=(0.0, 0.0, 0.0), pe="j")
    result = choose([pa, dwi, ap], scans, topup="auto")
    assert result.paths == [dwi] and not result.reverse_polarity


def test_topup_auto_uses_the_latest_b0_of_a_direction(tmp_path):
    dwi, pa, ap, scans = _073_session(tmp_path)
    pa2 = make_zip(tmp_path, "I4.zip", "Axial_MB_DTI_b0_PA", 1044999)
    scans[pa2] = b0_scan("j", acq=5000.0)
    result = choose([pa, pa2, dwi, ap], scans, topup="auto")
    assert result.paths == [pa2, dwi, ap]


def test_invalid_topup_value_is_rejected(tmp_path):
    with pytest.raises(ValueError):
        select_scans([], topup="always")


def test_original_order_is_preserved_for_kept_scans(tmp_path):
    dwi, pa, ap, scans = _073_session(tmp_path)
    assert choose([ap, dwi, pa], scans, topup="auto").paths == [ap, dwi, pa]


def test_unreadable_archive_is_left_in_the_list(tmp_path):
    a = make_zip(tmp_path, "I1.zip", "Axial_DTI_", 1)
    b = make_zip(tmp_path, "I2.zip", "Axial_DTI_", 2)

    def failing(path):
        raise RuntimeError("converter failed")

    assert select_scans([a, b], failing).paths == [a, b]


def test_a_session_with_no_diffusion_scan_is_left_alone(tmp_path):
    ap = make_zip(tmp_path, "I1.zip", "Axial_MB_DTI_b0_AP", 1)
    pa = make_zip(tmp_path, "I2.zip", "Axial_MB_DTI_b0_PA", 2)
    assert choose([ap, pa], {ap: b0_scan("j-"), pa: b0_scan("j")}).paths == [ap, pa]


def test_scan_with_too_few_diffusion_volumes_is_kept_but_warned_about(tmp_path, caplog):
    # 052-s-6305 style fragments: the best one still has fewer than 6 diffusion volumes.
    a = make_zip(tmp_path, "I1.zip", "Axial_MB_DTI", 1)
    b = make_zip(tmp_path, "I2.zip", "Axial_MB_DTI", 2)
    tiny = dict(shape=(116, 116, 81, 4), bvals=(0.0, 0.0, 1000.0, 2000.0), bvecs=(0.0,) * 6 + (1.0, 0.0, 0.0, 0.0, 1.0, 0.0))
    with caplog.at_level(logging.WARNING):
        result = choose([a, b], {a: scan(**tiny), b: scan(shape=(116, 116, 81, 2), bvals=(0.0, 0.0), bvecs=(0.0,) * 6)})
    assert result.paths == [a]
    assert "at least 6" in caplog.text


# ---- stacking scans that add new directions -------------------------------------------------------------------------
from leukoquant.utils.duplicate_dwi import _per_volume_bvecs

# b=1000 directions that are not in BVECS (BVECS holds the x, y, z axes and two diagonals)
NEW_BVALS = (0.0, 1000.0, 1000.0, 1000.0)
NEW_BVECS = (0.0, 0.0, 0.0, 0.577, 0.577, 0.577, 0.577, -0.577, 0.577, -0.577, 0.577, 0.577)


def complementary(**kw):
    return scan(shape=(116, 116, 81, 4), bvals=NEW_BVALS, bvecs=NEW_BVECS, **kw)


def test_scan_with_new_directions_on_the_same_grid_is_stacked(tmp_path):
    main = make_zip(tmp_path, "I1.zip", "Axial_DTI_", 10)
    extra = make_zip(tmp_path, "I2.zip", "Axial_DTI_", 11)
    assert choose([main, extra], {main: scan(), extra: complementary(acq=2000.0)}).paths == [main, extra]


def test_stacking_does_not_depend_on_the_order_given(tmp_path):
    main = make_zip(tmp_path, "I1.zip", "Axial_DTI_", 10)
    extra = make_zip(tmp_path, "I2.zip", "Axial_DTI_", 11)
    scans = {main: scan(), extra: complementary(acq=2000.0)}
    assert choose([extra, main], scans).paths == [extra, main]


def test_new_directions_on_another_slab_or_phase_encoding_are_not_stacked(tmp_path):
    main = make_zip(tmp_path, "I1.zip", "Axial_DTI_", 10)
    lower = make_zip(tmp_path, "I2.zip", "Axial_DTI_", 11)
    other_pe = make_zip(tmp_path, "I3.zip", "Axial_DTI_", 12)
    scans = {main: scan(), lower: complementary(origin=(116.0, 140.0, -190.0)), other_pe: complementary(pe="i")}
    assert choose([main, lower, other_pe], scans).paths == [main]


def test_new_directions_on_another_voxel_grid_are_not_stacked(tmp_path):
    main = make_zip(tmp_path, "I1.zip", "Axial_DTI_", 10)
    coarse = make_zip(tmp_path, "I2.zip", "Axial_DTI_", 11)
    assert choose([main, coarse], {main: scan(), coarse: complementary(zooms=(2.5, 2.5, 2.0))}).paths == [main]


def test_same_directions_with_a_different_shell_count_as_new(tmp_path):
    main = make_zip(tmp_path, "I1.zip", "Axial_DTI_", 10)
    other = make_zip(tmp_path, "I2.zip", "Axial_DTI_", 11)
    # the x axis again, but at b=2000: a new measurement, not a repeat
    shell2 = scan(shape=(116, 116, 81, 2), bvals=(0.0, 2000.0), bvecs=(0.0, 0.0, 0.0, 1.0, 0.0, 0.0))
    assert choose([main, other], {main: scan(), other: shell2}).paths == [main, other]


def test_directions_within_a_few_degrees_count_as_repeats(tmp_path):
    main = make_zip(tmp_path, "I1.zip", "Axial_DTI_", 10)
    other = make_zip(tmp_path, "I2.zip", "Axial_DTI_", 11)
    # x axis tilted by about 1 degree, opposite sign (b-vectors are axial: +v and -v are the same direction)
    near = scan(shape=(116, 116, 81, 2), bvals=(0.0, 1000.0), bvecs=(0.0, 0.0, 0.0, -1.0, 0.017, 0.0))
    assert choose([main, other], {main: scan(), other: near}).paths == [main]


def test_three_scans_stack_only_the_ones_that_add_something(tmp_path):
    main = make_zip(tmp_path, "I1.zip", "Axial_DTI_", 10)
    extra = make_zip(tmp_path, "I2.zip", "Axial_DTI_", 11)
    repeat = make_zip(tmp_path, "I3.zip", "Axial_DTI_", 12)
    scans = {main: scan(), extra: complementary(acq=2000.0), repeat: short_scan(acq=3000.0)}
    assert choose([main, extra, repeat], scans).paths == [main, extra]


def test_bvecs_are_read_per_volume():
    # dcm2niix writes 3 rows (x, y, z) of n values; the helper turns that into n (x, y, z) triplets.
    import numpy as np
    rows = np.array([[0.0, 1.0, 0.0, 0.5], [0.0, 0.0, 1.0, 0.5], [0.0, 0.0, 0.0, 0.7]])
    expected = (0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 0.7)
    assert _per_volume_bvecs(rows, 4) == expected
    assert _per_volume_bvecs(rows.T, 4) == expected  # already one row per volume
    assert _per_volume_bvecs(np.array([0.0, 0.0, 0.0]), 1) == (0.0, 0.0, 0.0)  # a single volume
