import numpy as np
import nibabel as nib
from pathlib import Path

import pytest

from leukoquant.utils.image_utils import main as image_utils_main
from leukoquant.utils.image_utils import reorient_to_ras, same_voxel_grid, select_first_frame


def test_select_first_frame_on_3d_image(tmp_path: Path):
    data = np.zeros((10, 10, 10), dtype=np.float32)
    affine = np.eye(4)
    img = nib.Nifti1Image(data, affine)
    src = tmp_path / "img_3d.nii.gz"
    dst = tmp_path / "out_3d.nii.gz"
    nib.save(img, str(src))

    select_first_frame(str(src), str(dst))

    out_img = nib.load(str(dst))
    assert out_img.ndim == 3
    assert out_img.shape == (10, 10, 10)


def test_select_first_frame_on_4d_image(tmp_path: Path):
    data = np.zeros((10, 10, 10, 2), dtype=np.float32)
    data[..., 0] = 1.0
    data[..., 1] = 2.0
    affine = np.eye(4)
    img = nib.Nifti1Image(data, affine)
    src = tmp_path / "img_4d.nii.gz"
    dst = tmp_path / "out_collapsed.nii.gz"
    nib.save(img, str(src))

    select_first_frame(str(src), str(dst))

    out_img = nib.load(str(dst))
    assert out_img.ndim == 3
    assert out_img.shape == (10, 10, 10)
    assert np.allclose(out_img.get_fdata(), 1.0)


def test_select_first_frame_inplace(tmp_path: Path):
    # In-place on 4D
    data = np.zeros((8, 8, 8, 3), dtype=np.float32)
    data[..., 0] = 5.0
    img = nib.Nifti1Image(data, np.eye(4))
    f = tmp_path / "inplace_4d.nii.gz"
    nib.save(img, str(f))

    select_first_frame(str(f), str(f))

    out_img = nib.load(str(f))
    assert out_img.ndim == 3
    assert out_img.shape == (8, 8, 8)
    assert np.allclose(out_img.get_fdata(), 5.0)

    # In-place on 3D (no-op)
    select_first_frame(str(f), str(f))
    out_img2 = nib.load(str(f))
    assert out_img2.ndim == 3
    assert np.allclose(out_img2.get_fdata(), 5.0)


def test_reorient_to_ras_inplace(tmp_path: Path):
    data = np.ones((6, 6, 6), dtype=np.float32)
    # create non-RAS affine
    affine = np.diag([-1, 1, 1, 1])
    img = nib.Nifti1Image(data, affine)
    f = tmp_path / "inplace_reorient.nii.gz"
    nib.save(img, str(f))

    reorient_to_ras(str(f), str(f))

    out_img = nib.load(str(f))
    assert nib.aff2axcodes(out_img.affine) == ("R", "A", "S")


def _write_psr_4d_t1(path: Path) -> np.ndarray:
    """A 4D T1 with a duplicate second frame stored P,S,R, like the raw EPAD sessions that broke the metrics step."""
    rng = np.random.default_rng(0)
    vol = rng.random((6, 7, 5)).astype(np.float32)
    vol[0, 0, 0] = -3.0  # a negative voxel: pure reorientation must leave it alone
    # Axis 0 runs towards Posterior, axis 1 towards Superior, axis 2 towards Right (axcodes P,S,R).
    affine = np.array([[0, 0, 1, 0], [-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=float)
    assert nib.aff2axcodes(affine) == ("P", "S", "R")
    nib.save(nib.Nifti1Image(np.stack([vol, vol], axis=-1), affine), str(path))
    return vol


def test_first_frame_then_canonical_gives_3d_ras_on_the_bamos_grid(tmp_path: Path):
    """The metrics stage must end up on the same grid as BaMoS (first frame, then RAS+)."""
    src, single, canon, bamos = (tmp_path / n for n in ("t1.nii.gz", "single.nii.gz", "canon.nii.gz", "bamos.nii.gz"))
    _write_psr_4d_t1(src)

    select_first_frame(str(src), str(single))
    reorient_to_ras(str(single), str(canon), clip_negatives=False)
    # What BaMoS is given: the same two steps, but with negatives clipped.
    reorient_to_ras(str(single), str(bamos))

    out, ref = nib.load(str(canon)), nib.load(str(bamos))
    assert out.ndim == 3
    assert nib.aff2axcodes(out.affine) == ("R", "A", "S")
    # Same voxel grid as the BaMoS input, so a lesion map from BaMoS passes the shape check.
    assert out.shape == ref.shape
    assert np.allclose(out.affine, ref.affine)


def test_reorient_without_clipping_keeps_negative_voxels(tmp_path: Path):
    src, single, canon = (tmp_path / n for n in ("t1.nii.gz", "single.nii.gz", "canon.nii.gz"))
    _write_psr_4d_t1(src)
    select_first_frame(str(src), str(single))

    reorient_to_ras(str(single), str(canon), clip_negatives=False)

    assert np.asarray(nib.load(str(canon)).dataobj).min() == -3.0


def test_reorient_without_clipping_copies_an_already_ras_image(tmp_path: Path):
    src, dst = tmp_path / "ras.nii.gz", tmp_path / "out.nii.gz"
    nib.save(nib.Nifti1Image(np.full((4, 4, 4), -1.0, dtype=np.float32), np.eye(4)), str(src))

    reorient_to_ras(str(src), str(dst), clip_negatives=False)

    assert src.read_bytes() == dst.read_bytes()


def _save(path: Path, shape, affine) -> str:
    nib.save(nib.Nifti1Image(np.zeros(shape, dtype=np.float32), affine), str(path))
    return str(path)


def test_same_voxel_grid_accepts_identical_and_rejects_flip_and_shape(tmp_path: Path):
    ras = np.diag([1.0, 1.0, 1.0, 1.0])
    las = np.diag([-1.0, 1.0, 1.0, 1.0])
    ref = _save(tmp_path / "ref.nii.gz", (6, 7, 8), ras)
    assert same_voxel_grid(_save(tmp_path / "same.nii.gz", (6, 7, 8), ras), ref)
    # Same shape, mirrored first axis: this is the LAS vs RAS case that used to pass silently.
    assert not same_voxel_grid(_save(tmp_path / "flip.nii.gz", (6, 7, 8), las), ref)
    # A DWI-sized image is never on the T1 grid.
    assert not same_voxel_grid(_save(tmp_path / "dwi.nii.gz", (3, 3, 2), ras), ref)


def test_same_grid_cli_exit_status(tmp_path: Path):
    ref = _save(tmp_path / "ref.nii.gz", (4, 4, 4), np.eye(4))
    same = _save(tmp_path / "same.nii.gz", (4, 4, 4), np.eye(4))
    other = _save(tmp_path / "other.nii.gz", (5, 4, 4), np.eye(4))

    with pytest.raises(SystemExit) as ok:
        image_utils_main(["--input", same, "--ref", ref, "--op", "same_grid"])
    assert ok.value.code == 0
    with pytest.raises(SystemExit) as diff:
        image_utils_main(["--input", other, "--ref", ref, "--op", "same_grid"])
    assert diff.value.code == 1
