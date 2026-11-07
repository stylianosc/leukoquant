import numpy as np
import nibabel as nib
from pathlib import Path

from leukoquant.utils.image_utils import select_first_frame, reorient_to_ras


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
