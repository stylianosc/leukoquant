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
