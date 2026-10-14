"""Regression tests for transform_trk.py's --convert-affine (header-only re-affine) step.

The atlas tracts (FreeSurfer trctrain/hcp/syn/*.display.trk) are re-labelled with the
MGH35 template header before their density map is built. nibabel returns streamline
points in RAS+mm, so the step must keep them in world space. Reading them as voxel
indices clamped every streamline to the volume boundary and produced empty atlas
tracts after resampling to subject space.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

nib = pytest.importorskip("nibabel")
pytest.importorskip("dipy")

SCRIPT = Path(__file__).resolve().parents[1] / "utils" / "transform_trk.py"

# MGH35-like grid: 1.5 mm voxels, x axis flipped (LAS), origin offset like the real template.
AFFINE = np.array(
    [
        [-1.5, 0.0, 0.0, 90.5],
        [0.0, 1.5, 0.0, -126.5],
        [0.0, 0.0, 1.5, -72.5],
        [0.0, 0.0, 0.0, 1.0],
    ]
)
SHAPE = (122, 146, 122)


def _write_reference(path: Path) -> None:
    nib.save(nib.Nifti1Image(np.zeros(SHAPE, dtype=np.float32), AFFINE), str(path))


def _write_trk(path: Path, streamlines_rasmm) -> None:
    """Save streamlines (RAS+mm) with a trk header on the same grid as the reference."""
    tractogram = nib.streamlines.Tractogram(streamlines_rasmm, affine_to_rasmm=np.eye(4))
    header = {
        nib.streamlines.Field.VOXEL_TO_RASMM: AFFINE,
        nib.streamlines.Field.VOXEL_SIZES: (1.5, 1.5, 1.5),
        nib.streamlines.Field.DIMENSIONS: SHAPE,
    }
    nib.streamlines.save(tractogram, str(path), header=header)


def test_convert_affine_preserves_world_coordinates(tmp_path):
    """Streamlines with negative RAS coordinates must come out unmoved, not clamped to a boundary."""
    # A streamline well inside the volume whose world coordinates include negative x and y.
    t = np.linspace(0, 1, 20)[:, None]
    sl = (np.array([-30.0, -40.0, 10.0]) + t * np.array([10.0, 20.0, 15.0])).astype(np.float32)

    in_trk = tmp_path / "atlas.display.trk"
    out_trk = tmp_path / "atlas_mgh35.trk"
    ref = tmp_path / "template.nii.gz"
    _write_trk(in_trk, [sl])
    _write_reference(ref)

    result = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--input-trk", str(in_trk),
            "--output-trk", str(out_trk),
            "--in-reference", str(ref),
            "--convert-affine",
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    # No point should have been clamped: the tract lies fully inside the grid.
    assert "Clamped" not in result.stderr

    out = nib.streamlines.load(str(out_trk))
    assert len(out.streamlines) == 1
    np.testing.assert_allclose(out.streamlines[0], sl, atol=1e-3)


def test_clamped_points_survive_save_and_density_map_on_a_flipped_60_slice_grid(tmp_path):
    """EPAD's FA grid has 60 slices and a flipped/oblique orientation. Points clamped to the field-of-view edge used
    to land at 59.5000x after the float32 .trk round trip, and dipy's density_map then indexed slice 60 of 60
    ("index 60 is out of bounds for axis 2 with size 60")."""
    from dipy.io.stateful_tractogram import Space
    from dipy.io.streamline import load_tractogram, save_tractogram
    from dipy.tracking.utils import density_map

    from leukoquant.utils.transform_trk import _make_sft_safe

    rng = np.random.default_rng(0)
    for trial in range(30):
        angle = rng.uniform(-0.4, 0.4)
        rot = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        affine = np.eye(4)
        affine[:3, :3] = rot @ np.diag([-1.8, 1.8, 2.0])
        affine[:3, 3] = [100.0, -90.0, -60.0]
        shape = (110, 110, 60)
        ref = tmp_path / f"ref{trial}.nii.gz"
        nib.save(nib.Nifti1Image(np.zeros(shape, dtype=np.float32), affine), str(ref))
        vox = np.stack([rng.uniform(5, 100, 40), rng.uniform(5, 100, 40), rng.uniform(40, 75, 40)], axis=1)
        points_ras = nib.affines.apply_affine(affine, vox)
        sft = _make_sft_safe([points_ras[:20], points_ras[20:]], str(ref), Space.RASMM)
        out = tmp_path / f"out{trial}.trk"
        save_tractogram(sft, str(out), bbox_valid_check=False)
        back = load_tractogram(str(out), str(ref), bbox_valid_check=False)
        density_map(back.streamlines, back.affine, back.dimensions)
