# ============================================================================
# Metric family: per-lesion metrics  (FAST implementation)
#
# Drop-in replacement for leukoquant/utils/lesion_metrics.py. Same public
# function (_compute_lesion_metrics), same output DataFrames (column names,
# column order, dtypes, values). Differences are purely in how the numbers are
# computed:
#
#   * every per-lesion operation works on the lesion's bounding box / its own
#     voxel list instead of the full ~10 M voxel volume;
#   * fractal dimension (the dominant cost) counts occupied boxes from the
#     lesion's voxel coordinates instead of histogramdd over the whole padded
#     volume 216 x 14 times per lesion;
#   * surface area is a vectorised exposed-face count instead of a Python
#     triple loop over every voxel of the volume;
#   * the original->dilated/eroded label bookkeeping uses bincount/unique on
#     lesion voxels instead of a Python loop over every voxel of the volume;
#   * per-lesion map statistics gather the lesion's voxels once (sorted by
#     label) instead of boolean-masking every map for every lesion;
#   * the per-lesion DataFrame is assembled from row dicts instead of ~10^4
#     DataFrame.loc[...] = value assignments.
#
# Legacy quirks are PRESERVED on purpose so results stay identical:
#   * "merged_after_dilation" is True whenever the dilated component holds more
#     than one lesion *voxel* (len() of a per-voxel list), not more than one lesion;
#   * "Remained Single" means exactly ONE eroded voxel; two or more eroded
#     voxels is reported as "Split into smaller lesions" even for a single
#     eroded component (the legacy list holds one entry per voxel);
#   * fractal-dimension box edges use the FULL (padded) volume shape, so
#     voxels beyond the last bin edge are dropped exactly like histogramdd did.
# ============================================================================
import os
import scipy
import scipy.stats
from scipy import ndimage
import skimage
import skimage.measure
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import nibabel as nib


# ----------------------------------------------------------------------------
# Fractal dimension
# ----------------------------------------------------------------------------

_FRACTAL_SCALES = None


def _fractal_scales():
    """The 14 box sizes the legacy implementation used (computed once)."""
    global _FRACTAL_SCALES
    if _FRACTAL_SCALES is None:
        scales = np.logspace(0, 3, num=15, endpoint=True, base=2)
        _FRACTAL_SCALES = np.array([a for a in scales if a != 0.0])
    return _FRACTAL_SCALES


def _count_occupied_boxes(coords, shape, scales, max_pad=5):
    """Occupied-box counts for every (x_pad, y_pad, z_pad) offset and scale.

    Equivalent to the legacy np.pad + np.histogramdd(...) loop, but computed
    from the lesion voxel coordinates only.

    Legacy semantics reproduced exactly: for a padded axis length L the bin
    edges are np.arange(0, L, scale); a voxel belongs to the bin found by
    searchsorted(edges, c, side="right"), a voxel exactly on the last edge goes
    to the last bin, and voxels beyond the last edge are dropped.

    Occupied cells are counted with a boolean occupancy bitmap restricted to the
    lesion's own extent (scatter + count_nonzero), which is linear in the number
    of lesion voxels (np.unique would sort them for each of the 216 x 14 cases).

    Returns an int array of shape (n_offsets, n_scales); offsets are ordered
    x_pad-major (x, then y, then z), the same order as the legacy loops.
    """
    coords = np.asarray(coords, dtype=np.int64)
    pads = range(0, max_pad + 1)
    n_pad = len(pads)
    counts = np.zeros((n_pad ** 3, len(scales)), dtype=np.int64)
    n_vox = len(coords)

    for si, scale in enumerate(scales):
        # per axis and padding: bin index per voxel (-1 = dropped), its extent
        ix = [[None] * n_pad for _ in range(3)]
        span = [[0] * n_pad for _ in range(3)]
        all_valid = [[True] * n_pad for _ in range(3)]
        for d in range(3):
            for p in pads:
                edges = np.arange(0, shape[d] + 2 * p, scale)
                c = coords[:, d] + p
                ncount = np.searchsorted(edges, c, side="right")
                ncount[c == edges[-1]] -= 1
                nbins = len(edges) - 1
                valid = (ncount >= 1) & (ncount <= nbins)
                idx = np.where(valid, ncount - 1, -1)
                all_valid[d][p] = bool(valid.all())
                if valid.any():
                    lo = int(idx[valid].min())
                    idx = np.where(valid, idx - lo, -1)         # re-base to the lesion's own extent
                    span[d][p] = int(idx[valid].max()) + 1
                ix[d][p] = idx.astype(np.int32)
        k = 0
        for px in pads:
            for py in pads:
                a_, b_ = ix[0][px], ix[1][py]
                for pz in pads:
                    c_ = ix[2][pz]
                    da, db, dc = span[0][px], span[1][py], span[2][pz]
                    if da == 0 or db == 0 or dc == 0:
                        counts[k, si] = 0
                    else:
                        if all_valid[0][px] and all_valid[1][py] and all_valid[2][pz]:
                            key = (a_ * db + b_) * dc + c_
                        else:
                            ok = (a_ >= 0) & (b_ >= 0) & (c_ >= 0)
                            key = (a_[ok] * db + b_[ok]) * dc + c_[ok]
                        occ = np.zeros(da * db * dc, dtype=bool)
                        occ[key] = True
                        counts[k, si] = np.count_nonzero(occ)
                    k += 1
    return counts


def _calculate_fractal_dimension(lesion_coords, volume_shape):
    """Box-counting fractal dimension of ONE lesion given its voxel coordinates.

    lesion_coords: (N, 3) integer voxel coordinates of the lesion.
    volume_shape:  shape of the FULL volume the legacy code padded (needed so
                   the box edges, and therefore the dropped-voxel behaviour, match).

    Only box sizes up to half the lesion's longest extent are used: a box larger than the lesion contains all of it, so
    the count is 1 whatever its shape, and including such sizes flattens the fitted slope (most lesions here are a few
    voxels across). At least three box sizes are needed; a lesion that is too small for that, with fewer than 6 voxels
    or with a longest extent under 4 voxels, has no defined dimension and gets 0.0.
    """
    coords = np.asarray(lesion_coords)
    if len(coords) < 6:
        return 0.0
    extent = float((coords.max(axis=0) - coords.min(axis=0) + 1).max())
    if extent < 4:
        return 0.0
    scales = _fractal_scales()
    scales = scales[scales <= extent / 2.0]
    if len(scales) < 3:
        return 0.0

    min_boxes = np.min(_count_occupied_boxes(coords, volume_shape, scales).T, 1)
    if np.any(min_boxes <= 0):
        return 0.0
    # linear fit, polynomial of degree 1; the fractal dimension is the OPPOSITE of the fitting coefficient
    coeffs = np.polyfit(np.log(scales), np.log(min_boxes), 1)
    return float(-1 * coeffs[0])


# ----------------------------------------------------------------------------
# Small helpers (unchanged semantics)
# ----------------------------------------------------------------------------

def _calculate_coalescence_score(lesion_volume_pct_tiv, volume_after_erosion, lesion_volume_mm3, lesion_after_erosion):
    """Ordinal coalescence score 0-4 of one lesion, from its size and how much of it survives one erosion.

    `lesion_volume_pct_tiv` is the lesion volume as a PERCENTAGE of the total intracranial volume (the size bands 0.05, 0.4
    and 1.35 are percentages). The caller used to pass the plain fraction (volume / TIV), about 100 times too small, which
    left every lesion in the lowest size band: only the scores 0 and 1 could occur.
    """
    lesion_volume_tiv_normalized_mm3 = lesion_volume_pct_tiv
    coalescence_score = None
    if lesion_after_erosion == "Vanished":
        coalescence_score = 0
    elif lesion_volume_tiv_normalized_mm3 < 0.05:
        coalescence_score = 1
    elif lesion_volume_tiv_normalized_mm3 >= 1.35:
        if volume_after_erosion / lesion_volume_mm3 >= 0.2:
            coalescence_score = 4
        else:
            coalescence_score = 3
    elif lesion_volume_tiv_normalized_mm3 >= 0.4:
        if volume_after_erosion / lesion_volume_mm3 >= 0.2:
            coalescence_score = 3
        else:
            coalescence_score = 2
    elif lesion_volume_tiv_normalized_mm3 >= 0.05:
        if volume_after_erosion / lesion_volume_mm3 >= 0.2:
            coalescence_score = 2
        else:
            coalescence_score = 1

    return coalescence_score


def _surface_area_from_mask(sub_mask, voxel_sizes=(1.0, 1.0, 1.0)):
    """Exposed-face surface area of a binary mask (vectorised).

    sub_mask must already carry a 1-voxel zero margin (or touch no edge it
    should treat as open) so that faces on the boundary of the array are
    counted: the legacy loop counted a face as exposed when the neighbour was
    background OR outside the volume, which is what a zero margin reproduces.

    Every exposed face separates one mask voxel from one background voxel, so
    the number of exposed faces along an axis equals the number of adjacent
    voxel pairs that differ along that axis.
    """
    dx, dy, dz = voxel_sizes
    nx = int(np.count_nonzero(sub_mask[1:, :, :] != sub_mask[:-1, :, :]))
    ny = int(np.count_nonzero(sub_mask[:, 1:, :] != sub_mask[:, :-1, :]))
    nz = int(np.count_nonzero(sub_mask[:, :, 1:] != sub_mask[:, :, :-1]))
    # Same per-face area values (and float32 rounding of the product) the
    # legacy loop added one by one.
    return nx * (dy * dz) + ny * (dx * dz) + nz * (dx * dy)


def _mesh_surface_area(sub_mask, voxel_sizes=(1.0, 1.0, 1.0)):
    """Surface area (mm^2) of a binary mask from a triangle mesh of its boundary (marching cubes at half a voxel).

    Counting exposed voxel faces measures the staircase of the voxel grid, about 1.5 times the true area for a smooth
    object, so sphericity and compactness came out biased. The mesh is the standard remedy (e.g. in radiomics). The mask
    must carry a zero margin. If a mesh cannot be built (degenerate shapes) the exposed-face count is used.
    """
    try:
        verts, faces, _, _ = skimage.measure.marching_cubes(
            np.asarray(sub_mask, dtype=np.float32), level=0.5, spacing=tuple(float(v) for v in voxel_sizes))
        return float(skimage.measure.mesh_surface_area(verts, faces))
    except (ValueError, RuntimeError):
        return float(_surface_area_from_mask(sub_mask, voxel_sizes))


def _calculate_surface_area(lesion_mask, voxel_sizes=(1.0, 1.0, 1.0)):
    """Surface area of a full-size binary lesion mask (kept for API parity)."""
    padded = np.pad((lesion_mask == 1), 1)
    return _surface_area_from_mask(padded, voxel_sizes)


def _compute_image_stats(values: np.ndarray) -> dict:
    """Summary statistics for a 1D array of values (NaNs for an empty array).

    Same statistics as the region summaries (p5 and p95 instead of min and max); the four percentiles are taken in one
    np.percentile call (one partition instead of four).
    """
    if len(values) == 0:
        return {k: np.nan for k in (
            "mean", "median", "std", "p5", "p95", "25th", "75th", "IQR",
            "peak_width", "kurtosis", "skew")}
    p5, p25, p75, p95 = np.percentile(values, [5, 25, 75, 95])
    return {
        "mean": np.mean(values),
        "median": np.median(values),
        "std": np.std(values),
        "p5": p5,
        "p95": p95,
        "25th": p25,
        "75th": p75,
        "IQR": p75 - p25,
        "peak_width": p95 - p5,
        "kurtosis": scipy.stats.kurtosis(values),
        "skew": scipy.stats.skew(values),
    }


def _compute_pca_eigenvalues(coords: np.ndarray, voxel_sizes: tuple = None):
    """Return PCA explained variance (3 values) for a set of coordinates, or NaNs.

    If `voxel_sizes` is provided the coordinates are scaled to physical units (mm)
    before computing PCA so that returned variances are in mm^2.
    """
    if coords.shape[0] < 3:
        return [np.nan, np.nan, np.nan]

    coords_mm = coords
    if voxel_sizes is not None:
        vs = np.asarray(voxel_sizes)
        coords_mm = coords * vs

    if PCA is not None:
        pca = PCA(n_components=3)
        pca.fit(coords_mm)
        return pca.explained_variance_
    return [np.nan, np.nan, np.nan]


# ----------------------------------------------------------------------------
# Main entry point
# ----------------------------------------------------------------------------

def _compute_lesion_metrics(
    subject: str,
    lesion_array: np.ndarray,
    tiv_volume: float,
    tract_binary_arrays: dict,
    images_dict: dict,
    voxel_sizes: tuple,
    tract_binary_idx: dict = None,
    ):
    """
    Compute per-lesion metrics (see legacy docstring for the output contract).

    lesion_array: 3D (binary) lesion mask
    tiv_volume: TIV volume (float)
    tract_binary_arrays: dict mapping tract name -> binarised tract mask (bool/uint8).
    images_dict: dict of {metric_name: 3D array}
    voxel_sizes: tuple of (x_size, y_size, z_size) in mm
    tract_binary_idx: optional {tract_name: sorted C-order flat indices of the
        binary tract voxels}. Callers that already built these once (the
        tract-level loop does) pass them to avoid recomputing 42 flatnonzero's.
    """
    shape = lesion_array.shape
    voxel_volume = np.prod(voxel_sizes)

    # Binarize lesion array in case it's not already binary
    lesion_array = (lesion_array > 0).astype(np.uint8)

    # Label connected lesions (26-connectivity, as before)
    struct = np.ones((3, 3, 3))
    lesion_labels, lesion_labels_count = ndimage.label(lesion_array, structure=struct)
    regionprops = skimage.measure.regionprops(lesion_labels)

    # Dilate and erode for coalescence (bounded to the lesion bounding box: both
    # operations are local, and the legacy default border_value=0 makes a bbox
    # crop with a 1-voxel margin exactly equivalent)
    nzr = np.nonzero(lesion_array)
    if len(nzr[0]):
        box = tuple(slice(max(int(a.min()) - 2, 0), min(int(a.max()) + 3, s)) for a, s in zip(nzr, shape))
    else:
        box = tuple(slice(0, s) for s in shape)
    sub_lesion = lesion_array[box]
    dil_sub = ndimage.binary_dilation(sub_lesion, iterations=1, structure=struct)
    ero_sub = ndimage.binary_erosion(sub_lesion, iterations=1, structure=struct)
    lab_dil_sub, n_dil = ndimage.label(dil_sub, structure=struct)
    lab_ero_sub, n_ero = ndimage.label(ero_sub, structure=struct)
    lab_sub = lesion_labels[box]

    # ---- original -> dilated / eroded label bookkeeping on lesion voxels only ----
    in_lesion = lab_sub != 0
    o_lab = lab_sub[in_lesion]
    d_lab = lab_dil_sub[in_lesion]
    e_lab = lab_ero_sub[in_lesion]

    # dilated: (number of lesion voxels inside each dilated component, new label per lesion)
    orig_to_dil = np.zeros(lesion_labels_count + 1, dtype=np.int64)
    orig_to_dil[o_lab] = d_lab
    # number of distinct lesions inside each dilated component (the old code counted lesion VOXELS, which is > 1 for every
    # lesion, so "merged" was always true)
    lesions_per_dil = np.bincount(orig_to_dil[1:], minlength=n_dil + 1)
    dil_comp_size = np.bincount(lab_dil_sub.ravel(), minlength=n_dil + 1)

    # eroded: number of eroded voxels per lesion, distinct eroded labels per lesion
    kept = e_lab != 0
    ero_voxel_count = np.bincount(o_lab[kept], minlength=lesion_labels_count + 1)
    ero_comp_size = np.bincount(lab_ero_sub.ravel(), minlength=n_ero + 1)
    pair = np.unique(o_lab[kept].astype(np.int64) * (n_ero + 1) + e_lab[kept])
    pair_orig = pair // (n_ero + 1)
    pair_ero = pair % (n_ero + 1)
    # largest eroded component per lesion (used only when >= 2 eroded voxels)
    max_ero_size = np.zeros(lesion_labels_count + 1, dtype=np.int64)
    np.maximum.at(max_ero_size, pair_orig, ero_comp_size[pair_ero])
    # the (single) eroded label of a lesion with exactly one eroded voxel
    single_ero_label = np.zeros(lesion_labels_count + 1, dtype=np.int64)
    single_ero_label[pair_orig] = pair_ero

    # ---- per-lesion voxel lists (C-order within lesion), bbox, sizes ----
    flat_idx = np.flatnonzero(lesion_labels)             # C-order, ascending
    flat_lab = lesion_labels.reshape(-1)[flat_idx]
    order = np.argsort(flat_lab, kind="stable")           # keeps ascending idx inside a label
    idx_sorted = flat_idx[order]
    lab_sorted = flat_lab[order]
    starts = np.searchsorted(lab_sorted, np.arange(1, lesion_labels_count + 2))
    bboxes = ndimage.find_objects(lesion_labels)

    # overlap of each lesion with each tract: counts[tract, label]
    tract_names = list(tract_binary_arrays.keys())
    counts = np.zeros((len(tract_names), lesion_labels_count + 1), dtype=np.int64)
    for ti, tname in enumerate(tract_names):
        tb = tract_binary_arrays[tname]
        if tb.shape != shape:
            raise ValueError(f"Tract array for {tname} has shape {tb.shape} but expected {shape}")
        if tract_binary_idx is not None and tname in tract_binary_idx:
            tidx = tract_binary_idx[tname]
        else:
            tidx = np.flatnonzero(tb > 0)
        counts[ti] = np.bincount(lesion_labels.reshape(-1)[tidx], minlength=lesion_labels_count + 1)

    lesion_dict = {}

    for lesion_label in range(1, lesion_labels_count + 1):
        sl_box = bboxes[lesion_label - 1]
        lo, hi = starts[lesion_label - 1], starts[lesion_label]
        lesion_volume_voxels = np.uint64(hi - lo)                      # np.sum(uint8 mask) -> uint64
        lesion_volume_mm3 = lesion_volume_voxels * voxel_volume
        lesion_volume_tiv_normalized_mm3 = lesion_volume_mm3 / tiv_volume if tiv_volume > 0 else 0

        # Skip very small lesions (< 4 voxels) to get convex full image.
        if lesion_volume_voxels < 4: continue
        # Skip 2D lesions (flat along one axis): extent of the bounding box == 1.
        if any(s.stop - s.start == 1 for s in sl_box): continue

        idx_c = idx_sorted[lo:hi]
        lesion_coords = np.stack(np.unravel_index(idx_c, shape), axis=1)   # == np.argwhere(label == L)

        try:
            region_prop_temp = regionprops[lesion_label - 1]
        except Exception as e:
            print(f"Warning: regionprops failed for lesion label {lesion_label} with error: {e}")
            continue

        # Surface area (count exposed faces) on the bbox + 1-voxel zero margin
        sub = np.pad(lab_sub_of(lesion_labels, sl_box) == lesion_label, 1)
        surface_area = _mesh_surface_area(sub, voxel_sizes=voxel_sizes)
        # Thickness: distance (mm) from the lesion's deepest voxel to the nearest voxel outside it.
        max_inscribed_radius_mm = float(ndimage.distance_transform_edt(sub, sampling=voxel_sizes).max())

        # Coalescence (after dilation)
        new_label = orig_to_dil[lesion_label]
        merged_after_dilation = False
        merged_after_dilation_volume = 0
        n_lesions_merged = 1
        if new_label != 0:
            n_lesions_merged = int(lesions_per_dil[new_label])
            merged_after_dilation = n_lesions_merged > 1
            merged_after_dilation_volume = np.int64(dil_comp_size[new_label]) * voxel_volume

        # Erosion
        n_eroded_vox = ero_voxel_count[lesion_label]
        volume_after_erosion_mm3 = 0.0
        if n_eroded_vox == 0:
            lesion_after_erosion = "Vanished"
        elif n_eroded_vox == 1:
            lesion_after_erosion = "Remained Single"
            volume_after_erosion_mm3 = np.int64(ero_comp_size[single_ero_label[lesion_label]]) * voxel_volume
        else:
            lesion_after_erosion = "Split into smaller lesions"
            volume_after_erosion_mm3 = np.int64(max_ero_size[lesion_label]) * voxel_volume

        # Coalescence score
        coalescence_score = _calculate_coalescence_score(
            lesion_volume_tiv_normalized_mm3 * 100.0, volume_after_erosion_mm3, lesion_volume_mm3, lesion_after_erosion)

        # Fractal dimension (from voxel coordinates; box edges from the full shape)
        fractal_dimension = _calculate_fractal_dimension(lesion_coords, shape)

        # PCA eigenvalues (compute on mm-scaled coordinates so eigenvalues are in mm^2)
        eigenvalues = _compute_pca_eigenvalues(lesion_coords, voxel_sizes)

        # Tract match: tract with the largest overlap (first tract wins ties,
        # like max() over the insertion-ordered dict did); 'No Tract' if none.
        overlaps = counts[:, lesion_label]
        ranked = np.argsort(-overlaps, kind="stable")[:3]       # largest overlap first; the first tract wins ties
        matches = [
            (tract_names[int(i)], float(overlaps[i]) / float(lesion_volume_voxels)) if overlaps[i] > 0 else ("No Tract", 0.0)
            for i in ranked
        ]
        matches += [("No Tract", 0.0)] * (3 - len(matches))
        most_common_tract_name = matches[0][0]

        overall = {
            # --- Location / identity metrics ---
            "lesion_label": lesion_label,
            "centre_of_mass_x_mm": round(region_prop_temp.centroid[0] * voxel_sizes[0], 2),
            "centre_of_mass_y_mm": round(region_prop_temp.centroid[1] * voxel_sizes[1], 2),
            "centre_of_mass_z_mm": round(region_prop_temp.centroid[2] * voxel_sizes[2], 2),

            # --- Tract / anatomical matching ---
            "tract_match": most_common_tract_name,
            "tract_match_overlap_fraction": matches[0][1],      # share of the lesion's voxels inside that tract
            "tract_match_2": matches[1][0],
            "tract_match_2_overlap_fraction": matches[1][1],
            "tract_match_3": matches[2][0],
            "tract_match_3_overlap_fraction": matches[2][1],

            # --- Volume & coalescence metrics (physical units: mm^3) ---
            "volume_mm3": lesion_volume_mm3,
            "volume_tiv_normalized_mm3": lesion_volume_tiv_normalized_mm3,
            "volume_merged_after_dilation_mm3": merged_after_dilation_volume if merged_after_dilation else 0,
            "volume_after_erosion_mm3": volume_after_erosion_mm3,
            "coalescence_score": coalescence_score,

            # --- Merge flag ---
            "merged_after_dilation": merged_after_dilation,
            "n_lesions_merged_after_dilation": n_lesions_merged,

            # --- Shape / geometry metrics ---
            "volume_convex_hull_mm3": getattr(region_prop_temp, "convex_area", np.nan) * voxel_volume,
            "convex_deficiency_mm3": (getattr(region_prop_temp, "convex_area", 0) - getattr(region_prop_temp, "area", 0)) * voxel_volume,
            "equivalent_diameter_mm": ( (6 * lesion_volume_mm3 / np.pi) ** (1/3) ) if lesion_volume_mm3 > 0 else np.nan,
            "surface_area": surface_area,
            "max_inscribed_radius_mm": max_inscribed_radius_mm,
            "compactness": surface_area ** 3 / lesion_volume_mm3 ** 2 if lesion_volume_mm3 > 0 else np.nan,
            "sphericity": np.pi ** (1 / 3) * (6 * lesion_volume_mm3) ** (2 / 3) / surface_area if surface_area > 0 else np.nan,
            "euler_number": getattr(region_prop_temp, "euler_number", np.nan),
            "solidity": getattr(region_prop_temp, "solidity", np.nan),
            "extent": getattr(region_prop_temp, "extent", np.nan),
            "fractal_dimension": fractal_dimension,

            # --- PCA / orientation metrics ---
            "eigenvalue_1": eigenvalues[0],
            "eigenvalue_2": eigenvalues[1],
            "eigenvalue_3": eigenvalues[2],
            "ratio_first_second_eigenvalue": eigenvalues[0] / eigenvalues[1] if eigenvalues[1] else np.nan,
            "ratio_second_third_eigenvalue": eigenvalues[1] / eigenvalues[2] if eigenvalues[2] else np.nan,
            "ratio_first_third_eigenvalue": eigenvalues[0] / eigenvalues[2] if eigenvalues[2] else np.nan,
        }

        # Lesion-Map metrics: stats of each input map inside the lesion. The
        # lesion's voxels were gathered once (idx_c, C-order, same order as the
        # legacy boolean mask), so each map costs one small fancy-index.
        image_metrics = {}
        coord_tuple = tuple(lesion_coords.T)
        for image_metric, image_data in images_dict.items():
            image_metrics[image_metric] = _compute_image_stats(image_data[coord_tuple])

        lesion_dict[lesion_label] = {"overall": overall}
        lesion_dict[lesion_label].update(image_metrics)

    # Flatten to DataFrame. The legacy code created each "<map>_<stat>" column
    # with one DataFrame.loc[mask, col] = value per (lesion, map, stat); the
    # resulting columns appear in (map, stat) order and take the dtype of the
    # assigned numpy scalar (float32 for float32 maps' mean/min/max, float64
    # for percentiles ...). Building them from the same numpy scalars keeps
    # both the order and the dtypes -- hence the CSV text -- identical.
    lesion_df_list = []
    for lesion_label, lesion_data in lesion_dict.items():
        row = dict(lesion_data["overall"])
        row["lesion_label"] = lesion_label
        lesion_df_list.append(row)

    lesion_df = pd.DataFrame(lesion_df_list)

    # add subject column if provided
    lesion_df["subject"] = subject

    if len(lesion_dict):
        extra_cols = {}
        for lesion_label, lesion_data in lesion_dict.items():
            for metric, metric_data in lesion_data.items():
                if metric == "overall":
                    continue
                for stat, value in metric_data.items():
                    extra_cols.setdefault(f"{metric}_{stat}", []).append(value)
        # Column dtype rule of the legacy `.loc[mask, new_col] = scalar`: a numpy
        # floating scalar keeps its dtype (float32 stays float32), anything
        # else (integer min/max of an int map, python floats) becomes float64.
        def _col_dtype(vals):
            v0 = np.asarray(vals[0])
            return v0.dtype if np.issubdtype(v0.dtype, np.floating) else np.float64
        extra_df = pd.DataFrame(
            {c: pd.Series(v, index=lesion_df.index, dtype=_col_dtype(v)) for c, v in extra_cols.items()}
        )
        lesion_df = pd.concat([lesion_df, extra_df], axis=1)

    # Reorder columns
    cols = lesion_df.columns.tolist()
    main_cols = ["subject", "lesion_label", "volume_mm3", "tract_match"]
    cols = main_cols + [x for x in cols if x not in main_cols]
    lesion_df = lesion_df[cols]

    # Aggregate per-tract: median of numeric metrics grouped by tract_match (include 'No Tract')
    if lesion_df.shape[0] == 0:
        tract_aggregated_lesion_df = pd.DataFrame(columns=["subject", "tract_name"])
    else:
        numeric_cols = [c for c in lesion_df.select_dtypes(include=[np.number]).columns if c != "lesion_label"]
        # compute medians while explicitly skipping NaNs
        grouped = (
            lesion_df.groupby("tract_match")[numeric_cols]
            .agg(lambda s: s.median(skipna=True))
            .reset_index()
        )
        # rename for clarity
        grouped = grouped.rename(columns={"tract_match": "tract"})
        # add subject column
        grouped.insert(0, "subject", subject)
        # Drop rows with all Nan values
        grouped = grouped.dropna(how="all", axis=0, subset=numeric_cols)
        # Drop columns with all values the same, but always keep 'subject' and 'tract'
        cols_to_keep = [col for col in grouped.columns if col in ("subject", "tract") or (grouped[col] != grouped[col].iloc[0]).any()]
        grouped = grouped[cols_to_keep]

        tract_cols = ["subject", "tract"] + [c for c in grouped.columns if c not in ("subject", "tract")]
        tract_aggregated_lesion_df = grouped[tract_cols]

    return lesion_df, tract_aggregated_lesion_df


def lab_sub_of(labels, sl_box):
    """Bounding-box view of the label image for one lesion."""
    return labels[sl_box]
