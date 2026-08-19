# ============================================================================
# Metric family: per-lesion metrics 
# ============================================================================
import os
import scipy
from scipy import ndimage
import skimage
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import nibabel as nib


def _count_boxes_multiscales(seg, scales):

    indices = np.where(seg > 0)

    Lx = seg.shape[0]
    Ly = seg.shape[1]
    Lz = seg.shape[2]

    pixels = np.asarray(indices).T

    Ns = []
    # looping over several scales
    for scale in scales:
        # computing the histogram
        H, edges = np.histogramdd(indices, bins=(np.arange(0, Lx, scale),
                                                    np.arange(0, Ly, scale),
                                                    np.arange(0, Lz, scale)))
        Ns.append(np.sum(H > 0))

    return Ns

def _calculate_fractal_dimension(lesion_mask):
    if np.sum(lesion_mask) < 6:
        return 0

    epsilon_max = 3
    epsilon_min = 0

    scales = np.logspace(epsilon_min, epsilon_max, num=15, endpoint=True, base=2)
    scales = np.array([a for a in scales if a != 0.0])

    offset_res = []
    for x_pad in range(0, 6):
        for y_pad in range(0, 6):
            for z_pad in range(0, 6):
                seg_pad = np.pad(lesion_mask, ((x_pad, x_pad), 
                                                (y_pad,y_pad),
                                                (z_pad, z_pad)), 
                                                mode='constant')
                Ns_temp = _count_boxes_multiscales(seg_pad, scales)
                offset_res.append(np.expand_dims(Ns_temp, 1))

    final_concat = np.concatenate(offset_res, 1)
    min_boxes = np.min(final_concat, 1)

    # linear fit, polynomial of degree 1
    log_scales = np.log(scales)
    log_min_boxes = np.log(min_boxes)

    # Remove pairs with nan values
    indices = np.where(np.isnan(log_scales) | np.isnan(log_min_boxes))
    log_scales = np.delete(log_scales, indices)
    log_min_boxes = np.delete(log_min_boxes, indices)

    coeffs = np.polyfit(np.log(scales), np.log(min_boxes), 1)

    # The Hausdorff dimension is -coeffs[0]
    # the fractal dimension is the OPPOSITE of the fitting coefficient
    return -1 * coeffs[0]

def _calculate_coalescence_score(lesion_volume_tiv_normalized_mm3, volume_after_erosion, lesion_volume_mm3, lesion_after_erosion):
    coalescence_score = None
    if lesion_after_erosion == "Vanished":
        coalescence_score = 0
    if lesion_volume_tiv_normalized_mm3 < 0.05:
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

def _calculate_surface_area(lesion_mask, voxel_sizes=(1.0, 1.0, 1.0)):
    x_dim, y_dim, z_dim = lesion_mask.shape
    dx, dy, dz = voxel_sizes
    surface_area = 0.0
    for x in range(x_dim):
        for y in range(y_dim):
            for z in range(z_dim):
                if lesion_mask[x, y, z] == 1:
                    # x-direction faces
                    for nx in [x-1, x+1]:
                        if nx < 0 or nx >= x_dim or lesion_mask[nx, y, z] == 0:
                            surface_area += dy * dz
                    # y-direction faces
                    for ny in [y-1, y+1]:
                        if ny < 0 or ny >= y_dim or lesion_mask[x, ny, z] == 0:
                            surface_area += dx * dz
                    # z-direction faces
                    for nz in [z-1, z+1]:
                        if nz < 0 or nz >= z_dim or lesion_mask[x, y, nz] == 0:
                            surface_area += dx * dy
    return surface_area

def _map_label_relationships(original_labels, transformed_labels, allow_multiple=False):
    """Build mappings between original and transformed label images.

    Returns (label_dict, orig_to_new) where label_dict maps transformed_label -> [original_labels]
    and orig_to_new maps original_label -> transformed_label (or list of transformed labels
    when allow_multiple is True).
    """
    label_dict = defaultdict(list)
    for orig, trans in zip(original_labels.ravel(), transformed_labels.ravel()):
        if orig != 0:
            label_dict[trans].append(orig)

    if allow_multiple:
        orig_to_new = defaultdict(list)
        for new_label, originals in label_dict.items():
            if new_label != 0:
                for o in originals:
                    orig_to_new[o].append(new_label)
    else:
        orig_to_new = {}
        for new_label, originals in label_dict.items():
            for o in originals:
                orig_to_new[o] = new_label

    return label_dict, orig_to_new

def _compute_image_stats(values: np.ndarray) -> dict:
    """Compute summary statistics for a 1D array of values, returning a dict.
    Handles empty arrays by returning NaNs.
    """
    return {
        "mean": np.mean(values) if len(values) > 0 else np.nan,
        "median": np.median(values) if len(values) > 0 else np.nan,
        "std": np.std(values) if len(values) > 0 else np.nan,
        "min": np.min(values) if len(values) > 0 else np.nan,
        "max": np.max(values) if len(values) > 0 else np.nan,
        "25th": np.percentile(values, 25) if len(values) > 0 else np.nan,
        "75th": np.percentile(values, 75) if len(values) > 0 else np.nan,
        "IQR": np.percentile(values, 75) - np.percentile(values, 25) if len(values) > 0 else np.nan,
        "peak_width": np.percentile(values, 95) - np.percentile(values, 5) if len(values) > 0 else np.nan,
        "kurtosis": scipy.stats.kurtosis(values) if len(values) > 0 else np.nan,
        "skew": scipy.stats.skew(values) if len(values) > 0 else np.nan,
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

def _compute_lesion_erosion_effect(orig_to_new_dict_eroded, lesion_label, lesion_labels_eroded, voxel_volume=1.0):
    eroded_labels = orig_to_new_dict_eroded[lesion_label] if lesion_label in orig_to_new_dict_eroded else []
    volume_after_erosion_mm3 = 0.0
    if len(eroded_labels) == 0:
        lesion_after_erosion = "Vanished"
    elif len(eroded_labels) == 1:
        lesion_after_erosion = "Remained Single"
        volume_after_erosion_mm3 = np.sum(np.where(lesion_labels_eroded == orig_to_new_dict_eroded[lesion_label][0], 1, 0)) * voxel_volume
    else:
        lesion_after_erosion = "Split into smaller lesions"
        volume_after_erosion_mm3 = np.max([np.sum(np.where(lesion_labels_eroded == x, 1, 0)) for x in orig_to_new_dict_eroded[lesion_label]]) * voxel_volume

    return lesion_after_erosion, volume_after_erosion_mm3
    
def _compute_lesion_metrics(
    subject: str,
    lesion_array: np.ndarray,
    tiv_volume: float,
    tract_binary_arrays: dict,
    images_dict: dict,
    voxel_sizes: tuple,
    ):
    """
    Compute per-lesion metrics and write to lesion_metrics.csv.
    lesion_array: 3D (binary) lesion mask
    tiv_volume: TIV volume (float)
    tract_binary_arrays: dict mapping tract name -> binarised tract mask (bool/uint8).
                         Callers passing these same arrays elsewhere should reuse this
                         exact dict so tract_match stays consistent across outputs.
    images_dict: dict of {metric_name: 3D array}
    voxel_sizes: tuple of (x_size, y_size, z_size) in mm
    """
    # Lesion voxel volume in mm^3
    voxel_volume = np.prod(voxel_sizes)

    # Binarize lesion array in case it's not already binary
    lesion_array = (lesion_array > 0).astype(np.uint8)

    # Label connected lesions
    lesion_labels, lesion_labels_count = ndimage.label(lesion_array, structure=np.ones((3, 3, 3)))
    # Save the image of labeled lesions for potential QC use (not currently saved to disk, but could be returned or saved as needed)
    regionprops = skimage.measure.regionprops(lesion_labels)

    # Dilate and erode for coalescence
    lesion_mask_dilated = ndimage.binary_dilation(lesion_array, iterations=1, structure=np.ones((3, 3, 3)))
    lesion_labels_dilated, _ = ndimage.label(lesion_mask_dilated, structure=np.ones((3, 3, 3))) 
    lesion_mask_eroded = ndimage.binary_erosion(lesion_array, iterations=1, structure=np.ones((3, 3, 3)))
    lesion_labels_eroded, _ = ndimage.label(lesion_mask_eroded, structure=np.ones((3, 3, 3)))

    # Map original to new labels after dilation/erosion
    label_dict_dilated, orig_to_new_dict_dilated = _map_label_relationships(
        lesion_labels, lesion_labels_dilated, allow_multiple=False
    )
    label_dict_eroded, orig_to_new_dict_eroded = _map_label_relationships(
        lesion_labels, lesion_labels_eroded, allow_multiple=True
    )

    lesion_dict = {}

    for lesion_label in range(1, lesion_labels_count + 1):

        lesion_mask_temp = (lesion_labels == lesion_label).astype(np.uint8)
        lesion_coords = np.argwhere(lesion_labels == lesion_label)

        lesion_volume_voxels = np.sum(lesion_mask_temp)
        lesion_volume_mm3 = lesion_volume_voxels * voxel_volume
        lesion_volume_tiv_normalized_mm3 = lesion_volume_mm3 / tiv_volume if tiv_volume > 0 else 0

        # Skip very small lesions (< 4 voxels) to get convex full image.
        if lesion_volume_voxels < 4: continue
        # Check if lesion is 2D (if all coords have same value in one dimension), if so skip to avoid regionprops errors.
        if np.any(np.ptp(lesion_coords, axis=0) == 0): continue

        try:
            region_prop_temp = regionprops[lesion_label - 1]
        except Exception as e:
            print(f"Warning: regionprops failed for lesion label {lesion_label} with error: {e}")
            continue

        # Surface area (count exposed faces)
        surface_area = _calculate_surface_area(lesion_mask_temp, voxel_sizes=voxel_sizes)

        # Coalescence (after dilation)
        new_label = orig_to_new_dict_dilated.get(lesion_label, None)
        merged_after_dilation = False
        merged_after_dilation_volume = 0
        if new_label is not None:
            merged_after_dilation = len(label_dict_dilated[new_label]) > 1
            merged_after_dilation_volume = np.sum(lesion_labels_dilated == new_label) * voxel_volume

        # Erosion
        lesion_after_erosion, volume_after_erosion_mm3 = _compute_lesion_erosion_effect(orig_to_new_dict_eroded, lesion_label, lesion_labels_eroded, voxel_volume)

        # Coalescence score
        coalescence_score = _calculate_coalescence_score(lesion_volume_tiv_normalized_mm3, volume_after_erosion_mm3, lesion_volume_mm3, lesion_after_erosion)

        # Fractal dimension
        fractal_dimension = _calculate_fractal_dimension(lesion_mask_temp)

        # PCA eigenvalues (compute on mm-scaled coordinates so eigenvalues are in mm^2)
        eigenvalues = _compute_pca_eigenvalues(lesion_coords, voxel_sizes)

        # Tract match
        tract_lesion_overlap_dict = {}
        for tract_name, tract_binary in tract_binary_arrays.items():
            if tract_binary.shape != lesion_array.shape:
                raise ValueError(f"Tract array for {tract_name} has shape {tract_binary.shape} but expected {lesion_array.shape}")

            overlap_vox = int(np.sum(lesion_mask_temp.astype(bool) & (tract_binary > 0)))
            overlap_volume_mm3 = overlap_vox * voxel_volume
            if overlap_volume_mm3 > 0:
                tract_lesion_overlap_dict[tract_name] = overlap_volume_mm3

        # Get tract with largest overlap, or 'No Tract' if no overlap
        if tract_lesion_overlap_dict:
            most_common_tract_name = max(tract_lesion_overlap_dict, key=lambda k: tract_lesion_overlap_dict[k])
        else:
            most_common_tract_name = "No Tract"
        

        # Build a single flat `overall` dict with commented category headers
        overall = {
            # --- Location / identity metrics ---
            "lesion_label": lesion_label,
            "centre_of_mass_x_mm": round(region_prop_temp.centroid[0] * voxel_sizes[0], 2),
            "centre_of_mass_y_mm": round(region_prop_temp.centroid[1] * voxel_sizes[1], 2),
            "centre_of_mass_z_mm": round(region_prop_temp.centroid[2] * voxel_sizes[2], 2),

            # --- Tract / anatomical matching ---
            "tract_match": most_common_tract_name,

            # --- Volume & coalescence metrics (physical units: mm^3) ---
            "volume_mm3": lesion_volume_mm3,
            "volume_tiv_normalized_mm3": lesion_volume_tiv_normalized_mm3,
            "volume_merged_after_dilation_mm3": merged_after_dilation_volume if merged_after_dilation else 0,
            "volume_after_erosion_mm3": volume_after_erosion_mm3,
            "coalescence_score": coalescence_score,

            # --- Merge flag ---
            "merged_after_dilation": merged_after_dilation,

            # --- Shape / geometry metrics ---
            "volume_convex_hull_mm3": getattr(region_prop_temp, "convex_area", np.nan) * voxel_volume,
            "convex_deficiency_mm3": (getattr(region_prop_temp, "convex_area", 0) - getattr(region_prop_temp, "area", 0)) * voxel_volume,
            "equivalent_diameter_mm": ( (6 * lesion_volume_mm3 / np.pi) ** (1/3) ) if lesion_volume_mm3 > 0 else np.nan,
            "surface_area": surface_area,
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

        # Lesion-Map metrics: compute stats for each input map (e.g. FA, MD) within the lesion
        image_metrics = {}
        for image_metric, image_data in images_dict.items():
            lesion_values = image_data[lesion_mask_temp == 1]
            image_metrics[image_metric] = _compute_image_stats(lesion_values)

        lesion_dict[lesion_label] = {"overall": overall}
        lesion_dict[lesion_label].update(image_metrics)

    # Flatten to DataFrame
    lesion_df_list = []
    for lesion_label, lesion_data in lesion_dict.items():
        lesion_data["overall"]["lesion_label"] = lesion_label
        lesion_df_list.append(lesion_data["overall"])

    lesion_df = pd.DataFrame(lesion_df_list)

    # add subject column if provided
    lesion_df["subject"] = subject

    for lesion_label, lesion_data in lesion_dict.items():
        for metric, metric_data in lesion_data.items():
            if metric == "overall":
                continue
            for stat, value in metric_data.items():
                lesion_df.loc[lesion_df["lesion_label"] == lesion_label, f"{metric}_{stat}"] = value

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