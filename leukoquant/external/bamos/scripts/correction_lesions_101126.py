# Modified by Stylianos Charalampous on 11/11/2025
# Modified by Stylianos Charalampous on 10/11/2026: each correction is now an independent object selected by
# the -corr names, and the choroid correction runs when selected (it never ran before).

import nibabel as nib
import pandas as pd
import os
import argparse
import numpy as np
import sys
from dataclasses import dataclass
from functools import cached_property
from typing import Callable, Optional
from Parsing_LesionFile_args import parse
from scipy.ndimage import binary_dilation
from scipy.ndimage import label

def create_ventricles(parcellation):
    ventr = np.where(np.logical_or(parcellation==52, parcellation==53),
                     np.ones_like(parcellation), np.zeros_like(parcellation))
    return ventr

def create_border_thalamusventr(parcellation):
    ventr = np.where(np.logical_or(parcellation == 52, parcellation == 53),
                     np.ones_like(parcellation), np.zeros_like(parcellation))
    thal = np.where(np.logical_or(parcellation == 60, parcellation == 61),
                     np.ones_like(parcellation), np.zeros_like(parcellation))
    border_thal = binary_dilation(ventr, iterations=1) * thal
    return border_thal


def slicewise_choroid(lesion, ventricles, border_thal):
    intersect = lesion * (ventricles+border_thal)

    indices = np.asarray(np.where(intersect>0)).T
    slices = np.unique(indices[:,2])
    choroid = np.zeros_like(lesion)
    choroid_temp = np.zeros_like(lesion)
    for sl in slices:
        full_slice = intersect[..., sl]
        bin_slice = np.asarray(full_slice, dtype=bool)
        labelled, numb_labels = label(bin_slice)
        for lab in range(1, numb_labels+1):
            lab_temp = np.where(labelled == lab, np.ones_like(labelled),
                                np.zeros_like(labelled))
            dil_bord = binary_dilation(lab_temp, iterations=1) - lab_temp
            ventrinter = dil_bord * ventricles[..., sl]
            ratio = np.sum(ventrinter)*1.0/np.sum(dil_bord)
            # print(ratio, sl, lab, np.sum(lab_temp))
            if ratio > 0.75:
                choroid[..., sl] = np.where(lab_temp, lesion[...,sl],
                                            choroid[...,sl])
            if ratio >= 0.25:
                choroid_temp[..., sl] = np.where(lab_temp, lesion[..., sl],
                                            choroid_temp[..., sl])
    aggreg_choroid = np.squeeze(np.sum(choroid, axis=2))
    bin_aggreg = np.where(aggreg_choroid>0, np.ones_like(aggreg_choroid),
                          np.zeros_like(aggreg_choroid))
    dilated_aggreg = binary_dilation(bin_aggreg, iterations=1)
    # print(np.sum(choroid), np.sum(choroid_temp))
    aggreg_3d = np.zeros_like(choroid_temp)
    for sl in slices:
        aggreg_3d[..., sl] = aggreg_choroid
        if np.sum(choroid_temp[..., sl]) > 0:
            bin_slice = np.where(choroid_temp[..., sl]>0, np.ones_like(
                aggreg_choroid), np.zeros_like(aggreg_choroid))
            labelled, numb_labels = label(bin_slice)
            for lab in range(1, numb_labels+1):
                temp_lab = np.where(labelled == lab, np.ones_like(bin_slice),
                                    np.zeros_like(bin_slice))

                if np.sum(dilated_aggreg * temp_lab) > 0:
                    choroid[...,sl] = np.where(temp_lab>0, choroid_temp[...,
                                                                        sl],
                                               choroid[...,sl])
    print('Choroid correction is' ,np.sum(choroid))
    return choroid, aggreg_3d

def ventr_corr(lesion, parcellation, dil=1):
    ventr = create_ventricles(parcellation)
    #struct = np.zeros((3,3,3))

    dil_ventr = binary_dilation(ventr,iterations=dil)
    corr_ventr  = lesion*dil_ventr
    print("Ventr correction is ", np.sum(corr_ventr))
    return corr_ventr

def cortical_corr(file, lesion, connected):
    df_parse = pd.DataFrame.from_dict(parse(file, 2))
    to_remove = np.zeros_like(lesion)
    if df_parse.shape[0] == 0:
        return to_remove
    df_toremove = df_parse[df_parse['ot']=='16']
    for lab in df_toremove['label']:
        label = int(lab)
        to_remove = np.where(connected == label, lesion, to_remove)
    print("Cortical 16 volume is ", np.sum(to_remove))
    return to_remove


def hypot1_corr(file, lesion, connected):
    df_parse = pd.DataFrame.from_dict(parse(file,2))
    to_remove = np.zeros_like(lesion)
    if df_parse.shape[0] == 0:
        return to_remove
    df_toremove = df_parse[df_parse['mahalmean1'].astype(float)<-5.5]
    for lab in df_toremove['label']:
        label = int(lab)
        to_remove = np.where(connected==label, lesion, to_remove)
    print("To remove because of T1 hypo is ", np.sum(to_remove))
    return to_remove


def cortical_sheet(lesion, parcellation, thresh_vol=10, ratio_thresh=0.55):
    shape = lesion.shape
    gm = (parcellation>98)
    correction_axial = np.zeros_like(lesion)
    structure = np.ones([3,3])
    for s in range(0, shape[2]):

        lesion_temp = lesion[:,:,s]
        correct_temp = np.zeros_like(lesion_temp)
        gm_temp = gm[:,:,s]
        bin_temp = (lesion_temp > 0.5)
        if np.sum(bin_temp) > thresh_vol:
            conn_temp, num_lab = label(bin_temp, structure=structure)
            for val in range(1, num_lab+1):
                seg_temp = (conn_temp==val).astype(float)
                value_temp = np.where(seg_temp, lesion_temp, np.zeros_like(
                    lesion_temp))
                if np.sum(seg_temp) > thresh_vol/2:
                    border_ext = binary_dilation(seg_temp) - seg_temp
                    gm_border = border_ext * gm_temp
                    ratio_gmborder = np.sum(gm_border) / np.sum(seg_temp)
                    if ratio_gmborder > ratio_thresh:
                        correct_temp += value_temp
            correction_axial[:,:,s] = correct_temp
    correction_coronal = np.zeros_like(lesion)
    for s in range(0, shape[1]):
        lesion_temp = lesion[:, s, :]
        correct_temp = np.zeros_like(lesion_temp)
        gm_temp = gm[:, s, :]
        bin_temp = (lesion_temp > 0.5)
        if np.sum(bin_temp) > thresh_vol:
            conn_temp, num_lab = label(bin_temp, structure=structure)
            for val in range(1, num_lab + 1):
                seg_temp = (conn_temp == val).astype(float)
                value_temp = np.where(seg_temp, lesion_temp, np.zeros_like(
                    lesion_temp))
                if np.sum(seg_temp) > thresh_vol/2:
                    border_ext = binary_dilation(seg_temp) - seg_temp
                    gm_border = border_ext * gm_temp
                    ratio_gmborder = np.sum(gm_border) / np.sum(seg_temp)
                    if ratio_gmborder > ratio_thresh:
                        correct_temp += value_temp
            correction_coronal[:, s, :] = correct_temp

    correction_sagittal = np.zeros_like(lesion)
    for s in range(0, shape[0]):
        lesion_temp = lesion[s, :, :]
        correct_temp = np.zeros_like(lesion_temp)
        gm_temp = gm[s, :, :]
        bin_temp = (lesion_temp > 0.5)
        if np.sum(bin_temp) > thresh_vol:
            conn_temp, num_lab = label(bin_temp, structure=structure)
            for val in range(1, num_lab + 1):
                seg_temp = (conn_temp == val).astype(float)
                value_temp = np.where(seg_temp, lesion_temp, np.zeros_like(
                    lesion_temp))
                if np.sum(seg_temp) > thresh_vol/2:
                    border_ext = binary_dilation(seg_temp) - seg_temp
                    gm_border = border_ext * gm_temp
                    ratio_gmborder = np.sum(gm_border) / np.sum(seg_temp)
                    if ratio_gmborder > ratio_thresh:
                        correct_temp += value_temp
            correction_sagittal[s, :, :] = correct_temp
    correction_fin = correction_coronal + correction_axial + correction_sagittal
    print("Sagittal, Axial, Coronal correction:%f %f %f" % (np.sum(
        correction_sagittal), np.sum(correction_axial), np.sum(
        correction_coronal)))
    correction_fin = np.minimum(lesion, correction_fin)
    return correction_fin




# ---------------------------------------------------------------------------------------------
# Corrections. Each one is independent: it reads the original lesion map, writes its own
# intermediate NIfTI (deleted at the end of main) and returns the voxel weights to subtract.
# A correction is selected when any of its option names is passed to -corr.
# ---------------------------------------------------------------------------------------------

@dataclass
class CorrectionInputs:
    """Everything a correction may need; files are read only when a correction asks for them."""
    lesion_data: np.ndarray
    affine: np.ndarray
    path_results: str
    parc_path: str
    connect_data: Optional[np.ndarray]
    label_path: str
    dil_ventr: int

    @cached_property
    def parc_data(self):
        return nib.load(self.parc_path).get_fdata()

    def save(self, data, filename):
        nib.save(nib.Nifti1Image(data, self.affine), os.path.join(self.path_results, filename))

    def require_connectivity(self, purpose):
        if self.connect_data is None:
            raise ValueError("Connectivity file must be provided for %s correction" % purpose)
        return self.connect_data


def _choroid(inp):
    ventricles = create_ventricles(inp.parc_data)
    border_thal = create_border_thalamusventr(inp.parc_data)
    choroid, aggreg_3d = slicewise_choroid(inp.lesion_data, ventricles, border_thal)
    inp.save(choroid, 'ChoroidCorr.nii.gz')
    inp.save(aggreg_3d, 'AggregChor.nii.gz')
    return choroid


def _cortical(inp):
    cortical = cortical_corr(inp.label_path, inp.lesion_data, inp.require_connectivity("cortical"))
    inp.save(cortical, 'CorticalCorr.nii.gz')
    return cortical


def _sheet(inp):
    sheet = cortical_sheet(inp.lesion_data, inp.parc_data)
    inp.save(sheet, 'SheetCorr.nii.gz')
    return sheet


def _hypot1(inp):
    hypo = hypot1_corr(inp.label_path, inp.lesion_data, inp.require_connectivity("hypo T1"))
    inp.save(hypo, 'HypoCorr.nii.gz')
    return hypo


def _ventr(inp):
    ventr = ventr_corr(inp.lesion_data, inp.parc_data, dil=inp.dil_ventr)
    inp.save(ventr, 'VentrCorr.nii.gz')
    return ventr


@dataclass(frozen=True)
class Correction:
    log_name: str        # name used in the "Lesion range after ..." log line
    options: frozenset   # -corr values that select this correction
    compute: Callable    # CorrectionInputs -> array to subtract

    def is_selected(self, requested):
        return not self.options.isdisjoint(requested)


# "all" means choroid, cortical and sheet (not hypot1 or ventr).
CORRECTIONS = [
    Correction('choroid', frozenset({'choroid', 'all', 'choroid-sheet', 'choroid-cortex'}),
               _choroid),
    Correction('cortical', frozenset({'all', 'cortex', 'choroid-cortex', 'cortical-sheet'}), _cortical),
    Correction('sheet', frozenset({'all', 'sheet', 'choroid-sheet', 'cortical-sheet'}), _sheet),
    Correction('hypo corr', frozenset({'hypot1'}), _hypot1),
    Correction('ventr', frozenset({'ventr'}), _ventr),
]


def main(argv):

    parser = argparse.ArgumentParser(description='Correction of lesions')
    parser.add_argument('-les', dest='lesion', metavar='lesion file',
                        type=str, required=True,
                        help='Path to lesion file')
    parser.add_argument('-parc', dest='parc', action='store',
                        default='', type=str,
                        help='path to parcellation file')
    parser.add_argument('-connect', dest='connect', action='store',
                        default='',
                        type=str)
    parser.add_argument('-dil', dest='dil_ventr', action='store', type=int,
                       default=1, help='Number of dilation iterations for '
                                       'ventricle correction')
    parser.add_argument('-label', dest='label', action='store', type=str,
                        default='')
    parser.add_argument('-corr',nargs='+', dest='corr', action='store',type=str,
                        default='all',choices=['choroid', 'cortex',
                                                'sheet', 'choroid-cortex',
                                                'choroid-sheet',
                                                'cortical-sheet',
                                                'all','hypot1','ventr'])

    try:
        args = parser.parse_args()
    except argparse.ArgumentTypeError:
        print('correction_lesions_101126.py \
            -les <lesion file> \
            -parc <parcellation file> \
            -connect <connectivity file> \
            -dil <dilation iterations> \
            -label <label file> \
            -corr <correction type> \
            -id <subject id>')

        sys.exit(2)


    lesion_nii = nib.load(args.lesion)
    # Read the stored voxel values WITHOUT applying the NIfTI scl_slope/scl_inter.
    # BaMoS lesion maps hold weights in [0, 1], but BaMoS builds up to this release copied the
    # input scan's slope into the lesion map header; get_fdata() applies that slope and, because
    # the result is saved below with a fresh header, would bake the wrong scaling permanently into
    # CorrectLesion.nii.gz (EPAD maps came out as 2.94-5.87 instead of 0.5-1).
    # nibabel moves the header scaling onto the array proxy (the header's own slope reads as None).
    slope, inter = lesion_nii.dataobj.slope, lesion_nii.dataobj.inter
    if (slope not in (None, 1.0)) or (inter not in (None, 0.0)):
        print('WARNING: lesion header has scl_slope=%s scl_inter=%s; ignoring it and using the '
              'raw stored values' % (slope, inter))
    lesion_data = np.asanyarray(lesion_nii.dataobj.get_unscaled()).astype(np.float32)
    if np.max(lesion_data) > 1.0 + 1e-3:
        print('WARNING: raw lesion values exceed 1 (max %s); the lesion weights are expected in '
              '[0, 1], so this map has a problem that is not the header scaling' % np.max(lesion_data))
    print('Lesion range before corr ',np.min(lesion_data), np.max(lesion_data))
    path_results = os.path.split(args.lesion)[0]
    connect_nii = None
    if args.connect != '':
        connect_nii = nib.load(args.connect)
    # Selected corrections are all computed on the ORIGINAL lesion map and subtracted afterwards.
    # The -corr option takes several names (e.g. "-corr choroid cortex sheet"); a bare string is
    # also accepted because that is the argparse default.
    requested = {args.corr} if isinstance(args.corr, str) else set(args.corr)
    inputs = CorrectionInputs(lesion_data=lesion_data, affine=lesion_nii.affine,
                              path_results=path_results, parc_path=args.parc,
                              connect_data=None if connect_nii is None else connect_nii.get_fdata(),
                              label_path=args.label, dil_ventr=args.dil_ventr)
    corrections = [c for c in CORRECTIONS if c.is_selected(requested)]
    results = [(c, c.compute(inputs)) for c in corrections]
    for correction, correction_data in results:
        lesion_data -= correction_data
        print('Lesion range after %s is' % correction.log_name, np.min(lesion_data), np.max(lesion_data))

    lesion_data = np.where(lesion_data < 0, 0.0, lesion_data)
    print('Lesion range after clip', np.min(lesion_data), np.max(lesion_data))

    nii_lesion = nib.Nifti1Image(lesion_data, lesion_nii.affine)
    nib.save(nii_lesion, os.path.join(path_results,
                                      'CorrectLesion.nii.gz'))
    
    # Delete intermediate files if they were created
    intermediate_files = ['ChoroidCorr.nii.gz', 'AggregChor.nii.gz', 'CorticalCorr.nii.gz', 'SheetCorr.nii.gz', 'HypoCorr.nii.gz', 'VentrCorr.nii.gz']
    for filename in intermediate_files:
        filepath = os.path.join(path_results, filename)
        if os.path.exists(filepath):
            os.remove(filepath)

if __name__ == '__main__':
    main(sys.argv[1:])