#!/bin/bash
# Modified by Stylianos Charalampous on 08/05/2026
# set -x

#### WARNING !!!! THIS SCRIPT ASSUMES THAT THE FLAIR IMAGE WILL BE REGISTERED TO THE T1 IMAGE. ALL RESULTS SHOULD BE VISUALISED IN THE T1 SPACE. ALL SCANS ARE PUT IN THE RAS ORIENTATION BY DEFAULT

ArrayModalities=("T1" "T2" "FLAIR" "PD" "SWI")
ListCorr=(102 103 187 188 47 48 49)
AW=0
JC=0
MRE=250
OptSP=0
OptCL=2
OptOL=0
OptWMI=2 # If >0, indicates that a higher level of sensitivity to lesion should be considered from voxels coming from WMI
TVS=1
arrayMod=("T1" "FLAIR")
arrayModNumber=(1 3)
# ArtefactArray=(101 102 139 140 187 188 167 168)
# ArtefactArray=(5 12 32 33 39 40 50 51 62 63 64 65 72 73 74 48 49 101 102 117 118 139 140 171 172 187 188 167 168 202 203 181 182 117 118 123 124 133 134 155 156 181 182 185 186 201 202 203 204 207 208)
# ListCorr=(5 12 32 33 39 40 50 51 62 63 64 65 72 73 74 48 49 101 102 117 118 139 140 171 172 187 188 167 168 202 203 181 182 117 118 123 124 133 134 155 156 181 182 185 186 201 202 203 204 207 208)
# #ArtefactArray=(101 102 139 140 187 188 167 168 39 40 117 118 123 124 125 126 133 134 137 138 141 142 147 148 155 156 181 182 185 186 187 188 201 202 203 204 207 208)

ListCorr=(101 102 103 104 105 106 187 188 47 48 49 32 33 173 174)
ArtefactArray=(5 12 24 31 39 40 50 51 62 63 64 65 72 73 74 48
	49 101 102 105 106 139 140 187 188 167 168 39 40 117 118 123 124 125 126 133 134 137 138 141 142 147 148 155 156 181 182 185 186 187 188 201 202 203 204 207 208)

#RuleFileName=/cluster/project0/SegBiASM/DataToTryBaMoS/GenericRule_CSF.txt
#NameGMatrix=/cluster/project0/SegBiASM/DataToTryBaMoS/GMatrix4_Low3.txt

# Current File directory of the script
current_path=$(dirname "$(readlink -f "$0")")

RuleFileName=${current_path}/GenericRule_CSF.txt
NameGMatrix=${current_path}/GMatrix4_Low3.txt

# RuleFileName=/mounts/auto/xnat/pipelines/BiASM/GenericRule_CSF.txt

# NameGMatrix=/mounts/auto/xnat/pipelines/BiASM/GMatrix4_Low3.txt
#Mem=7.9
OW=0.01
#OW=0.01
#OW=0.10
#OW=0.50
#OW=0.001
#OW=0.99
JC=1
Opt=${OptCross}
OptVL=1
OptSP=1
OptTA=1
OptTxt=Test
Opt=TA
# Flags for levels of correction
flag_NoCerrebellum=0
flag_furtherCorr=1
flag_corrDGM=1

PathReg=/leukoquant/leukoquant/external/niftyreg/bin
PathSeg=/leukoquant/leukoquant/external/bamos/seg-apps/bin
#PathFSL=$(which fslmaths | xargs dirname 2>/dev/null || echo "/opt/conda/envs/bamos/bin")
PathICBM=/leukoquant/leukoquant/external/bamos/ICBM_Priors

if [ $# -lt 5 ] || [ $# -gt 11 ]; then
	echo ""
	echo "*******************************************************************************"
	echo "Usage: $0 ID ImageFLAIR ImageT1 GIF_results_path PathResults jump_start opt space scratch_dir"
	echo "*******************************************************************************"
	echo "ID: Subject ID"
	echo "ImageFLAIR: Path to FLAIR image"
	echo "ImageT1: Path to T1 image"
	echo "GIF_results_path: Path to GIF results"
	echo "PathResults: Path to results directory"
	echo "jump_start: Whether to jump start the processing"
	echo "opt: Optimization method"
	echo "space: Space to process in - 1 for T1 space, 2 for FLAIR space"
	echo "scratch_dir: Path to scratch directory"
	echo "*******************************************************************************"
	echo ""
	exit
fi

ID=$1
ImageFLAIR=$2
ImageT1=$3
GIF_results_path=$4
PathResults=$5
JUMP_START=$6
Opt=$7
Space=$8
PathScratch=$9/BaMoS_${ID}

echo "ID: $ID"
echo "ImageFLAIR: $ImageFLAIR"
echo "ImageT1: $ImageT1"
echo "GIF_results_path: $GIF_results_path"
echo "PathResults: $PathResults"
echo "jump_start: $JUMP_START"
echo "opt: $Opt"
echo "space: $Space"
echo "scratch_dir: $PathScratch"

# Folders declaration
PN=${ID}
ModalitiesTot=T1FLAIR
PathGIF=${GIF_results_path}

# Binaries declaration
seg_maths_bin=${PathSeg}/seg_maths
seg_biasm_bin=${PathSeg}/Seg_BiASM
seg_analysis_bin=${PathSeg}/Seg_Analysis
reg_aladin_bin=${PathReg}/reg_aladin
reg_resample_bin=${PathReg}/reg_resample
reg_f3d_bin=${PathReg}/reg_f3d

# Files declarion
# Scratch files declaration
gif_parcellation_file=${PathScratch}/GIF_${ID}_Parcellation.nii.gz
gif_segmentation_file=${PathScratch}/GIF_${ID}_Segmentation.nii.gz
gif_gm_file=${PathScratch}/GIF_GM_${ID}.nii.gz
gif_dgm_file=${PathScratch}/GIF_DGM_${ID}.nii.gz
gif_dgm_bis_file=${PathScratch}/GIF_DGM_${ID}_bis.nii.gz
gif_wm_file=${PathScratch}/GIF_WM_${ID}.nii.gz
gif_wmi_file=${PathScratch}/GIF_WMI_${ID}.nii.gz
gif_prior_scratch=${PathScratch}/GIF_prior_${ID}.nii.gz
gif_tiv_scratch=${PathScratch}/GIF_${ID}_TIV.nii.gz
mask_file=${PathScratch}/Mask_${ID}.nii.gz
gif_b1_file=${PathScratch}/GIF_${ID}_B1.nii.gz
flair_file=${PathScratch}/FLAIR_${ID}.nii.gz
t1_file=${PathScratch}/T1_${ID}.nii.gz
aff_flair_to_t1=${PathScratch}/Aff_FLAIRtoT1.txt
aff_t1_to_flair=${PathScratch}/Aff_T1toFLAIR.txt
aff_transform_file=${PathScratch}/${ID}_AffTransf.txt
cpp_file=${PathScratch}/${ID}_cpp.nii.gz
gif_cgm_file=${PathScratch}/GIF_CGM_${ID}.nii.gz
gif_wmi_bis_file=${PathScratch}/GIF_WMI_${ID}_bis.nii.gz
gif_brainstem_file=${PathScratch}/GIF_Brainstem_${ID}.nii.gz
icbm_dgm_file=${PathScratch}/ICBM_DGM_${ID}.nii.gz
icbm_cgm_file=${PathScratch}/ICBM_CGM_${ID}.nii.gz
icbm_icsf_file=${PathScratch}/ICBM_ICSF_${ID}.nii.gz
icbm_ecsf_file=${PathScratch}/ICBM_ECSF_${ID}.nii.gz
icbm_gm_file=${PathScratch}/ICBM_GM_${ID}.nii.gz
icbm_csfs_file=${PathScratch}/ICBM_CSFs_${ID}.nii.gz
artefacts_file=${PathScratch}/${ID}_Artefacts.nii.gz
artefacts_map_file=${PathScratch}/Artefacts_${ID}.nii.gz
biasm_txt_file=${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.txt
biasm_nii_file=${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.nii.gz
insula_file=${PathScratch}/Insula_${ID}.nii.gz
correction_file=${PathScratch}/Correction_${ID}.nii.gz
lesion_corr_file=${PathScratch}/Lesion_${ID}_corr.nii.gz
# Lesion and SPCP files
primary_lesions_file=${PathScratch}/PrimaryLesions_${ID}.nii.gz
secondary_lesions_file=${PathScratch}/SecondaryLesions_${ID}.nii.gz
primary_spcp_file=${PathScratch}/PrimarySPCP_${ID}_${OptTxt}.nii.gz
secondary_spcp_file=${PathScratch}/SecondarySPCP_${ID}_${OptTxt}.nii.gz
merged_lesion_spcp_file=${PathScratch}/MergedLesion_${ID}_SPCP_${OptTxt}.nii.gz
merged_lesion_corr_file=${PathScratch}/MergedLesion_${ID}_${OptTxt}_corr.nii.gz
# Ventricle and SPCP intermediate files
ventricle_lining_file=${PathScratch}/VentricleLining.nii.gz
temp1_file=${PathScratch}/Temp1.nii.gz
potential_inter_ventr_file=${PathScratch}/PotentialInterVentr.nii.gz
expanded_ventr_lin_file=${PathScratch}/ExpandedVentrLin.nii.gz
potential_cp_file=${PathScratch}/PotentialCP.nii.gz
seg_dgm_file=${PathScratch}/SegDGM.nii.gz
potential_sp3_file=${PathScratch}/PotentialSP3.nii.gz
potential_spcp_file=${PathScratch}/PotentialSPCP.nii.gz
potential_spcp2_file=${PathScratch}/PotentialSPCP2.nii.gz
potential_spcp3_file=${PathScratch}/PotentialSPCP3.nii.gz
wmdil_file=${PathScratch}/WMDil.nii.gz
wmdil2_file=${PathScratch}/WMDil2.nii.gz
ventr_lin_file=${PathScratch}/VentrLin.nii.gz
ventr_lin_sp_file=${PathScratch}/VentrLinSP.nii.gz
ventr1_file=${PathScratch}/Ventr1.nii.gz
# DGM and WM correction files
dgm_file=${PathScratch}/${ID}_DGM.nii.gz
wmmin_in_file=${PathScratch}/WMminIn_${ID}.nii.gz
infra_dgm_file=${PathScratch}/InfraDGM_${ID}.nii.gz
# ICBM registration intermediate files
icbm_aladin_file=${PathScratch}/${ID}_ICBM_aladin.nii.gz
icbm_f3d_file=${PathScratch}/${ID}_ICBM_f3d.nii.gz
# BiASM output and GIF prior files
biasm_g_file=${PathScratch}/${ModalitiesTot}_BiASMG_${ID}.nii.gz
gif_csf_file=${PathScratch}/GIF_CSF_${ID}.nii.gz
gif_out_file=${PathScratch}/GIF_Out_${ID}.nii.gz

# GIF files declaration
tiv_pattern=(${PathGIF}/*TIV*)
tiv_file=${tiv_pattern[0]}
gif_prior_pattern=(${PathGIF}/*prior*)
gif_prior_file=${gif_prior_pattern[0]}



mkdir -p "$PathResults"

NameScript=${PathResults}/ScriptBaMoS_${PN}_${Opt}.sh
echo \#\!/bin/bash >${NameScript}
echo "mkdir -p ${PathScratch} " >>${NameScript}
ChangePathString="-inChangePath ${PathScratch}/"

echo "cp ${PathResults}/* ${PathScratch}/." >>${NameScript}

StringPriorsIO=" "

echo "echo \"Reorientation of T1 TIV mask and binarisation\"" >>${NameScript}

echo "cp ${tiv_file} ${mask_file}" >>${NameScript}
#echo "${PathFSL}/fslreorient2std ${PathGIF}/*TIV* ${PathScratch}/Mask_${ID}.nii.gz" >> ${NameScript}
echo "${seg_maths_bin} ${mask_file} -bin -odt char ${gif_b1_file}" >>${NameScript}


echo "echo \"Reorientation of GIF priors and creation of GIF priors array\"" >>${NameScript}
echo "cp ${gif_prior_file} ${gif_prior_scratch}" >>${NameScript}
#echo "${PathFSL}/fslreorient2std ${gif_prior_file} ${gif_prior_scratch}" >> ${NameScript}

# FIX #30: registration now lives exclusively inside the if/else branches below.
# The unconditional T1 copy + FLAIR->T1 reg block that previously ran before
# this if/else (causing a duplicate registration when Space=1) has been removed.
if ((${Space} > 1)); then
	echo "cp ${ImageFLAIR} ${flair_file}" >>${NameScript}
	echo "${seg_maths_bin} ${flair_file} -odt float ${flair_file}" >>${NameScript}

	echo "echo \"Alignment of T1 image to FLAIR\"" >>${NameScript}
	echo "${reg_aladin_bin} -ref ${flair_file} -flo ${ImageT1} -aff ${aff_t1_to_flair} -res ${t1_file} -rigOnly" >>${NameScript}
	echo "${seg_maths_bin} ${t1_file} -odt float ${t1_file}" >>${NameScript}
	for p in Segmentation prior TIV Parcellation; do
        temp_p_file_pattern=(${PathGIF}/*${p}*)
        temp_p_file=${temp_p_file_pattern[0]}
        temp_gif_file=${PathScratch}/GIF_${ID}_${p}.nii.gz

		echo "${reg_resample_bin} -inter 1 -flo ${temp_p_file} -res ${temp_gif_file} -ref ${flair_file} -aff ${aff_t1_to_flair}" >>${NameScript}
		echo "${reg_resample_bin} -inter 1 -flo ${temp_p_file} -res ${temp_gif_file} -ref ${flair_file} -aff ${aff_t1_to_flair}" >>${NameScript}
	done

else
	echo "echo \"Reorientation to standard space of T1 image\"" >>${NameScript}
	echo "cp ${ImageT1} ${t1_file}" >>${NameScript}
	echo "${seg_maths_bin} ${t1_file} -odt float ${t1_file}" >>${NameScript}
	echo "echo \"Alignment of FLAIR image to T1\"" >>${NameScript}
	echo "${reg_aladin_bin} -ref ${t1_file} -flo ${ImageFLAIR} -aff ${aff_flair_to_t1} -res ${flair_file} -rigOnly" >>${NameScript}
	echo "${seg_maths_bin} ${flair_file} -odt float ${flair_file}" >>${NameScript}
	for p in Segmentation prior TIV Parcellation; do
        gif_file_id_pattern=${PathScratch}/GIF_${ID}_${p}.nii.gz
		gif_file_pattern_id=${PathScratch}/GIF_${p}_${ID}.nii.gz
        temp_file_pattern=(${PathGIF}/*${p}*)
        temp_file_p=${temp_file_pattern[0]}
		echo "cp ${temp_file_p} ${gif_file_id_pattern}" >>${NameScript}
		echo "cp ${temp_file_p} ${gif_file_pattern_id}" >>${NameScript}
	done
fi

echo "${seg_maths_bin} ${gif_tiv_scratch} -odt char ${gif_b1_file}" >>${NameScript}
echo "${seg_maths_bin} ${gif_tiv_scratch} -odt char ${gif_b1_file}" >>${NameScript}

#  If GIF not performed in FLAIR space, need to resample parcellation to FLAIR space
#    echo "echo \"Reorientation of GIF parcellation\"" >> ${NameScript}
#    echo "cp ${PathGIF}/*Parcellation* ${gif_parcellation_file} " >> ${NameScript}
#echo "${PathFSL}/fslreorient2std ${PathGIF}/*Parcellation* ${gif_parcellation_file} " >> ${NameScript}

echo "echo \"Reorientation of GIF segmentation\"" >>${NameScript}
#echo "cp ${PathGIF}/*Segmentation* ${PathScratch}/GIF_Segmentation_${ID}.nii.gz " >> ${NameScript}
#  echo "${PathFSL}/fslreorient2std ${PathGIF}/*Segmentation* ${PathScratch}/GIF_Segmentation_${ID}.nii.gz " >> ${NameScript}

echo "find ${PathScratch} -maxdepth 1 -type f ! -name 'GIF_*' -exec cp {} ${PathResults}/. \;" >>${NameScript}

PriorsArray=("Out" "CSF" "CGM" "WMI" "DGM" "Brainstem")
array_Priors=()

for ((p = 0; p < 6; p++)); do
    gif_prior_scratch_temp=${PathScratch}/GIF_${PriorsArray[p]}_${ID}.nii.gz
	echo " ${seg_maths_bin} ${PathScratch}/GIF_prior* -tp ${p} ${gif_prior_scratch_temp}" >>${NameScript}
done

echo "echo \"Creation of artefact map based on T1 parcellation and registration to FLAIR\"" >>${NameScript}

if ((${#ArtefactArray[@]} > 0)); then
	stringAddition="${PathScratch}/${ID}_ArtConstruction_0.nii.gz "
	for ((i = 0; i < ${#ArtefactArray[@]}; i++)); do
		Value=${ArtefactArray[i]}
		ValueMin=$(echo "$Value - 0.5" | bc -l)
		ValueMax=$(echo "$Value + 0.5" | bc -l)
		echo "${seg_maths_bin} ${gif_parcellation_file} -thr $ValueMin -uthr $ValueMax -bin ${PathScratch}/${ID}_ArtConstruction_${i}.nii.gz" >>${NameScript}
		stringAddition="${stringAddition} -add ${PathScratch}/${ID}_ArtConstruction_${i}.nii.gz "
	done
	echo "${seg_maths_bin} ${stringAddition} -bin ${artefacts_file}" >>${NameScript}
fi

echo "rm ${PathScratch}/*_ArtConstruction_*" >>${NameScript}

if [ "$JUMP_START" -eq 0 ]; then
	echo "echo \"Registration of ICBM template and creation of ICBM atlases\"" >>${NameScript}
	# Pass explicit -res paths so niftyreg writes its resampled outputs into
	# scratch. Without -res, both binaries default to ./outputResult.nii(.gz)
	# which pollutes the Snakemake workdir.
	echo "${reg_aladin_bin} -ref ${t1_file} -flo ${PathICBM}/ICBM_Template.nii.gz -aff ${aff_transform_file} -res ${icbm_aladin_file}" >>${NameScript}
	echo "${reg_f3d_bin} -ref ${t1_file} -flo ${PathICBM}/ICBM_Template.nii.gz -aff ${aff_transform_file} -cpp ${cpp_file} -res ${icbm_f3d_file}" >>${NameScript}

    cpp_file_pattern=(${PathScratch}/*cpp*)
	echo "cp ${cpp_file_pattern[0]} ${PathResults}/." >>${NameScript}

	for p in CGM DGM ECSF ICSF Out WM; do
		echo "${reg_resample_bin} -ref ${t1_file} -flo ${PathICBM}/ICBM_${p}.nii.gz -cpp ${cpp_file} -res ${PathScratch}/ICBM_${p}_${ID}.nii.gz" >>${NameScript}
	done

	echo "${seg_maths_bin} ${icbm_dgm_file} -add ${icbm_cgm_file} ${icbm_gm_file}" >>${NameScript}

	echo "${seg_maths_bin} ${icbm_icsf_file} -add ${icbm_ecsf_file} ${icbm_csfs_file}" >>${NameScript}

	echo "${seg_maths_bin} ${icbm_dgm_file} -bin -mul ${gif_dgm_file} ${gif_dgm_bis_file}" >>${NameScript}

	echo "${seg_maths_bin} ${gif_cgm_file} -add ${gif_dgm_bis_file} ${gif_gm_file}" >>${NameScript}
	echo "${seg_maths_bin} ${gif_dgm_file} -sub ${gif_dgm_bis_file} -add ${gif_wmi_file} ${gif_wmi_bis_file}" >>${NameScript}

	echo "${seg_maths_bin} ${gif_wmi_bis_file} -add ${gif_brainstem_file} ${gif_wm_file}" >>${NameScript}

	arrayImage=()
	arrayModNumber=()
	echo "echo \"Modalities to put are ${arrayMod[0]}\"" >>${NameScript}
	for ((m = 0; m < ${#arrayMod[@]}; m++)); do
		for ((pos = 0; pos < ${#ArrayModalities[@]}; pos++)); do
			TmpModa=${ArrayModalities[pos]}
			TmptestModa=${arrayMod[m]}
			Subtracted="${TmptestModa/$TmpModa/}"
			#echo ${#Subtracted}
			if ((${#Subtracted} < ${#TmptestModa})); then
				# echo "Testing between $TmptestModa $TmpModa"
				FinModa=$((pos + 1))
				arrayModNumber=(${arrayModNumber[*]} $FinModa)
			fi
		done
		# FIX #35: use [@] instead of [*] to preserve separate array elements
		arrayImage=("${arrayImage[@]}" "${PathScratch}/${arrayMod[m]}_${ID}.nii.gz")
	done

	array_Priors=("${gif_gm_file}" "${gif_wm_file}" "${gif_csf_file}" "${gif_out_file}")
	echo "echo \"Segmentation SegBiASM\"" >>${NameScript}
	#${PathSeg}/Seg_BiASM -in 2 ${arrayImage[0]} ${arrayImage[1]} -priors 4 ${array_Priors[*]} -mask ${gif_b1_file} -out 2 ${PathScratch}/${arrayMod[0]}${arrayMod[1]}_BiASM_${ID}_${Opt}.nii.gz ${PathScratch}/${arrayMod[0]}${arrayMod[1]}_BiASMG_${ID}.nii.gz -txt_out ${PathScratch}/${arrayMod[0]}${arrayMod[1]}_BiASM_${ID}_${Opt}.txt -bc_order 3 -CovPriors 8 -BFP 1 -maxRunEM ${MRE} -AtlasSmoothing 1 1 -AtlasWeight 1 ${AW}  -SMOrder 0 -KernelSize 3 -PriorsKept 5 -VLkappa 1.5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW ${OW} -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr ${JC} -progMod 0 -priorDGM ${icbm_dgm_file} -TypicalityAtlas ${OptTA} ${StringPriorsIO}

	echo "${seg_biasm_bin} -VLkappa 3 -in 2 ${arrayImage[*]} -priors 4 ${array_Priors[*]} -mask ${gif_b1_file} -out 2 ${biasm_nii_file} ${biasm_g_file} -txt_out ${biasm_txt_file} -bc_order 3 -CovPriors 8 -BFP 1 -maxRunEM ${MRE} -AtlasSmoothing 1 1 -AtlasWeight 1 ${AW}  -SMOrder 0 -KernelSize 3 -PriorsKept 5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW ${OW} -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr ${JC} -progMod 0 -priorDGM ${icbm_dgm_file} -TypicalityAtlas ${OptTA} -scratchDir ${PathScratch}" >>${NameScript}

	# Preserve MRFOut_*: it is auto-generated by Seg_BiASM when -GMRF is set
	# (Seg_BiASM.cpp:1170) and is read by the subsequent Seg_Analysis call via
	# the .txt tree file's recorded MRF path. The original glob MRF* was too
	# aggressive and deleted MRFOut_*, causing Seg_Analysis to fail with
	# "nifti_image_read: failed to find header file for .../MRFOut_*.nii.gz".
	echo "find ${PathScratch} -maxdepth 1 -name 'BG*' -delete" >>${NameScript}
	echo "find ${PathScratch} -maxdepth 1 -name 'MRF*' -not -name 'MRFOut_*' -delete" >>${NameScript}
    t1_flair_pattern=(${PathScratch}/T1FLAIR*)
    data_t1_flair_pattern=(${PathScratch}/Data*T1FLAIR*)
	echo "cp ${t1_flair_pattern[0]} ${PathResults}/. " >>${NameScript}
	echo "cp ${data_t1_flair_pattern[0]} ${PathResults}/. " >>${NameScript}

	echo "echo \"Lesion segmentation\"" >>${NameScript}
	echo "${seg_analysis_bin}  -inTxt2 ${biasm_txt_file} ${biasm_nii_file} -mask ${gif_b1_file} -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa 2 1 3 -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${icbm_icsf_file} -inPriorsDGM ${icbm_dgm_file} -inPriorsCGM ${icbm_cgm_file} -inPriorsECSF ${icbm_ecsf_file} -TO 1 -ParcellationIn ${gif_parcellation_file} -typeVentrSeg ${TVS} -Simple -Secondary 60 -juxtaCorr ${JC}  -SP ${OptSP} -LevelCorrection ${OptCL} ${ChangePathString} -LesWMI ${OptWMI} -Neigh 18 -scratchDir ${PathScratch}" >>${NameScript}

fi

if [ "$JUMP_START" -eq 1 ]; then

	for p in CGM DGM ECSF ICSF Out WM; do
		echo "${reg_resample_bin} -ref ${t1_file} -flo ${PathICBM}/ICBM_${p}.nii.gz -cpp ${cpp_file} -res ${PathScratch}/ICBM_${p}_${ID}.nii.gz" >>${NameScript}
	done

	PriorsArray=("Out" "CSF" "CGM" "WMI" "DGM" "Brainstem")
	array_Priors=()

	for ((p = 0; p < 6; p++)); do
		echo " ${seg_maths_bin} ${PathScratch}/GIF_prior* -tp ${p} ${PathScratch}/GIF_${PriorsArray[p]}_${ID}.nii.gz" >>${NameScript}
	done

	echo "${seg_maths_bin} ${icbm_dgm_file} -add ${icbm_cgm_file} ${icbm_gm_file}" >>${NameScript}

	echo "${seg_maths_bin} ${icbm_icsf_file} -add ${icbm_ecsf_file} ${icbm_csfs_file}" >>${NameScript}

	echo "${seg_maths_bin} ${icbm_dgm_file} -bin -mul ${gif_dgm_file} ${gif_dgm_bis_file}" >>${NameScript}

	echo "${seg_maths_bin} ${gif_cgm_file} -add ${gif_dgm_bis_file} ${gif_gm_file}" >>${NameScript}
	echo "${seg_maths_bin} ${gif_dgm_file} -sub ${gif_dgm_bis_file} -add ${gif_wmi_file} ${gif_wmi_bis_file}" >>${NameScript}

	echo "${seg_maths_bin} ${gif_wmi_bis_file} -add ${gif_brainstem_file} ${gif_wm_file}" >>${NameScript}
fi

# FIX #34: The stray unconditional Seg_Analysis call that previously appeared here
# (between the two fi statements) has been removed. It ran regardless of JUMP_START,
# causing a missing MRFOut file error when JUMP_START=1 (Seg_BiASM not yet run).
# The correct call is inside JUMP_START=0 above; the definitive final call is below.

echo "cp ${PathScratch}/LesionCorrected*WS3WT3WC1* ${primary_lesions_file}" >>${NameScript}
echo "cp ${PathScratch}/SecondaryCorrected*WS3WT3WC1* ${secondary_lesions_file}" >>${NameScript}

echo "cp ${PathScratch}/Primary* ${PathResults}/. " >>${NameScript}
echo "cp ${PathScratch}/Secondary* ${PathResults}/. " >>${NameScript}

# FIX #35: use rm -f so missing intermediate files don't cause fatal errors
echo "rm -f ${PathScratch}/DataR* ${PathScratch}/LesionT* ${PathScratch}/Summ*" >>${NameScript}

echo "echo \"Correction for septum pellucidum\"" >>${NameScript}
echo "${seg_maths_bin} ${gif_parcellation_file} -thr 65.5 -uthr 67.5 -bin ${ventricle_lining_file} " >>${NameScript}
echo "${seg_maths_bin} ${gif_parcellation_file}  -equal 52 -bin -euc -uthr 0 -abs ${temp1_file}" >>${NameScript}
echo "${seg_maths_bin} ${gif_parcellation_file}  -equal 53 -bin -euc -uthr 0 -abs -add ${temp1_file} -uthr 5 -bin ${potential_inter_ventr_file}" >>${NameScript}

echo "${seg_maths_bin} ${ventricle_lining_file} -dil 1 ${expanded_ventr_lin_file}" >>${NameScript}
echo "${seg_maths_bin} ${gif_parcellation_file}  -thr 49.5 -uthr 53.5 -bin -sub ${expanded_ventr_lin_file} -thr 0 -bin ${potential_cp_file}" >>${NameScript}

#${seg_maths_bin} ${PathGIF}/*Parcellation*  -thr 122.5 -uthr 124.5 -bin ${PathScratch}/Temp4.nii.gz

echo "${seg_maths_bin} ${gif_segmentation_file} -tp 4 -thr 0.5 -bin -sub ${PathScratch}/VentricleLining* -thr 0 ${seg_dgm_file}" >>${NameScript}
echo "${seg_maths_bin} ${gif_parcellation_file}  -equal 87 -mul -1 -add ${potential_inter_ventr_file} -sub ${seg_dgm_file} -thr 0 ${potential_sp3_file}" >>${NameScript}
echo "${seg_maths_bin}  ${potential_sp3_file} -add ${potential_cp_file} ${potential_spcp_file}" >>${NameScript}
echo "${seg_maths_bin} ${PathScratch}/ICBM_ICSF* -thr 0.3 -bin -mul ${potential_spcp_file} -add ${potential_sp3_file} -thr 0.2  -bin  ${potential_spcp2_file}" >>${NameScript}

echo "${seg_maths_bin} ${gif_parcellation_file} -thr 83.5 -uthr 84.5 -bin -dil 5 ${wmdil_file}" >>${NameScript}
echo "${seg_maths_bin} ${gif_parcellation_file} -thr 91.5 -uthr 92.5 -bin -dil 5 ${wmdil2_file}" >>${NameScript}
echo "${seg_maths_bin} ${potential_spcp2_file} -sub ${wmdil_file} -sub ${wmdil2_file} -thr 0 ${potential_spcp3_file} " >>${NameScript}

echo "${seg_maths_bin} ${PathScratch}/PrimaryLesions_* -sub ${potential_spcp3_file} -thr 0 ${primary_spcp_file}" >>${NameScript}
echo "${seg_maths_bin} ${PathScratch}/SecondaryLesions_* -sub ${potential_spcp3_file} -thr 0 ${secondary_spcp_file}" >>${NameScript}
echo "${seg_maths_bin} ${PathScratch}/PrimarySPCP*${OptTxt}* -merge 1 4 ${PathScratch}/SecondarySPCP_*${OptTxt}* -tmax ${merged_lesion_spcp_file}" >>${NameScript}
# FIX #31: was reading from PathResults (file not yet copied there); corrected to PathScratch
echo "${seg_maths_bin} ${merged_lesion_spcp_file} -sub ${artefacts_file} -thr 0 ${merged_lesion_spcp_file}" >>${NameScript}

echo "echo \"Refined correction for third ventricle\"" >>${NameScript}
echo "${seg_maths_bin} ${gif_parcellation_file}  -sub ${gif_parcellation_file} ${artefacts_map_file} " >>${NameScript}

for ((i = 0; i < ${#ListCorr[@]}; i++)); do
	echo "${seg_maths_bin} ${gif_parcellation_file}  -equal ${ListCorr[i]} -bin -add ${artefacts_map_file} " >>${NameScript}
done
echo "${seg_maths_bin} ${gif_parcellation_file}  -thr 65.5 -uthr 67.5 ${ventr_lin_file} " >>${NameScript}
echo "${seg_maths_bin} ${gif_parcellation_file}  -equal 87 -add ${ventr_lin_file} ${ventr_lin_file} " >>${NameScript}
echo "${seg_maths_bin} ${gif_parcellation_file}  -equal 47 -euc  -abs -uthr 5 -bin -mul ${ventr_lin_file} ${ventr_lin_sp_file} " >>${NameScript}

# FIX #33: removed double-slash typo (was ${PathScratch}//Ventr1.nii.gz)
echo "${seg_maths_bin} ${gif_parcellation_file}  -equal 52 -bin -euc  -abs  ${ventr1_file} " >>${NameScript}

echo "${seg_maths_bin} ${gif_parcellation_file}  -equal 53 -bin -euc -abs -uthr 5 -sub  ${ventr1_file} -thr -5 -uthr 5 -bin -mul  ${gif_parcellation_file}  -equal 52 -bin -euc -abs -uthr 5 -mul ${ventr_lin_file}  -add ${ventr_lin_sp_file} ${ventr_lin_sp_file}" >>${NameScript}
echo "${seg_maths_bin} ${ventr_lin_sp_file} -add ${artefacts_map_file} -mul -1 -add ${merged_lesion_spcp_file} -thr 0 ${merged_lesion_corr_file} " >>${NameScript}

echo "${seg_analysis_bin} -LesWMI ${OptWMI} ${OptInfarcts} -inLesCorr ${merged_lesion_corr_file}  -inTxt2 ${biasm_txt_file} ${biasm_nii_file} -mask ${gif_b1_file} -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa ${#arrayModNumber[@]} ${arrayModNumber[*]} -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${icbm_icsf_file} -inPriorsDGM ${icbm_dgm_file} -inPriorsCGM ${icbm_cgm_file} -inPriorsECSF ${icbm_ecsf_file} -TO 1 -juxtaCorr 1 -SP ${OptSP} -LevelCorrection ${OptCL} -inArtefact ${artefacts_file} ${ChangePathString} -ParcellationIn ${gif_parcellation_file} -typeVentrSeg 1 -outWM 1 -outConnect 1 -Neigh 6 -scratchDir ${PathScratch}" >>${NameScript}

echo "rm -f ${PathScratch}/LesionWeigh* ${PathScratch}/Binary* ${PathScratch}/WMDil* ${PathScratch}/WMCard* ${PathScratch}/LesionInit* ${PathScratch}/DataR* ${PathScratch}/DataT* ${PathScratch}/Summ* ${PathScratch}/LesSegHard* ${PathScratch}/Check* ${PathScratch}/BinaryNIV*" >>${NameScript}

# ${seg_maths_bin} ${gif_parcellation_file}  -equal 52 -bin -euc  -abs  ${ventr1_file}

# ${seg_maths_bin} ${gif_parcellation_file}  -equal 53 -bin -euc -abs -uthr 5 -sub  ${ventr1_file} -thr -5 -uthr 5 -bin -mul  ${gif_parcellation_file}  -equal 52 -bin -euc -abs -uthr 5 -mul ${ventr_lin_file}  -add ${ventr_lin_sp_file} ${ventr_lin_sp_file}
# ${seg_maths_bin} ${ventr_lin_sp_file} -add ${artefacts_map_file} -mul -1 -add ${PathScratch}/MergedLesion_${ID}_SPCP.nii.gz -thr 0 ${PathScratch}/MergedLesion_${ID}_corr.nii.gz

# echo "${seg_maths_bin} ${ventr_lin_sp_file} -add ${artefacts_map_file} -mul -1 -add ${PathScratch}/SecondarySPCP_${ID}.nii.gz -thr 0 ${PathScratch}/SecondarySPCP_${ID}_corr.nii.gz " >> ${NameScript}

# echo "${seg_maths_bin} ${ventr_lin_sp_file} -add ${artefacts_map_file} -mul -1 -add ${PathScratch}/PrimarySPCP_${ID}.nii.gz -thr 0 ${PathScratch}/PrimarySPCP_${ID}_corr.nii.gz " >> ${NameScript}

# echo "${seg_maths_bin} ${PathScratch}/PrimarySPCP_${ID}_corr.nii.gz -merge 1 4 ${PathScratch}/SecondarySPCP_${ID}_corr.nii.gz -tmax ${PathScratch}/MergedLesionSPCP_${ID}_corr.nii.gz" >> ${NameScript}

# echo "echo \"Optional second level of correction\"" >> ${NameScript}
# if ((flag_furtherCorr==1))
# then
#     echo "${seg_analysis_bin} -inArtefact ${artefacts_map_file} -inTxt2 ${biasm_txt_file} ${biasm_nii_file} -correct -connect -inLesCorr ${PathScratch}/MergedLesionSPCP_${ID}_corr.nii.gz  -mask ${gif_b1_file} -inPriorsICSF ${icbm_icsf_file} -inPriorsDGM ${icbm_dgm_file} -inPriorsCGM ${icbm_cgm_file} -inPriorsECSF ${icbm_ecsf_file} -outConnect 1 -outWM 1 -ParcellationIn ${gif_parcellation_file}  -inChangePath ${PathScratch}/ -inModa 2 1 3 -WeightedSeg 3 3 1 -LesWMI ${OptWMI} -Neigh 6 -scratchDir ${PathScratch}" >> ${NameScript}

# fi

if ((flag_corrDGM == 1)); then
	Array=(24 31 32 56 57 58 59 60 61 76 77 37 38)

	# Get segmentation of DGM from GIF

	if [ ! -f ${dgm_file} ]; then
		if ((${#Array[@]} > 0)); then
			stringAddition="${PathScratch}/${ID}_DGMConstruction_0.nii.gz "
			for ((k = 0; k < ${#Array[@]}; k++)); do
				Value=${Array[k]}
				echo "${seg_maths_bin} ${gif_parcellation_file} -equal ${Value} -bin ${PathScratch}/${ID}_DGMConstruction_${k}.nii.gz " >>${NameScript}
				stringAddition="${stringAddition} -add ${PathScratch}/${ID}_DGMConstruction_${k}.nii.gz "
			done
			echo "${seg_maths_bin} ${stringAddition} -bin -odt char ${dgm_file} " >>${NameScript}
			echo "rm -f ${PathScratch}/*DGMConstruction* " >>${NameScript}
		fi
	fi
	if [ ! -f ${PathScratch}/WMminIn* ]; then
		echo "${seg_maths_bin} ${gif_parcellation_file} -thr 81.5 -uthr 83.5 ${insula_file} " >>${NameScript}
		echo "${seg_maths_bin} ${gif_parcellation_file} -thr 89.5 -uthr 91.5 -add ${insula_file} -bin  ${insula_file} " >>${NameScript}
		echo "${seg_maths_bin} ${gif_parcellation_file} -thr 95.5 -uthr 97.5 -add ${insula_file} -bin  ${insula_file} " >>${NameScript}
		echo "${seg_maths_bin} ${gif_parcellation_file} -thr 79 -uthr 98 -bin -sub ${insula_file} ${wmmin_in_file}" >>${NameScript}
	fi
	echo "${seg_maths_bin} ${gif_parcellation_file} -thr 23 -uthr 45 -bin ${infra_dgm_file}" >>${NameScript}
	echo "${seg_maths_bin} ${PathScratch}/${ID}_DGM* -add ${PathScratch}/Infra* -add ${PathScratch}/Insula* -bin -dil 1 -sub ${wmmin_in_file} -thr 0 ${correction_file}" >>${NameScript}

	echo "${seg_maths_bin} ${PathScratch}/LesionMahal* -tp 2 -uthr 4 -bin -mul ${PathScratch}/Corr*Mer* -mul ${correction_file} -mul -1 -add ${PathScratch}/Corr*Mer* -thr 0 ${lesion_corr_file}" >>${NameScript}

	echo "${seg_analysis_bin} -LesWMI ${OptWMI}  -inLesCorr ${lesion_corr_file} -inTxt2 ${biasm_txt_file} ${biasm_nii_file} -mask ${gif_b1_file} -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa 2 1 3 -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${icbm_icsf_file} -inPriorsDGM ${icbm_dgm_file} -inPriorsCGM ${icbm_cgm_file} -inPriorsECSF ${icbm_ecsf_file} -TO 1 -juxtaCorr 1 -SP ${OptSP} -LevelCorrection ${OptCL} -inArtefact ${artefacts_file} ${ChangePathString} -ParcellationIn ${gif_parcellation_file} -typeVentrSeg 1 -outWM 1 -outConnect 1 -Neigh 6 -scratchDir ${PathScratch}" >>${NameScript}

fi

echo "rm -f ${PathScratch}/LesionWeigh* ${PathScratch}/Binary* ${PathScratch}/WMDil* ${PathScratch}/WMCard* ${PathScratch}/ICBM* ${PathScratch}/LesionInit* ${PathScratch}/DataR* ${PathScratch}/DataT* ${PathScratch}/Summ* ${PathScratch}/LesSegHard* ${PathScratch}/Check* ${PathScratch}/BinaryNIV*" >>${NameScript}
echo "cp ${PathScratch}/*Co* ${PathResults}/. 2>/dev/null || true" >>${NameScript}
echo "cp ${PathScratch}/*Co*.txt ${PathResults}/. 2>/dev/null || true" >>${NameScript}
echo "cp ${PathScratch}/LesionMahal* ${PathResults}/. 2>/dev/null || true" >>${NameScript}
echo "cp ${PathScratch}/Txt* ${PathResults}/. 2>/dev/null || true" >>${NameScript}
echo "cp ${PathScratch}/Out* ${PathResults}/. 2>/dev/null || true" >>${NameScript}
echo "cp ${PathScratch}/Autho* ${PathResults}/. 2>/dev/null || true" >>${NameScript}
echo "cp ${PathScratch}/*Infar* ${PathResults}/. 2>/dev/null || true" >>${NameScript}

echo "cp ${PathScratch}/*Artefacts.nii.gz ${PathResults}/." >>${NameScript}

# Run BaMoS
chmod +x ${NameScript}
${NameScript}