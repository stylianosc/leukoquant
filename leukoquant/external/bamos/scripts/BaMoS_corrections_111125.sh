# Modified by Stylianos Charalampous on 11/11/2025

while getopts "o:r:g:j:i:" OPT; do
    case $OPT in
        o) output_folder="${OPTARG}"
           ;;
        r) bamos_results_folder="${OPTARG}"
           ;;
        g) gif_results_folder="${OPTARG}"
           ;;
        i) ID="${OPTARG}"
           ;;
    esac
done

script_path=${output_folder}/script_BaMoS_corrections.sh

if [ ! -d "$output_folder" ]; then
    mkdir -p "$output_folder"
fi

exec > $script_path 2>&1

echo \#\!/bin/sh

#LesionFile in the format: Correct*1Les*
echo 'bamos_lesion_file_path=$(find "'$bamos_results_folder'" -maxdepth 1 -path "*Correct*1Les*" )'

#Connectivity file in the format: Connect*1Les*
echo 'bamos_connectivity_file_path=$(find "'$bamos_results_folder'" -maxdepth 1 -path "*Connect*1Les*" )'

#Label file in the format: Txt*1Les*
echo 'bamos_label_file_path=$(find "'$bamos_results_folder'" -maxdepth 1 -path "*Txt*1Les*")'

echo 'gif_parcellation_file_path=$(find "'$gif_results_folder'" -maxdepth 1 -path "*NeuroMorph_Parcellation*" )'

echo 'python_cmd="/usr/bin/env python"'
echo 'bamos_corr_cmd="$python_cmd /leukoquant/leukoquant/external/bamos/scripts/correction_lesions_111125.py -les ${bamos_lesion_file_path} \
                                                                                                                        -parc ${gif_parcellation_file_path} \
                                                                                                                        -connect ${bamos_connectivity_file_path} \
                                                                                                                        -dil 1 \
                                                                                                                        -label ${bamos_label_file_path} \
                                                                                                                        -corr choroid cortex sheet \
                                                                                                                        -id ${ID}"'
                                                                                           
                                                                                            
echo '$bamos_corr_cmd'

exec > ${output_folder}/log.txt 2>&1
