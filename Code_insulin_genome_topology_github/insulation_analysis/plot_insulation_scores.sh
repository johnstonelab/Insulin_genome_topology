#!/bin/bash

# Define the directories
BED_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_4_insulation/ins_reg_and_all_genes_bed"
BIGWIG_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_4_insulation/Mic8_1_4_bw"
MATRIX_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_4_insulation/for_paper_TSS_TES_Ctl/matrix_ctl"

PLOT_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_4_insulation/for_paper_TSS_TES_Ctl/plot_ctl"

# Load modules
module load gcc/6.2.0 python/2.7.12
module load deeptools/3.0.2

#loop through the file and compute matrix
for BED_FILE in "${BED_DIR}"/all_annotation_liver.bed; do
    # Extract the basename for the BED file
    FILE_NAME=$(basename "${BED_FILE}" .bed)

    echo "${FILE_NAME}"

    # Define the output file name
    MATRIX_FILE="${MATRIX_DIR}/matrix_${FILE_NAME}_TSS_10k.gz"

    # Run computeMatrix
    computeMatrix reference-point \
        --referencePoint TSS \
        -R "${BED_FILE}" \
        -S "${BIGWIG_DIR}"/Ctl_score.bigwig \
        -p 6 \
        -a 10000 -b 10000 \
        -o "${MATRIX_FILE}"
done


#Plot an aggregated plot
for MATRIX in "${MATRIX_DIR}"/*_TSS_10k.gz; do
    # Extract the basename
    NAME=$(basename "${MATRIX}" .gz)

    echo "${NAME}"

    # Define the output file name
    PLOT_FILE="${PLOT_DIR}/${NAME}.png"

    # plot the data
    plotProfile -m "${MATRIX}" \
        -out "${PLOT_FILE}" \
        --perGroup  --plotTitle "" \
        -z "" \
        --refPointLabel "TSS" \
        --colors darkblue
done


#Plot a heatmap
for MATRIX in "${MATRIX_DIR}"/*_TSS_10k.gz; do
    # Extract the basename
    NAME=$(basename "${MATRIX}" .gz)

    echo "${NAME}"

    # Define the output file name
    PLOT_FILE="${PLOT_DIR}/${NAME}_heatmap_blues.png"

    # plot the heatmap

    plotHeatmap -m "${MATRIX}" \
        -out "${PLOT_FILE}" \
        --colorMap Blues \
        --missingDataColor 0.5 \
        --whatToShow 'heatmap and colorbar'
done