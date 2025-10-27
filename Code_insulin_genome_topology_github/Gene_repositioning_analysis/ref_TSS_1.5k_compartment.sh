#!/bin/bash

# Define the directory containing the BED files
BED_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed"
BIGWIG_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/Mic5_2_2_bigwig_adjusted"
MATRIX_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic5_2_3_2_2_reg_ref_TSS/matrix_1.5k"

PLOT_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic5_2_3_2_2_reg_ref_TSS/plot_1.5k"

# Load modules
module load gcc/6.2.0 python/2.7.12
module load deeptools/3.0.2

#Use a for loop to compute matrix
for BED_FILE in "${BED_DIR}"/*.bed; do
    # Extract the basename for the BED file
    FILE_NAME=$(basename "${BED_FILE}" .bed)

    echo "${FILE_NAME}"

    # Define the output file name
    MATRIX_FILE="${MATRIX_DIR}/matrix_${FILE_NAME}_TSS_10k.gz"

    # Run computeMatrix
    computeMatrix reference-point \
        --referencePoint TSS \
        -R "${BED_FILE}" \
        -S "${BIGWIG_DIR}"/*.bigwig \
        -p 6 \
        -a 1500 -b 1500 \
        -o "${MATRIX_FILE}"
done

#plot data

for MATRIX in "${MATRIX_DIR}"/*.gz; do
    # Extract the basename
    NAME=$(basename "${MATRIX}" .gz)

    echo "${NAME}"

    # Define the output file name
    PLOT_FILE="${PLOT_DIR}/${NAME}.png"

    # plot
    plotProfile -m "${MATRIX}" \
        -out "${PLOT_FILE}" \
        --perGroup  --plotTitle "" \
        -z "" \
        --refPointLabel "TSS" \
        --colors darkblue orange
done
