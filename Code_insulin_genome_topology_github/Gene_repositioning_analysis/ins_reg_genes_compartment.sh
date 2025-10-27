#!/bin/bash

# Define the directory containing the BED files
BED_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed"
BIGWIG_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/Mic5_2_2_bigwig_adjusted"
OUTPUT_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/matrix/ins_reg_matrix"

# Load modules
module load gcc/6.2.0  python/2.7.12
module load deeptools/3.0.2

#Generate matrix files using a for loop
for BED_FILE in ${BED_DIR}/*.bed; do
  # Extract the basename for the BED file (e.g., "chr1" from "chr1.bed")
  FILE_NAME=$(basename "${BED_FILE}" .bed)
  
  echo "${FILE_NAME}"

  # Define the output file name
  OUTPUT_FILE="${OUTPUT_DIR}/matrix_${FILE_NAME}_scale.gz"

  # Run computeMatrix
  computeMatrix scale-regions \
    -R "${BED_FILE}" \
    -S ${BIGWIG_DIR}/*.bigwig \
    -p 6 \
    --regionBodyLength 20000 \
    -a 10000 -b 10000 \
    -o "${OUTPUT_FILE}"
done

##Plot the compartment

MATRIX_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/matrix/ins_reg_matrix"
OUTPUT_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/plots/ins_reg_plot"

#Generate png files using a for loop

for MATRIX in ${MATRIX_DIR}/*.gz; do
  # Extract the basename 
  NAME=$(basename "${MATRIX}" .gz)
  
  echo "${NAME}"

  # Define the output file name
  OUTPUT_FILE="${OUTPUT_DIR}/${NAME}.png"

  # plot the matrix
  plotProfile -m "${MATRIX}" \
	-out "${OUTPUT_FILE}" \
	--perGroup  --plotTitle "" \
	-z "" \
	--startLabel "TSS" \
	--endLabel "TES" \
	--colors darkblue orange
done


