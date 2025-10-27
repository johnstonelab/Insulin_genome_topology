#!/bin/bash

#For upregulated genes
# Define the directory containing the BED files
MATRIX_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic5_6_bw_subtract/matrix_Mic5_6/matrix_Mic5_6_2_2_up_only"
PLOT_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic5_6_bw_subtract/plots_same_y_Mic5622"

# Load modules
module load gcc/6.2.0 python/2.7.12
module load deeptools/3.0.2

# Generate plots using a for loop
for MATRIX in ${MATRIX_DIR}/*.gz; do
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
	--startLabel "TSS" \
	--endLabel "TES" \
	--colors red \
    --yMin -0.0004 \
    --yMax 0.00095 \
    --legendLocation none
done


#For downregulated genes
# Define the directory containing the BED files
MATRIX_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic5_6_bw_subtract/matrix_Mic5_6/matrix_Mic5_6_2_2_down_only"
PLOT_DIR="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic5_6_bw_subtract/plots_same_y_Mic5622"

# Load modules
module load gcc/6.2.0 python/2.7.12
module load deeptools/3.0.2

#Generate plots
for MATRIX in ${MATRIX_DIR}/*.gz; do
  # Extract the basename
  NAME=$(basename "${MATRIX}" .gz)

  echo "${NAME}"

  # Define the output file name
  PLOT_FILE="${PLOT_DIR}/${NAME}.png"

  # Plot
  plotProfile -m "${MATRIX}" \
	-out "${PLOT_FILE}" \
	--perGroup  --plotTitle "" \
	-z "" \
	--startLabel "TSS" \
	--endLabel "TES" \
	--colors black \
    --yMin -0.0005 \
    --yMax 0.00025 \
    --legendLocation none
done
