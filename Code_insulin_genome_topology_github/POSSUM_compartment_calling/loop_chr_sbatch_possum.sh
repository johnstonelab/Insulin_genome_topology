#!/bin/bash

#Assign command line arguments to variables
CHROM_SIZE_FILE=$1
HIC_FILE=$2
OUTPUT_DIRECTORY=$3
BIN_SIZE=$4

#Get the name of HiC file
base=$(basename "$HIC_FILE" .hic)

#Loop over each line in the chrom.size file, and assign the first column to CHR_NAME and the second column to CHR_LENGTH.

while read -r CHR_NAME CHR_LENGTH; do
 #construct the output file name
 FILE_OUT="${OUTPUT_DIRECTORY}/${base}_${BIN_SIZE}_${CHR_NAME}_output"

 echo "${CHR_NAME}"
 echo "${CHR_LENGTH}"
 echo "${base}"

 #Submit the job, the sbatch_submit_possum.sh should be in the same directory as this .sh
 sbatch sbatch_submit_possum.sh "$CHR_LENGTH" "$HIC_FILE" "$CHR_NAME" "$FILE_OUT" "$BIN_SIZE"

done < "$CHROM_SIZE_FILE"
