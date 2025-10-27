#!/bin/bash

#SBATCH -c 2
#SBATCH -t 0-5:00
#SBATCH -p short
#SBATCH --mem=50G

module load gcc/6.2.0 R/4.1.1

cd /n/data1/bidmc/medicine/kahn/lab/xyl/Packages_XYL/POSSUMM/EigenVector/R

#run the script
Rscript /n/data1/bidmc/medicine/kahn/lab/xyl/Packages_XYL/POSSUMM/EigenVector/R/eigFromHicRscript.R \
-n KR \
-s "$1" \
"$2" \
"$3" \
"$4" \
"$5"
