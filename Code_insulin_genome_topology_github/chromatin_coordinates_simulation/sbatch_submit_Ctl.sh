#!/bin/bash

#SBATCH -c 20
#SBATCH -t 4-23:00
#SBATCH -p highmem
#SBATCH --mem=650G

module load gcc/9.2.0 R/4.1.2

#enter the directory

cd /n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_4_flam/Mic5_4_2_clear_genes/Mic5_4_2_Ctl_analysis

#run the script
Rscript /n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_4_flam/Mic5_4_2_clear_genes/Mic5_4_2_Ctl_analysis/Ctl_simulation_R.R
