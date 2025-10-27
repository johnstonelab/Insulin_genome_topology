mm10_chrom_ordered_sizes <- read.csv("/n/data1/bidmc/medicine/kahn/lab/xyl/chr_size/mm10_ordered/mm10.chrom.ordered.sizes", sep = "\t", header = FALSE)

colnames(mm10_chrom_ordered_sizes) <- c("chr","length")

mm10_chrom_ordered_sizes_noYM <- mm10_chrom_ordered_sizes[1:20,]

save.image("/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_4_flam/Mic5_4_2_clear_genes/Mic5_4_2_Ctl_analysis/Ctl_chr.RData")

library(FLAMINGOrLite)
library(Matrix)
library(devtools)

flam_list_Ctl <- list()

for (i in 1:nrow(mm10_chrom_ordered_sizes_noYM)){
  chr_name_mm10 <- as.character(mm10_chrom_ordered_sizes_noYM$chr[i])
  chr_length <- mm10_chrom_ordered_sizes_noYM$length[i]
  
  print(chr_name_mm10)
  print(chr_length)
  
  #Generate flamingo
  flam_tmp = flamingo_main(
    hic_data='/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_4_flam/no_chr_FLAMINGO_hic_2nd/merge_nochr_Ctl_MicroC_intra.hic',
    file_format='hic',
    domain_res=1e6,
    frag_res=5000,
    chr_name=chr_name_mm10,
    normalization="KR",
    nThread=20)
  
  #Add the flam_tmp to the list
  flam_list_Ctl[[chr_name_mm10]] <- flam_tmp
  
  #Generate output path
  opt_path <- paste0("/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_4_flam/Mic5_4_2_clear_genes/Mic5_4_2_Ctl_analysis/Ctl_",chr_name_mm10,".vtk")
  print(opt_path)

  #write the output
  write.vtk(flam_tmp[,c("x","y","z")],flam_tmp$start,chr_name_mm10,opt_path)
}

save.image("/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_4_flam/Mic5_4_2_clear_genes/Mic5_4_2_Ctl_analysis/Ctl_flam_output.RData")

