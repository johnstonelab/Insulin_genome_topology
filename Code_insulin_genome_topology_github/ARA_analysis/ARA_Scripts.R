#Aggregated Region Analysis (ARA)

library(GENOVA)

#Input .hic files
merge_Ctl_MicroC <- load_contacts(signal_path = '/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_2_POSSUMM_Nova/merged_hic_files/with_chr_merged_hic_POSSUMM/merge_Ctl_MicroC.hic',
                                     sample_name = "Ctl",
                                     resolution = 1250,
                                     balancing = 'KR', # this is the default
                                     colour = "blue")

merge_Ins_MicroC <- load_contacts(signal_path = '/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_2_POSSUMM_Nova/merged_hic_files/with_chr_merged_hic_POSSUMM/merge_Ins_MicroC.hic',
                                  sample_name = "Ins",
                                  resolution = 1250,
                                  balancing = 'KR', # this is the default
                                  colour = "orange")


#ARA analysis for all PRO-seq defined TSSs.
Mic8_1_5_all_TSS_bed <- list()

Mic8_1_5_all_TSS_bed[["all_annotation_liver_TSS_500_with_strand"]] = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_5_ARA_strand/for_paper_Mic8_1_5_all_genes/Mic8_1_5_for_paper_all_genes_bed/all_annotation_liver_TSS_500_with_strand.bed', header = FALSE)

Mic8_1_5_all_TSS_ara_list <- list()

for (i in seq_along(Mic8_1_5_all_TSS_bed)) {
  ara_name <- paste0("ara","_",names(Mic8_1_5_all_TSS_bed)[i])
  print(ara_name)
  
  Mic8_1_5_all_TSS_ara_list[[ara_name]] <- ARA(list(merge_Ctl_MicroC), 
                                       bed = Mic8_1_5_all_TSS_bed[[i]][,1:3], 
                                       strand = Mic8_1_5_all_TSS_bed[[i]][,6],
                                       size_bin = 61)
}

visualise(Mic8_1_5_all_TSS_ara_list[["ara_all_annotation_liver_TSS_500_with_strand"]])

plot_tmp <- visualise(Mic8_1_5_all_TSS_ara_list[["ara_all_annotation_liver_TSS_500_with_strand"]])

library(ggplot2)

ggsave(plot_tmp, file="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_5_ARA_strand/all_TSS_Mic8_1_5_all_genes/ara_ctl_all_genes.png", width = 6, height = 6, dpi = 600)


#Perform ARA analysis for regulated genes
#Read the bed files
Mic8_1_5_bed <- list()

Mic8_1_5_bed[["order_4_at_leastplus_1.585_TSS_500_with_strand"]] = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_5_ARA_strand/Mic8_1_5_TSS_500_bed/order_4_at_leastplus_1.585_TSS_500_with_strand.bed', header = FALSE)

Mic8_1_5_bed[["order_5_at_leastplus_1_TSS_500_with_strand"]] = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_5_ARA_strand/Mic8_1_5_TSS_500_bed/order_5_at_leastplus_1_TSS_500_with_strand.bed', header = FALSE)

Mic8_1_5_bed[["order_6_at_leastplus_0.585_TSS_500_with_strand"]] = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_5_ARA_strand/Mic8_1_5_TSS_500_bed/order_6_at_leastplus_0.585_TSS_500_with_strand.bed', header = FALSE)

Mic8_1_5_bed[["order_9_at_leastminus_0.585_TSS_500_with_strand"]] = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_5_ARA_strand/Mic8_1_5_TSS_500_bed/order_9_at_leastminus_0.585_TSS_500_with_strand.bed', header = FALSE)

Mic8_1_5_bed[["order_10_at_leastminus_1_TSS_500_with_strand"]] = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_5_ARA_strand/Mic8_1_5_TSS_500_bed/order_10_at_leastminus_1_TSS_500_with_strand.bed', header = FALSE)

Mic8_1_5_bed[["order_11_at_leastminus_1.585_TSS_500_with_strand"]] = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_5_ARA_strand/Mic8_1_5_TSS_500_bed/order_11_at_leastminus_1.585_TSS_500_with_strand.bed', header = FALSE)

Mic8_1_5_bed[["True_NONE_TSS_500_with_strand"]] = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_5_ARA_strand/Mic8_1_5_TSS_500_bed/True_NONE_TSS_500_with_strand.bed', header = FALSE)

#perform the ARA analysis and compare insulin and control groups
Mic8_1_5_ara_list <- list()

for (i in seq_along(Mic8_1_5_bed)) {
  ara_name <- paste0("ara","_",names(Mic8_1_5_bed)[i])
  print(ara_name)
  
  Mic8_1_5_ara_list[[ara_name]] <- ARA(list(merge_Ctl_MicroC, merge_Ins_MicroC), 
                                       bed = Mic8_1_5_bed[[i]][,1:3], 
                                       strand = Mic8_1_5_bed[[i]][,6],
                                       size_bin = 61)
}


library(ggplot2)

# plot the graph
# Loop through all the data frames in the list and save each plot
for (i in seq_along(Mic8_1_5_ara_list)) {
  # Get the name of the data frame
  df_name <- names(Mic8_1_5_ara_list)[i]
  
  # Generate the file name for the plot
  plot_file_path <- paste0("/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_5_ARA_strand/Mic8_1_5_plot_2nd_same_y/", df_name, "_same_y.png")
  
  # Generate the plot
  plot_tmp <- visualise(Mic8_1_5_ara_list[[i]], colour_lim_contrast = c(-0.1, 0.12))
  
  ggsave(plot_tmp, file=plot_file_path)
}