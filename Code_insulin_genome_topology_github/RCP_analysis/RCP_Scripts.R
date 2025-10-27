#Relative Contact Probability (RCP) analysis

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


#RCP for specific regions
#Read the bed files.
order_1_at_leastplus_3 = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed/order_1_at_leastplus_3.bed', header = FALSE)

order_2_at_leastplus_2.58 = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed/order_2_at_leastplus_2.58.bed', header = FALSE)

order_3_at_leastplus_2 = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed/order_3_at_leastplus_2.bed', header = FALSE)

order_4_at_leastplus_1.585 = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed/order_4_at_leastplus_1.585.bed', header = FALSE)

order_5_at_leastplus_1 = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed/order_5_at_leastplus_1.bed', header = FALSE)

order_6_at_leastplus_0.585 = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed/order_6_at_leastplus_0.585.bed', header = FALSE)

order_9_at_leastminus_0.585 = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed/order_9_at_leastminus_0.585.bed', header = FALSE)

order_10_at_leastminus_1 = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed/order_10_at_leastminus_1.bed', header = FALSE)

order_11_at_leastminus_1.585 = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed/order_11_at_leastminus_1.585.bed', header = FALSE)

order_12_at_leastminus_2 = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed/order_12_at_leastminus_2.bed', header = FALSE)

order_13_at_leastminus_2.58 = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed/order_13_at_leastminus_2.58.bed', header = FALSE)

order_14_at_leastminus_3 = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed/order_14_at_leastminus_3.bed', header = FALSE)

True_NONE_reg = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_3_comp_deeptools/bed_files/ins_reg_bed/True_NONE_reg.bed', header = FALSE)


RCP_out = RCP(list(merge_Ctl_MicroC, merge_Ins_MicroC),
              bedlist = list("order_1_at_leastplus_3" = order_1_at_leastplus_3,
                             "order_2_at_leastplus_2.58" = order_2_at_leastplus_2.58,
                             "order_3_at_leastplus_2" = order_3_at_leastplus_2,
                             "order_4_at_leastplus_1.585" = order_4_at_leastplus_1.585,
                             "order_5_at_leastplus_1" = order_5_at_leastplus_1,
                             "order_6_at_leastplus_0.585" = order_6_at_leastplus_0.585,
                             "order_9_at_leastminus_0.585" = order_9_at_leastminus_0.585,
                             "order_10_at_leastminus_1" = order_10_at_leastminus_1,
                             "order_11_at_leastminus_1.585" = order_11_at_leastminus_1.585,
                             "order_12_at_leastminus_2" = order_12_at_leastminus_2,
                             "order_13_at_leastminus_2.58" = order_13_at_leastminus_2.58,
                             "order_14_at_leastminus_3" = order_14_at_leastminus_3,
                             "True_NONE_reg" = True_NONE_reg))


#Save the data frame
RCP_out_smooth_df <- RCP_out$smooth

write.table(RCP_out_smooth_df, file="/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_2_RCP_ins_reg/RCP_out_smooth_df.csv", row.names = FALSE, col.names = TRUE, sep= "\t", quote = FALSE)