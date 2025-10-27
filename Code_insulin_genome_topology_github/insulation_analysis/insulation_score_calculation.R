#Insulation and Cross-Point Interaction (CPI) analysis

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

#call the insulation score genome-widely
insulation_1250_window5 <- insulation_score(
  list(merge_Ctl_MicroC, merge_Ins_MicroC),
  window = 5)

insulation_Ctl <- insulation_1250_window5$insula_score[,c("chrom","start","end","Ctl")]
insulation_Ins <- insulation_1250_window5$insula_score[,c("chrom","start","end","Ins")]

#generate bedgraph
write.table(insulation_Ctl, "/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_4_insulation/Mic8_1_4_bedgraph/Ctl_score.bedgraph",row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

write.table(insulation_Ins, "/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_4_insulation/Mic8_1_4_bedgraph/Ins_score.bedgraph",row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

#The bedgraph files were then converted to bigwig files on HMS O2 cluster computer.

