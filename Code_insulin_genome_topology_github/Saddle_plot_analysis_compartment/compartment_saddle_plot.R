#Saddle Analysis for compartment interaction

library(GENOVA)

#Input .hic files (resolution 50k)
merge_Ctl_MicroC_50k <- load_contacts(signal_path = '/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_2_POSSUMM_Nova/merged_hic_files/with_chr_merged_hic_POSSUMM/merge_Ctl_MicroC.hic',
                                  sample_name = "Ctl",
                                  resolution = 50000,
                                  balancing = 'KR', # this is the default
                                  colour = "blue")

merge_Ins_MicroC_50k <- load_contacts(signal_path = '/n/data1/bidmc/medicine/kahn/lab/xyl/Mic5_2_2_POSSUMM_Nova/merged_hic_files/with_chr_merged_hic_POSSUMM/merge_Ins_MicroC.hic',
                                  sample_name = "Ins",
                                  resolution = 50000,
                                  balancing = 'KR', # this is the default
                                  colour = "orange")

#it needs H3K27Ac ChIP bed files to call compartment A/B
H3K27acPeaks = read.delim('/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_1_loop_TAD_MicroC/R_Mic8_1/Mic8_1_3_compt_low_res/H3K27Ac_ChIPseq_bed/66697_peaks.bed',
                          header = FALSE)

#Compute compartment scores
CS_out = compartment_score(list(merge_Ctl_MicroC_50k, merge_Ins_MicroC_50k),
                           bed = H3K27acPeaks)

#Saddle analysis
saddle_out = saddle(list(merge_Ctl_MicroC_50k, merge_Ins_MicroC_50k),
                    CS_discovery = CS_out,
                    bins = 50)

visualise(saddle_out)

#compartment interaction strength analysis
CSS <- quantify(saddle_out)

compared <- tidyr::spread(unique(CSS[,-c(3,4)]), key = 'exp', value = 'strength')

with(compared, plot(Ctl, Ins, xlim = c(0,4), ylim = c(0,4), pch = 20))
abline(a = 0, b = 1, lty = 1)

