#Differential analysis of enhancer activities in insulin and control groups.

ctl_F_norm_enhancer <- read.csv("/Users/liux/Documents/Lab_project_experiments/Mic4_9_2_Dif_enhancer_activity/Dif_enhancer_R/bwAvebed_output/ctl_F._norm_enhancer.tab", sep = "\t", header = FALSE)

ctl_R_norm_enhancer <- read.csv("/Users/liux/Documents/Lab_project_experiments/Mic4_9_2_Dif_enhancer_activity/Dif_enhancer_R/bwAvebed_output/ctl_R._norm_enhancer.tab", sep = "\t", header = FALSE)

ins_F_norm_enhancer <- read.csv("/Users/liux/Documents/Lab_project_experiments/Mic4_9_2_Dif_enhancer_activity/Dif_enhancer_R/bwAvebed_output/ins_F._norm_enhancer.tab", sep = "\t", header = FALSE)

ins_R_norm_enhancer <- read.csv("/Users/liux/Documents/Lab_project_experiments/Mic4_9_2_Dif_enhancer_activity/Dif_enhancer_R/bwAvebed_output/ins_R._norm_enhancer.tab", sep = "\t", header = FALSE)

colnames(ctl_F_norm_enhancer) <- c("name","size","covered","sum","mean0","mean")
colnames(ctl_R_norm_enhancer) <- c("name","size","covered","sum","mean0","mean")
colnames(ins_F_norm_enhancer) <- c("name","size","covered","sum","mean0","mean")
colnames(ins_R_norm_enhancer) <- c("name","size","covered","sum","mean0","mean")

#mean0 stands for the average over bases with non-covered bases counting as zeroes
enh_analysis <- data.frame(ctl_F_norm_enhancer$name)
colnames(enh_analysis) <- "name"

enh_analysis$ctl_F_mean0 <- ctl_F_norm_enhancer$mean0[match(enh_analysis$name, ctl_F_norm_enhancer$name)]
enh_analysis$ctl_R_mean0 <- ctl_R_norm_enhancer$mean0[match(enh_analysis$name, ctl_R_norm_enhancer$name)]

enh_analysis$ins_F_mean0 <- ins_F_norm_enhancer$mean0[match(enh_analysis$name, ins_F_norm_enhancer$name)]
enh_analysis$ins_R_mean0 <- ins_R_norm_enhancer$mean0[match(enh_analysis$name, ins_R_norm_enhancer$name)]

enh_analysis$log2_rela_FnR <- log2((enh_analysis$ins_F_mean0 + enh_analysis$ins_R_mean0 + 0.01)/(enh_analysis$ctl_F_mean0 + enh_analysis$ctl_R_mean0 + 0.01))

#plot the distribution of enh_analysis$log2_rela_FnR
hist(enh_analysis$log2_rela_FnR, breaks = seq(floor(min(enh_analysis$log2_rela_FnR)), ceiling(max(enh_analysis$log2_rela_FnR)), by = 0.25), ylab="Count", col="grey", xlab="log2_relative_transcription", main = NA)


####Add enhancer location
enhancer_peaks_location_liver <- read.csv("/Users/liux/Documents/Lab_project_experiments/Mic4_9_2_Dif_enhancer_activity/name_added_true_enhancer_peaks_liver.bed", sep = "\t", header = FALSE)

colnames(enhancer_peaks_location_liver) <- c("chr","start","end","name")

enh_analysis$chr <- enhancer_peaks_location_liver$chr[match(enh_analysis$name, enhancer_peaks_location_liver$name)]
enh_analysis$start <- enhancer_peaks_location_liver$start[match(enh_analysis$name, enhancer_peaks_location_liver$name)]
enh_analysis$end <- enhancer_peaks_location_liver$end[match(enh_analysis$name, enhancer_peaks_location_liver$name)]


#write bedgraph
write.table(enh_analysis[,c("chr","start","end","log2_rela_FnR")],"/Users/liux/Documents/Lab_project_experiments/Mic4_9_2_Dif_enhancer_activity/Dif_enhancer_R/bedgraph/enh_rel_activity.bedgraph", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

#bedgraph files were then converted to bigwig files on HMS O2 cluster computer.