###anchor annotation and E-P loop identification

###Step1: read loop files from mustache
loop_r1250 <- list()

loop_r1250[["Ctl_loop"]] <- read.delim("/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_4_loop_mustache/mustach_output/reso_1250/merge_Ctl_MicroC_r1250_loop.tsv", header=TRUE, sep="\t")

loop_r1250[["Ins_loop"]] <- read.delim("/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_4_loop_mustache/mustach_output/reso_1250/merge_Ins_MicroC_r1250_loop.tsv", header=TRUE, sep="\t")

View(loop_r1250[["Ctl_loop"]])


###Step2: Merge control and insulin loops
library(dplyr)

#generate loop datasets that only contain coordinates
loop_r1250[["Ctl_loop_coord_only"]] <- loop_r1250[["Ctl_loop"]][,1:6]
loop_r1250[["Ins_loop_coord_only"]] <- loop_r1250[["Ins_loop"]][,1:6]

loop_r1250[["loop_combined_coord_only"]] <- bind_rows(loop_r1250[["Ctl_loop_coord_only"]], loop_r1250[["Ins_loop_coord_only"]]) %>%
  distinct()

View(loop_r1250[["loop_combined_coord_only"]])


###Step3: Name every loop
loop_r1250[["loop_combined_coord_only"]]$Loop_Name <- paste("loop", seq_len(nrow(loop_r1250[["loop_combined_coord_only"]])), sep="_")

#Generate bed files for both anchor 1 and anchor 2
write.table(loop_r1250[["loop_combined_coord_only"]][,c(1,2,3,7)],"/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_4_loop_mustache/Mic8_4_3_loop_annotation/Mic8_4_3_R/anchors_loop_combined/loop_r1250_combined_anchor1.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)

write.table(loop_r1250[["loop_combined_coord_only"]][,c(4,5,6,7)],"/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_4_loop_mustache/Mic8_4_3_loop_annotation/Mic8_4_3_R/anchors_loop_combined/loop_r1250_combined_anchor2.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)


###Step4: Assign whether an anchor is an enhancer or promoter.
#This step was performed on Shell
#see anchor_enhancer_promoter_shell.txt

anchor_annotation <- list()
#anchor1
anchor_annotation[["anchor1_TSS_enh_annot"]] <- read.delim("/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_4_loop_mustache/Mic8_4_3_loop_annotation/Mic8_4_3_R/anchors_loop_combined_annotated/anchor1_with_TSS_enhancers.tsv", header = FALSE)

colnames(anchor_annotation[["anchor1_TSS_enh_annot"]]) <- c("bin_chr","bin_start","bin_end","loop_name","chr","start","end","anchor_annotation")

View(anchor_annotation[["anchor1_TSS_enh_annot"]])

#anchor2
anchor_annotation[["anchor2_TSS_enh_annot"]] <- read.delim("/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_4_loop_mustache/Mic8_4_3_loop_annotation/Mic8_4_3_R/anchors_loop_combined_annotated/anchor2_with_TSS_enhancers.tsv", header = FALSE)

colnames(anchor_annotation[["anchor2_TSS_enh_annot"]]) <- c("bin_chr","bin_start","bin_end","loop_name","chr","start","end","anchor_annotation")

View(anchor_annotation[["anchor2_TSS_enh_annot"]])

#Generate the annotated loop file
loop_combined_anchor_annotated <- loop_r1250[["loop_combined_coord_only"]]

loop_combined_anchor_annotated$anchor1_annot <- anchor_annotation[["anchor1_TSS_enh_annot"]]$anchor_annotation[match(loop_combined_anchor_annotated$Loop_Name, anchor_annotation[["anchor1_TSS_enh_annot"]]$loop_name)]

loop_combined_anchor_annotated$anchor2_annot <- anchor_annotation[["anchor2_TSS_enh_annot"]]$anchor_annotation[match(loop_combined_anchor_annotated$Loop_Name, anchor_annotation[["anchor2_TSS_enh_annot"]]$loop_name)]


###Step5: Filter out enhancer-promoter pairs.
#get loops that both anchors are annotated
loop_annotated_analysis <- list()

library(dplyr)

loop_annotated_analysis[["loop_combined_anchor_both_annotated"]] <- loop_combined_anchor_annotated %>%
  filter(!is.na(anchor1_annot) & !is.na(anchor2_annot))

View(loop_annotated_analysis[["loop_combined_anchor_both_annotated"]])

#P-P loop
loop_annotated_analysis[["loop_P_P"]] <- loop_annotated_analysis[["loop_combined_anchor_both_annotated"]][grepl("^ENSMUS", loop_annotated_analysis[["loop_combined_anchor_both_annotated"]]$anchor1_annot) & grepl("^ENSMUS", loop_annotated_analysis[["loop_combined_anchor_both_annotated"]]$anchor2_annot), ]

#E-E loop
loop_annotated_analysis[["loop_E_E"]] <- loop_annotated_analysis[["loop_combined_anchor_both_annotated"]][grepl("^enh", loop_annotated_analysis[["loop_combined_anchor_both_annotated"]]$anchor1_annot) & grepl("^enh", loop_annotated_analysis[["loop_combined_anchor_both_annotated"]]$anchor2_annot), ]

#E-P loop
loop_annotated_analysis[["loop_E_P"]] <- loop_annotated_analysis[["loop_combined_anchor_both_annotated"]][(grepl("^ENSMUS", loop_annotated_analysis[["loop_combined_anchor_both_annotated"]]$anchor1_annot) & grepl("^enh", loop_annotated_analysis[["loop_combined_anchor_both_annotated"]]$anchor2_annot)) |
                    (grepl("^enh", loop_annotated_analysis[["loop_combined_anchor_both_annotated"]]$anchor1_annot) & grepl("^ENSMUS", loop_annotated_analysis[["loop_combined_anchor_both_annotated"]]$anchor2_annot)), ]


#####Step6: analyze E-P loop. Add the fold change of enhancer activities and gene regulation.
loop_E_P_analysis <- loop_annotated_analysis[["loop_E_P"]]

library(stringr)

loop_E_P_analysis <- loop_E_P_analysis %>%
  mutate(enh_id = if_else(str_starts(anchor1_annot, "enh"), anchor1_annot, anchor2_annot))

loop_E_P_analysis <- loop_E_P_analysis %>%
  mutate(gene_id = if_else(str_starts(anchor1_annot, "ENSM"), anchor1_annot, anchor2_annot))

#Read relative activity of enhancers and relative expression of genes.
enh_gene_expression <- list()

enh_gene_expression[["Mic4_9_2_enh_analysis"]] <- read.delim("/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_4_loop_mustache/Mic8_4_3_loop_annotation/enh_gene_expression_Mic4_9_2_Mic4_4/Mic4_9_2_enh_analysis.tsv", header=TRUE, sep="\t")

enh_gene_expression[["Mic4_4_annotation_proseq"]] <- read.delim("/n/data1/bidmc/medicine/kahn/lab/xyl/Mic_project_others_from_2024/Mic8_4_loop_mustache/Mic8_4_3_loop_annotation/enh_gene_expression_Mic4_9_2_Mic4_4/Mic4_4_annotation_proseq.tsv", header=TRUE, sep="\t")

##Add relative expression levels to loop_E_P_analysis
loop_E_P_analysis$enh_log2_rela <- enh_gene_expression[["Mic4_9_2_enh_analysis"]]$log2_rela_FnR[match(loop_E_P_analysis$enh_id, enh_gene_expression[["Mic4_9_2_enh_analysis"]]$name)]

loop_E_P_analysis$gene_log2_rela <- enh_gene_expression[["Mic4_4_annotation_proseq"]]$log2FC[match(loop_E_P_analysis$gene_id, enh_gene_expression[["Mic4_4_annotation_proseq"]]$gene_id)]


##Draw relative expression of paired enhancers and promoters identified by E-P loop analysis 
library(ggplot2)

ggplot(loop_E_P_analysis, aes(x = enh_log2_rela, y = gene_log2_rela)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", bins = 30, contour = TRUE) +
  scale_fill_distiller(palette = 'RdYlBu') +
  labs(title = "", x = "log2(enhancer_activity)", y = "log2(gene_expression") +
  theme_minimal()

#linear regression
model <- lm(enh_log2_rela ~ gene_log2_rela, data = loop_E_P_analysis)
summary(model)
