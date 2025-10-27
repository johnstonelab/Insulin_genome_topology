##Gene distribution analysis in Fig. S2

#Generate the density of TSSs and regulated TSSs on the genome
#input chromosome sizes
mm10_chrom_ordered_sizes <- read.csv("/Users/liux/Documents/Lab_project_experiments/Gm4-30_HiC_WT_CTCFKO_liver/mm10_UCSC_chr_size/ordered/mm10.chrom.ordered.sizes", sep = "\t", header = FALSE)
colnames(mm10_chrom_ordered_sizes) <- c("chr","length")

#Write a function to calculate gene or tss density. Inputs are: TSS_file, bin_size, chrom_size_file, output_path
#The outputs should be: annotation_proseq_bin_density, and a bedgraph file saved at the output_path

create_bin_density_and_bedgraph <- function(file, bin_size, chrom_size_file, output_path) {
  
  #initiate
  bin_reso <- file[,c("gene_id","TSS","chr")]
  
  #assign the TSSs to the bins of genome
  bin_reso$bin_start <- ((bin_reso$TSS - 1) %/% bin_size) * bin_size + 1
  bin_reso$bin_end <- bin_reso$bin_start + bin_size - 1
  
  #merge the bins of genome (generate a column whose name is numbers)
  library(dplyr)
  bin_reso_density <- bin_reso %>%
    group_by(chr, bin_start, bin_end) %>%
    summarise(numbers = n(), .groups = 'drop')
  
  #density calculation
  bin_reso_density$density <- bin_reso_density$numbers/(bin_size/1000000)
  
  #add chrom size column
  bin_reso_density$chr_size <- chrom_size_file$length[match(bin_reso_density$chr,chrom_size_file$chr)]
  
  #adjust bin_end
  bin_reso_density <- bin_reso_density %>%
    mutate(bin_end = ifelse(bin_end > chr_size, chr_size, bin_end))
  
  #make scientific number as FALSE
  bin_reso_density$bin_start <- format(bin_reso_density$bin_start, scientific = FALSE)
  bin_reso_density$bin_end <- format(bin_reso_density$bin_end, scientific = FALSE)
  
  #Write bedgraph
  print(output_path)
  write.table(bin_reso_density[,c("chr","bin_start","bin_end","density")], output_path, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  
  #output the data frame
  return(bin_reso_density[,c("chr","bin_start","bin_end","numbers","density")])
}

 
#Use the function on different resolution
all_PROseq_TSS_bin_1mb_density <- create_bin_density_and_bedgraph(annotation_proseq, 1000000, mm10_chrom_ordered_sizes, "/Users/liux/Documents/Lab_project_experiments/Mic4_4_RNA_PRO_cor/Mic4_4_RNA_PRO_cor/Mic4-7_density_of_genes/all_PROseq_TSS_bin_1mb_density.bedgraph")

all_PROseq_TSS_bin_0.5mb_density <- create_bin_density_and_bedgraph(annotation_proseq, 500000, mm10_chrom_ordered_sizes, "/Users/liux/Documents/Lab_project_experiments/Mic4_4_RNA_PRO_cor/Mic4_4_RNA_PRO_cor/Mic4-7_density_of_genes/all_PROseq_TSS_bin_0.5mb_density.bedgraph")

all_PROseq_TSS_bin_0.1mb_density <- create_bin_density_and_bedgraph(annotation_proseq, 100000, mm10_chrom_ordered_sizes, "/Users/liux/Documents/Lab_project_experiments/Mic4_4_RNA_PRO_cor/Mic4_4_RNA_PRO_cor/Mic4-7_density_of_genes/all_PROseq_TSS_bin_0.1mb_density.bedgraph")

all_PROseq_TSS_bin_0.05mb_density <- create_bin_density_and_bedgraph(annotation_proseq, 50000, mm10_chrom_ordered_sizes, "/Users/liux/Documents/Lab_project_experiments/Mic4_4_RNA_PRO_cor/Mic4_4_RNA_PRO_cor/Mic4-7_density_of_genes/all_PROseq_TSS_bin_0.05mb_density.bedgraph")

all_PROseq_TSS_bin_0.01mb_density <- create_bin_density_and_bedgraph(annotation_proseq, 10000, mm10_chrom_ordered_sizes, "/Users/liux/Documents/Lab_project_experiments/Mic4_4_RNA_PRO_cor/Mic4_4_RNA_PRO_cor/Mic4-7_density_of_genes/all_PROseq_TSS_bin_0.01mb_density.bedgraph")


#Distribution analysis on differential expressed genes defined by PRO-seq
#annotation_proseq[which(annotation_proseq$pvalue<0.05 & abs(annotation_proseq$log2FC) >= 0.321),]

sig_PROseq_TSS_bin_1mb_density <- create_bin_density_and_bedgraph(annotation_proseq[which(annotation_proseq$pvalue<0.05 & abs(annotation_proseq$log2FC) >= 0.321),], 1000000, mm10_chrom_ordered_sizes, "/Users/liux/Documents/Lab_project_experiments/Mic4_4_RNA_PRO_cor/Mic4_4_RNA_PRO_cor/Mic4-7_density_of_genes/sig_PROseq_TSS_bin_1mb_density.bedgraph")

sig_PROseq_TSS_bin_0.5mb_density <- create_bin_density_and_bedgraph(annotation_proseq[which(annotation_proseq$pvalue<0.05 & abs(annotation_proseq$log2FC) >= 0.321),], 500000, mm10_chrom_ordered_sizes, "/Users/liux/Documents/Lab_project_experiments/Mic4_4_RNA_PRO_cor/Mic4_4_RNA_PRO_cor/Mic4-7_density_of_genes/sig_PROseq_TSS_bin_0.5mb_density.bedgraph")

sig_PROseq_TSS_bin_0.1mb_density <- create_bin_density_and_bedgraph(annotation_proseq[which(annotation_proseq$pvalue<0.05 & abs(annotation_proseq$log2FC) >= 0.321),], 100000, mm10_chrom_ordered_sizes, "/Users/liux/Documents/Lab_project_experiments/Mic4_4_RNA_PRO_cor/Mic4_4_RNA_PRO_cor/Mic4-7_density_of_genes/sig_PROseq_TSS_bin_0.1mb_density.bedgraph")

sig_PROseq_TSS_bin_0.05mb_density <- create_bin_density_and_bedgraph(annotation_proseq[which(annotation_proseq$pvalue<0.05 & abs(annotation_proseq$log2FC) >= 0.321),], 50000, mm10_chrom_ordered_sizes, "/Users/liux/Documents/Lab_project_experiments/Mic4_4_RNA_PRO_cor/Mic4_4_RNA_PRO_cor/Mic4-7_density_of_genes/sig_PROseq_TSS_bin_0.05mb_density.bedgraph")

sig_PROseq_TSS_bin_0.01mb_density <- create_bin_density_and_bedgraph(annotation_proseq[which(annotation_proseq$pvalue<0.05 & abs(annotation_proseq$log2FC) >= 0.321),], 10000, mm10_chrom_ordered_sizes, "/Users/liux/Documents/Lab_project_experiments/Mic4_4_RNA_PRO_cor/Mic4_4_RNA_PRO_cor/Mic4-7_density_of_genes/sig_PROseq_TSS_bin_0.01mb_density.bedgraph")

#bedgraph files were converted to bigwig files on HMS O2 cluster computer.