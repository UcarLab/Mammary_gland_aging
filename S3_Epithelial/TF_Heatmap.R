##10172022
##Making 
##
Myo <- read.csv("C:/Users/angarb/Box/Aging mammary gland paper/scRNA-seq/Pseudobulk_DE_analysis/DESeq2/pairwise/Myoepithelial_mm10_18mths_vs_mm10_3mths_all_genes.csv")
AV <- read.csv("C:/Users/angarb/Box/Aging mammary gland paper/scRNA-seq/Pseudobulk_DE_analysis/DESeq2/pairwise/Luminal-AV_mm10_18mths_vs_mm10_3mths_all_genes.csv")
HS <- read.csv("C:/Users/angarb/Box/Aging mammary gland paper/scRNA-seq/Pseudobulk_DE_analysis/DESeq2/pairwise/Luminal-HS_mm10_18mths_vs_mm10_3mths_all_genes.csv")

TF_List <- read.delim("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure2 epithelial/DA/TF_list.txt")

library(dplyr)
Myo_Subset <- subset(Myo,Myo$gene %in% TF_List$Gene)
AV_Subset <- subset(AV,AV$gene %in% TF_List$Gene)
HS_Subset <- subset(HS,HS$gene %in% TF_List$Gene)

TF <- merge(Myo_Subset[,c("gene","log2FoldChange")], HS_Subset[,c("gene","log2FoldChange")], by = "gene", suffixes = c(".Myo",".HS"), all=T)
TF_all <- merge(TF, AV_Subset[,c("gene","log2FoldChange")], by = "gene", all=T)
rownames(TF_all) <- TF_all[,1]
TF_all  <- TF_all[,-1]
TF_all <- mutate_all(TF_all, function(x) as.numeric(as.character(x)))

TF_pdj <- merge(Myo_Subset[,c("gene","padj")], HS_Subset[,c("gene","padj")], by = "gene", suffixes = c(".Myo",".HS"), all=T)
TF_padj_all <- merge(TF_pdj, AV_Subset[,c("gene","padj")], by = "gene", all=T)

library(dplyr)
TF_all[is.na(TF_all)] <- 0

library(pheatmap)
pheatmap(TF_all,
         cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames=T,
         fontsize = 8, na_col = "grey90")

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure2 epithelial/DA/TF_Epi.pdf", height=10, width=10)

myBreaks <- c(seq(min(TF_all), 0, length.out=ceiling(50/2) + 1), 
              seq(max(TF_all)/50, max(TF_all), length.out=floor(50/2)))
my.colors <- c(colorRampPalette(c("blue", "white","red")))

pheatmap(TF_all,
         cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames=T,
         fontsize = 8, color = colorRampPalette(c("navy", "white", "red"))(50), na_col = "grey90", breaks = myBreaks)
dev.off()

################

LumHS_Positive <- read.delim("C:/Users/angarb/Box/Aging mammary gland paper/scRNA-seq/singlecell_DE_analysis/DE_sc_log0_cutoff/Luminal-HS_18_vs_3_DE_pos_log0_LR.txt")
LumHS_Negative <- read.delim("C:/Users/angarb/Box/Aging mammary gland paper/scRNA-seq/singlecell_DE_analysis/DE_sc_log0_cutoff/Luminal-HS_18_vs_3_DE_neg_log0_LR.txt")

LumHS_SC <- rbind(LumHS_Positive, LumHS_Negative)

LumAV_Positive <- read.delim("C:/Users/angarb/Box/Aging mammary gland paper/scRNA-seq/singlecell_DE_analysis/DE_sc_log0_cutoff/Luminal-Av_18_vs_3_DE_pos_log0_LR.txt")
LumAV_Negative <- read.delim("C:/Users/angarb/Box/Aging mammary gland paper/scRNA-seq/singlecell_DE_analysis/DE_sc_log0_cutoff/Luminal-AV_18_vs_3_DE_neg_log0_LR.txt")

LumAV_SC <- rbind(LumAV_Positive, LumAV_Negative)

Myo_Positive <- read.delim("C:/Users/angarb/Box/Aging mammary gland paper/scRNA-seq/singlecell_DE_analysis/DE_sc_log0_cutoff/Myoepithelial_18_vs_3_DE_pos_log0_LR.txt")
Myo_Negative <- read.delim("C:/Users/angarb/Box/Aging mammary gland paper/scRNA-seq/singlecell_DE_analysis/DE_sc_log0_cutoff/Myoepithelial_18_vs_3_DE_neg_log0_LR.txt")

Myo_SC <- rbind(Myo_Positive, Myo_Negative)

Myo_Subset_SC <- subset(Myo_SC,Myo_SC$gene %in% TF_List$Gene)
AV_Subset_SC <- subset(LumAV_SC,LumAV_SC$gene %in% TF_List$Gene)
HS_Subset_SC <- subset(LumHS_SC,LumHS_SC$gene %in% TF_List$Gene)

TF_SC <- merge(Myo_Subset_SC[,c("gene","avg_log2FC")], HS_Subset_SC[,c("gene","avg_log2FC")], by = "gene", suffixes = c(".Myo",".HS"), all=T)
TF_all_SC <- merge(TF_SC, AV_Subset_SC[,c("gene","avg_log2FC")], by = "gene", all=T)
rownames(TF_all_SC) <- TF_all_SC[,1]
TF_all_SC  <- TF_all_SC[,-1]
TF_all_SC <- mutate_all(TF_all_SC, function(x) as.numeric(as.character(x)))
TF_all_SC [is.na(TF_all_SC )] <- 0

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure2 epithelial/DA/TF_Epi_SC.pdf", height=10, width=10)
myBreaks <- c(seq(min(TF_all_SC), 0, length.out=ceiling(50/2) + 1), 
              seq(max(TF_all_SC)/50, max(TF_all_SC), length.out=floor(50/2)))
my.colors <- c(colorRampPalette(c("blue", "white","red")))
pheatmap(TF_all_SC,
         cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames=T,
         fontsize = 8, color = colorRampPalette(c("navy", "white", "red"))(50), na_col = "grey90", breaks = myBreaks)
dev.off()

TF_pdj_SC <- merge(Myo_Subset_SC[,c("gene","p_val_adj")], HS_Subset_SC[,c("gene","p_val_adj")], by = "gene", suffixes = c(".Myo",".HS"), all=T)
TF_padj_all_SC <- merge(TF_pdj_SC, AV_Subset_SC[,c("gene","p_val_adj")], by = "gene", all=T)
