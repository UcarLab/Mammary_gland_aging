library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape)
library(tidyr)
library(cowplot)
library(RColorBrewer)
library(ggpubr)

##Luminal-HS
normalized_counts <- read.csv("~/path/to Luminal-HS counts", sep = '\t')
sig_genes <- read.csv("~/path/to/Luminal-HS diff genes", sep = '\t')

sig_genes_down <- c("Lef1","Mt2","Mt1","Fgg","Sox9","Cd37","Crispld2","Tox2","Gadd45a","Gimap6")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

sig_genes_up <- c("B2m","Epb41l3","Tmprss6","Cpe","Cpeb1","Pygl","Socs2","Igfals","Arhgap26","Stk32a")
sig_genes_up <- as.data.frame(sig_genes_up)
colnames(sig_genes_up) <- "genes"
sig_genes_up$Direction <- "Opening"

sig_genes <- rbind(sig_genes_down, sig_genes_up)
rownames(sig_genes) <- sig_genes$genes

sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_genes$genes)

#sig_norm <- sig_norm[,c(1,8:13,2:7)]
sig_norm <- sig_norm[,c(1,11:13,8:10,5:7,2:4)]

# Run pheatmap using the metadata data frame for the annotation
heat_colors <- rev(brewer.pal(9, "RdBu"))
df <- colnames(sig_norm[,2:length(colnames(sig_norm))])
df <- as.data.frame(df)
colnames(df) <- "group_id"
df$group_id[1:6] <- "Young"
df$group_id[7:12] <- "Old"
rownames(df) <- colnames(sig_norm[,2:length(colnames(sig_norm))])

# df_row <- as.data.frame(sig_genes$Role)
# rownames(df_row) <- rownames(sig_genes)
# colnames(df_row) <- "Role"

rownames(sig_norm) <- sig_norm$gene
sig_norm <- sig_norm[,-1]
sig_norm_Luminal-HS <- sig_norm

##Luminal-AV
normalized_counts <- read.csv("~/path/to Luminal-AV counts", sep = '\t')
sig_genes <- read.csv("~/path/to/Luminal-AV diff genes", sep = '\t')

sig_genes_down <- c("Fndc4","Spry4","Slc38a1","Masp1","Cdhr1","Zdfb2","Sema6a","Col11a1","Crispid2","Areg")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

sig_genes_up <- c("Rnase1","Il18r","Fam189a2","Mal","Aspa","Pde1c","Arg1","Hp","Nkd2","Cryba4")
sig_genes_up <- as.data.frame(sig_genes_up)
colnames(sig_genes_up) <- "genes"
sig_genes_up$Direction <- "Opening"

sig_genes <- rbind(sig_genes_down, sig_genes_up)
rownames(sig_genes) <- sig_genes$genes

sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_genes$genes)

#sig_norm <- sig_norm[,c(1,8:13,2:7)]
sig_norm <- sig_norm[,c(1,11:13,8:10,5:7,2:4)]

# Run pheatmap using the metadata data frame for the annotation
heat_colors <- rev(brewer.pal(9, "RdBu"))
df <- colnames(sig_norm[,2:length(colnames(sig_norm))])
df <- as.data.frame(df)
colnames(df) <- "group_id"
df$group_id[1:6] <- "Young"
df$group_id[7:12] <- "Old"
rownames(df) <- colnames(sig_norm[,2:length(colnames(sig_norm))])

# df_row <- as.data.frame(sig_genes$Role)
# rownames(df_row) <- rownames(sig_genes)
# colnames(df_row) <- "Role"

rownames(sig_norm) <- sig_norm$gene
sig_norm <- sig_norm[,-1]
sig_norm_Luminal-AV <- sig_norm

##Myoepithelial
normalized_counts <- read.csv("~/path/to Myoepithelial counts", sep = '\t')
sig_genes <- read.csv("~/path/to/Myoepithelial diff genes", sep = '\t')

sig_genes_down <- c("Slit2","Wnt6","Anxa3","Crlf1","Nkd2","C530043A13Rik","Nnat","Eln","Col4a2","Ccnd2")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

sig_genes_up <- c("Cp","Bbox1","Tspan8","Prrx1","Sncg","Atp1b1","Sfrp2","Fli1","Snhg11","Brinp3")
sig_genes_up <- as.data.frame(sig_genes_up)
colnames(sig_genes_up) <- "genes"
sig_genes_up$Direction <- "Opening"

sig_genes <- rbind(sig_genes_down, sig_genes_up)
rownames(sig_genes) <- sig_genes$genes

sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_genes$genes)

#sig_norm <- sig_norm[,c(1,8:13,2:7)]
sig_norm <- sig_norm[,c(1,11:13,8:10,5:7,2:4)]

# Run pheatmap using the metadata data frame for the annotation
heat_colors <- rev(brewer.pal(9, "RdBu"))
df <- colnames(sig_norm[,2:length(colnames(sig_norm))])
df <- as.data.frame(df)
colnames(df) <- "group_id"
df$group_id[1:6] <- "Young"
df$group_id[7:12] <- "Old"
rownames(df) <- colnames(sig_norm[,2:length(colnames(sig_norm))])

# df_row <- as.data.frame(sig_genes$Role)
# rownames(df_row) <- rownames(sig_genes)
# colnames(df_row) <- "Role"

rownames(sig_norm) <- sig_norm$gene
sig_norm <- sig_norm[,-1]
sig_norm_Myoepithelial <- sig_norm


rownames(sig_norm_Luminal-HS) <- paste(rownames(sig_norm_Luminal-HS), "_Luminal-HS", sep  = '')
rownames(sig_norm_Luminal-AV) <- paste(rownames(sig_norm_Luminal-AV), "_Luminal-AV", sep  = '')
rownames(sig_norm_Myoepithelial) <- paste(rownames(sig_norm_Myoepithelial), "_Myoepithelial", sep  = '')

sig_norm <- rbind(sig_norm_Luminal-AV, sig_norm_Luminal-HS, sig_norm_Myoepithelial)
#sig_norm <- as.data.frame(t(sig_norm))
ann_colors = list(group_id = c(Young = "lightgray", Old = "black"), CellType =c(Luminal-HS = "springgreen3", Myoepithelial = "purple", Luminal-HS = "deepskyblue"))

genes <- as.data.frame(rownames(sig_norm))
genes <- data.frame(do.call('rbind', strsplit(as.character(genes$`rownames(sig_norm)`),'_',fixed=TRUE)))
colnames(genes) <- c("genes","CellType")
rownames(genes) <- rownames(sig_norm)
Celltype <- as.data.frame(genes$CellType)
Genes <- as.data.frame(genes$genes)
rownames(Genes) <- rownames(sig_norm)
rownames(Celltype) <- rownames(sig_norm)
colnames(Genes) <- "Genes"
colnames(Celltype) <- "CellType"


pdf("~/path/to/pdf", height = 15, width = 5)
pheatmap(sig_norm,
         color = heat_colors,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         treeheight_row = 0, treeheight_col = 0,
         annotation_col = df,
         annotation_row = Celltype,
         labels_row = Genes$Genes,
         fontsize_col = 13,
         scale = "row",
         gaps_col = 6,
         gaps_row = c(20,40))
dev.off()
