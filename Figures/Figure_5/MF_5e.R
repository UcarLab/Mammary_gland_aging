library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape)
library(tidyr)
library(cowplot)
library(RColorBrewer)
library(ggpubr)

##CD8 GZMM
normalized_counts <- read.csv("~/path/to CD8 GZMM counts", sep = '\t')
sig_genes <- read.csv("~/path/to/CD8 GZMM diff genes", sep = '\t')

sig_genes_down <- c("Jun","Fos","Jund","Cxcr5","Ccr9","Irf4","Tcf7","Cd27","Lef1")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

sig_genes_up <- c("Gzmk","Gzmm","Gzmb","Prf1","Nkg7","S100a4","S100a6","S100a10","S100a13","Ccr2","Cxcr6","Ccr5","Ccl5","Il7r","Ctla2a","Cdk4","Fasl", "Eomes") 
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
sig_norm_Gzmm <- sig_norm

##CD8 GZMK
normalized_counts <- read.csv("~/path/to CD8 GZMK counts", sep = '\t')
sig_genes <- read.csv("~/path/to/CD8 GZMK diff genes", sep = '\t')

sig_genes_down <- c("Bach2","Foxp1","Ccr7","Sell","Jun","Ifi209","Ifi213","Ctla2b")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

sig_genes_up <- c("Casp1","S100a4","S100a6","Ccl5","Ccl3","Gzmk","Ifng","Nkg7","Fasl","Socs3","Dusp5")
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
sig_norm_Gzmk <- sig_norm

##Memory CD4
normalized_counts <- read.csv("~/path/to Mem CD4 counts", sep = '\t')
sig_genes <- read.csv("~/path/to/Mem CD4 diff genes", sep = '\t')

sig_genes_down <- c("Jun","Fos","Jund","Cxcr5","Ccr9","Irf4","Tcf7","Cd27","Lef1")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

sig_genes_up <- c("Gzmk","Gzmm","Gzmb","Prf1","Nkg7","S100a4","S100a6","S100a10","S100a13","Ccr2","Cxcr6","Ccr5","Ccl5","Il7r","Ctla2a","Cdk4","Fasl", "Eomes") 
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
sig_norm_MemCD4 <- sig_norm

##Gamma Delta
normalized_counts <- read.csv("~/path/to GD counts", sep = '\t')
sig_genes <- read.csv("~/path/to/GD diff genes", sep = '\t')

sig_genes_down <- c("Lgals3", "Il4", "Icos", "Cd27", "Cd28", "Dusp10")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

sig_genes_up <- c("Hspa1b", "Jag1", "Klf9", "Cxcr3", "Tbx21", "Large1", "Serpina3g", "Ccl5","Pglyrp1","Asb2")
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
sig_norm_GD <- sig_norm

rownames(sig_norm_Cd44) <- paste(rownames(sig_norm_MemCD4), "_Cd44", sep  = '')
rownames(sig_norm_Gzmk) <- paste(rownames(sig_norm_Gzmk), "_Gzmk", sep  = '')
rownames(sig_norm_Gzmm) <- paste(rownames(sig_norm_Gzmm), "_Gzmm", sep  = '')
rownames(sig_norm_GD) <- paste(rownames(sig_norm_GD), "_GD", sep  = '')

sig_norm <- rbind(sig_norm_MemCD4, sig_norm_Gzmk, sig_norm_Gzmm, sig_norm_GD)
#sig_norm <- as.data.frame(t(sig_norm))
ann_colors = list(group_id = c(Young = "lightgray", Old = "black"), CellType =c(Cd44 = "firebrick4", Gzmk = "purple", Gzmm = "dodgerblue3", GD = "goldenrod2"))

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


pdf("~/path/to/pdf", height = 20, width = 5)
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
         gaps_row = c(31,50,77))
dev.off()
