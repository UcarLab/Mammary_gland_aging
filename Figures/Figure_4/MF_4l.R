library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape)
library(tidyr)
library(cowplot)
library(RColorBrewer)
library(ggpubr)

##Fib C0
normalized_counts <- read.csv("~/path/to Fib C0 counts", sep = '\t')
sig_genes <- read.csv("~/path/to/Fib C0 diff genes", sep = '\t')

sig_genes_down <- c("Htra1","Serpinh1","Ptn","Galnt16","Col1a1","Sema3c","Aoc3")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

sig_genes_up <- c("Sept4","Cdkn2a","Cdkn1a","Lsamp","Slc43a3","Crabp1")
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
sig_norm_C0 <- sig_norm

##Fib C2
normalized_counts <- read.csv("~/path/to Fib C2 counts", sep = '\t')
sig_genes <- read.csv("~/path/to/Fib C2 diff genes", sep = '\t')

sig_genes_down <- c("Slc16a1","Sfrp2","Slc36a2","Sdc1","Itm2a","Cxcl10","Gpt2","Fn1","Car3","Meg3")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

sig_genes_up <- c("Gm14964","Cd300lg","Nrxn2","Csprs","Ndufa4l2","Smoc2","Tmem204","Slco2b1","Cdkn1a")
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
sig_norm_C2 <- sig_norm

##Fib C3
normalized_counts <- read.csv("~/path/to Fib C3 counts", sep = '\t')
sig_genes <- read.csv("~/path/to/Fib C3 diff genes", sep = '\t')

sig_genes_down <- c("Clec3b","Olfml2b","Sparc","Itm2a","Basp1","Meg3")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

sig_genes_up <- c("Fabp4","Acta2","Galnt15","Gdf10","mt-Atp6")
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
sig_norm_C3 <- sig_norm

##Fib C4
normalized_counts <- read.csv("~/path/to Fib C4 counts", sep = '\t')
sig_genes <- read.csv("~/path/to/Fib C4 diff genes", sep = '\t')

sig_genes_down <- c("Il1r2", "Ggt5", "Adam23", "Chrdl1", "Ctsh", "Ccdc80","Spon1",'Igsf10',"Cpxm1","Itm2a")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

sig_genes_up <- c("C13002612Rik", "Dmkn", "Hmcn2", "Ntn4", "Fdps", "Il34", "Serpine2", "Duox1","Crip2","Adamtsi3")
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
sig_norm_C4 <- sig_norm

rownames(sig_norm_C0) <- paste(rownames(sig_norm_C0), "_C0", sep  = '')
rownames(sig_norm_C2) <- paste(rownames(sig_norm_C2), "_C2", sep  = '')
rownames(sig_norm_C3) <- paste(rownames(sig_norm_C3), "_C3", sep  = '')
rownames(sig_norm_C4) <- paste(rownames(sig_norm_C4), "_C4", sep  = '')

sig_norm <- rbind(sig_norm_C0, sig_norm_C2, sig_norm_C3, sig_norm_C4)
#sig_norm <- as.data.frame(t(sig_norm))
ann_colors = list(group_id = c(Young = "lightgray", Old = "black"), CellType =c(C0 = "springgreen4", C3 = "purple", C4 = "turquoise", C2 = "goldenrod2"))

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
         gaps_row = c(13,23,34))
dev.off()
