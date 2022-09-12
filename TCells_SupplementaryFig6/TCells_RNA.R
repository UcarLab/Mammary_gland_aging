setwd("~/Desktop/Mice_BC/FinalFigurePanels")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape)
library(readxl)
library(purrr)
library(cinaR)
library(cinaRgenesets)
library(tidyr)
library(stringr)
library(scater)
library(Seurat)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(ggpubr)

TCells <- readRDS("./TCells_RNA.rds")
Idents(TCells) <- "NewAnnotations"
TCells <- RunUMAP(TCells, dims = 1:30, reduction = "har30", reduction.name="humap33", reduction.key = "humap3key")
TCells <- subset(TCells, idents = "Doublets", invert = T)
my_levels <- c("Naive CD4 LEF1+","CD4 CCR7+ JUN+","Memory CD4 CD44+", "CD8 CCR7+ ", "CD8 CCR9+", "CD8 GZMK+", "CD8 GZMM +","CD8 CD69+ ISG15+", "IL7R+ ICOS+","NK Cells")
my_levels <- rev(my_levels)
Idents(TCells) <- "NewAnnotations"

pdf("FeaturePlots_ITGAE_PDCD1_CTLA4_LAG3.pdf", height = 10, width = 20)
FeaturePlot(TCells, features = c("Gzmk","Il17a","Itgae","Pdcd1","Gzmm", "Ccl5", "Ctla4", "Lag3"), reduction = "humap33", cols = c("lightgrey", "red3"),label = F, ncol = 4, keep.scale="all", slot = "scale.data")
dev.off()

#GD <- subset(TCells, idents='CD8 GZMK+')


A <- GD@assays$RNA@data
A <- as.data.frame(A)
A1 <- A[rownames(A) %in% c("Pdcd1","Ctla4","Itgae"),]
A1 <- as.data.frame(t(A1))

a = A1[A1$Itgae > 0,]
b = A1[A1$Pdcd1 > 0,]
c = A1[A1$Ctla4 > 0,]

pdf("CorrelationPlot_GZMK.pdf", height = 7, width = 7)
plot(A1$Itgae, A1$Pdcd1, main="Itgae vs Pdcd1",xlab="Itgae expression", ylab="Pdcd1 expression", pch=19)
plot(A1$Itgae, A1$Ctla4, main="Itgae vs Ctla4",xlab="Itgae expression", ylab="Ctla4 expression", pch=19)
dev.off()
TCells@meta.data$NewAnnotations <- factor(x = TCells@meta.data$NewAnnotations, levels = my_levels)

FeaturePlot(TCells, features = c("Ccl3"), reduction = "humap33", cols = c("lightgrey", "red3"),label = F, ncol = 1, keep.scale="all", slot = "scale.data")
TCells_Old <- subset(TCells, idents="mm10_18mths")
TCells_Young <- subset(TCells, idents="mm10_3mths")


p1 <- FeaturePlot(TCells_Old, features = c("Il17a", "Itgae", "Tnf", "Jag1"), reduction = "humap33", cols = c("lightgrey", "red3"),label = F, ncol = 4, keep.scale="all", slot = "scale.data")
p2 <- FeaturePlot(TCells_Young, features = c("Il17a", "Itgae", "Tnf", "Jag1"), reduction = "humap33", cols = c("lightgrey", "red3"),label = F, ncol = 4, keep.scale="all", slot = "scale.data")
p1
p2

DiffGenes <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/Pseudobulk/Combined/IL7R_ICOS.txt", sep = '\t')


##FIGURE 1A
#colr <- c("orangered1","deepskyblue","lightblue2","darkblue","purple","dodgerblue3", "goldenrod2", "firebrick4", "lightpink", "springgreen4")
colr <- c("springgreen4","goldenrod2","darkblue","dodgerblue3","purple","lightblue2", "deepskyblue", "firebrick4","orangered1", "lightpink")
p1 <- DimPlot(TCells, reduction = "humap33", group.by = "NewAnnotations", cols = colr,label = T, label.size = 5)
#pdf("./TCells_RNA_Fig1A_Poster.pdf", height = 6, width = 8)
p1
dev.off()

##FIGURE 1B
Idents(TCells) <- "orig.ident"
OldTCells <- subset(TCells, idents = "mm10_18mths")
YoungTCells <- subset(TCells, idents = "mm10_3mths")
p1 <- DimPlot(OldTCells, reduction = "humap33", group.by = "NewAnnotations", cols = rev(colr),label = F, label.size = 8)
p2 <- DimPlot(YoungTCells, reduction = "humap33", group.by = "NewAnnotations", cols = rev(colr),label = F, label.size = 8)
pdf("./TCells_RNA_Fig1B.pdf", height = 4, width = 6)
p1
p2
dev.off()

##FIGURE 1C
cd_genes <- list('TCells'=c('Cd3e','Cd4','Cd8a'),
                 'Naive T'=c('Sell','Il7r','Lef1', 'Ccr7', "Ccr9"),
                 'Memory T'=c('Cd44', "Itgae"), 
                 'Cytotoxic T'=c( 'Gzma', 'Gzmk', 'Gzmm'),
                 'Activated T'=c('Cd69',  'Icos', 'Jun'),
                 'Interferon'=c('Isg15',"Irf1","Irf7"), 
                 'Gamma Delta'=c("Il17a",'Trdc',"Trgv2"),
                 'NK'=c('Ncr1','Nkg7', "Klrc1"),
                 'Other' = c('Pdcd1','Ccl5','Ctla4','Eomes'))
p1 <- DotPlot(TCells,features = cd_genes, cols = c("lightgrey","red"), group.by = "NewAnnotations", dot.scale = 4, cluster.idents = F) + RotatedAxis() + theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
pdf("./TCells_RNA_Fig1C_All.pdf", height = 8, width = 18)
p1
dev.off()

##FIGURE 1D
Meta <- TCells@meta.data
Meta <- Meta[,c(1,4,29)]
Meta <- Meta %>% group_by_all() %>% summarise(COUNT = n())
Meta$replicate <- gsub("M3","M18",Meta$replicate)
NewMeta = cast(Meta, replicate + NewAnnotations ~ orig.ident, sum) 
NewMeta$Ratio <- log2(NewMeta$mm10_18mths/NewMeta$mm10_3mths)
NewMeta <- NewMeta[,c(1,2,5)]
NewMeta$NewAnnotations <- as.character(NewMeta$NewAnnotations)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep1","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep2","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep3","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep4","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep5","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep6","Base",0)
NewMeta$NewAnnotations <- as.factor(NewMeta$NewAnnotations)
NewMeta$Ratio <- as.numeric(NewMeta$Ratio)

#colr1 <- c("orangered1","deepskyblue","lightblue2","darkblue","purple","dodgerblue3", "goldenrod2", "firebrick4", "lightpink", "springgreen4")
colr <- c("black","springgreen4","goldenrod2","darkblue","dodgerblue3","purple","lightblue2", "deepskyblue", "firebrick4","orangered1", "lightpink")

p1 <- ggplot(NewMeta,aes(x = reorder(NewAnnotations, -Ratio), y = Ratio, fill = NewAnnotations)) + 
  geom_boxplot()+
  scale_fill_manual(values=(colr)) +
  xlab("CellTypes") +
  ylab("log2 FC (Old/Young)") +
  ggtitle("Cell Type Proportions") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  stat_compare_means(method = "anova", label.y = 10)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = "Base", hide.ns = T, label.y = 8)   
p1                                                                            
pdf("./TCells_RNA_Fig1D_Poster_Pvals.pdf", height = 7, width = 8)
p1
dev.off()

##FIGURE 1E
FeaturePlot(TCells, features = c("Ccl3"), reduction = "humap3", cols = c("lightgrey", "red3"),label = F, ncol = 1, keep.scale="all", slot = "scale.data")
#p1 <- FeaturePlot(TCells, features = c("Gzmk","Pdcd1","Ctla4", "Lag3","Cd44","Icos","Il17a", "Ccl5","Itgae","Eomes"), reduction = "humap33", cols = c("lightgrey", "red3"),label = F, ncol = 5, keep.scale="all", slot = "scale.data")
pdf("./TCells_RNA_Fig1E_Poster.pdf", height = 8, width = 15)
p1
dev.off()

##FIGURE 1F
A <- read_xlsx("~/Desktop/Mice_BC/Subclustering_TCells/Counts2.xlsx")
colnames(A) <- c("CellType", "RNA_Closing", "RNA_Opening", "ATAC_Closing", "ATAC_Opening")
A <- A[-1,]

Closing <- A[,c(1,2)]
Opening <- A[,c(1,3)]
colnames(Opening) <- c("CellType","Counts")
colnames(Closing) <- c("CellType","Counts")
Opening$Direction <- "Upregulated"
Closing$Direction <- "Downregulated"
Closing$Counts <- as.numeric(Closing$Counts) * -1
FinalCounts <- rbind(Opening, Closing)
long_DF <- FinalCounts %>% gather("Type", "Count", Counts)
long_DF$CellType <- factor(long_DF$CellType, levels=c("NaiveCD4_LEF1","CD4_CCR7_JUN","MemoryCD4_CD44", "CD8_CCR7", "CD8_CCR9", "CD8_GZMK", "CD8_GZMM","CD8_CD69_ISG15", "IL7R_ICOS","NK"))
long_DF$Count <- as.numeric(long_DF$Count)

pdf("./TCells_RNA_Fig1F.pdf", height = 10, width = 14)
ggplot(long_DF[order(long_DF$Type, decreasing = F),], aes(y=Count, x=CellType, fill = Direction)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("lightgrey","black")) +
  xlab("Gene Counts") +
  theme(axis.text =element_text(size=16),axis.title=element_text(size=16,face="bold")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Cell Types")+
  guides(fill=guide_legend(title="Direction"))
dev.off()

##FIGURE 1G
ImmMod <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/Pseudobulk/Combined/vp2008_hpea-list.csv")
WP <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/Pseudobulk/Combined//wp_hpea-list.csv")
ImmMod$Module <- "ImmuneModules"
WP$Module <- "Wikipathways"

ImmMod <- ImmMod[grep("Inflammation|Cytotoxic cells", ImmMod$module.name), ]
WP <- WP[grep("WP437|WP2704|WP382|WP69|WP585|WP49|WP205|WP3404|WP366|WP2112|WP75|WP1838|WP707|WP2679|WP395",WP$module.name),]

ImmMod <- ImmMod[(ImmMod$adj.p <= 0.1),]
WP <- WP[(WP$adj.p <= 0.1),]
HPEA.table.filtered <- rbind(ImmMod, WP)
df.plot <- HPEA.table.filtered


filter.pathways <- TRUE
fdr.cutoff <- 0.1

if(filter.pathways){
  if (sum(df.plot$adj.p < fdr.cutoff) == 0){
    stop("You can't filter because there are no pathways to be displayed!")
  }
  df.plot <- subset(df.plot, adj.p < fdr.cutoff)
}

df.plot$Status <- factor(df.plot$Status, levels=c("Up", "Down"))
df.plot$contrast <- factor(df.plot$contrast, levels=c("CD4_CCR7_JUN", "NaiveCD4_LEF1","MemoryCD4_CD44","CD8_CCR7","CD8_CCR9","CD8_GZMK","CD8_GZMM","CD8_CD69_ISG15","IL7R_ICOS","NK"))

# create ggplot
color_values <- color_values[c(4,7)]
plot.dot <- ggplot2::ggplot(df.plot,
                            ggplot2::aes(x = contrast,
                                         y = module.name,
                                         size = ifelse(adj.p < fdr.cutoff, -log(adj.p), NA),
                                         color = Status))

plot.dot <- plot.dot + ggplot2::geom_point()

plot.dot <- plot.dot + facet_grid(scales='free', space = "free")
plot.dot <- plot.dot + ggplot2::labs(x = "Contrast",
                                     y = "Pathways",
                                     color = "Sign",
                                     size = "-log10(adj.p)",
                                     caption = paste0("FDR < ", fdr.cutoff))
plot.dot <- plot.dot + ggplot2::scale_color_manual(values = color_values)
plot.dot <- plot.dot + theme(axis.text =element_text(size=13),
                             axis.title=element_text(size=15,face="bold"))
plot.dot <- plot.dot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot.dot <- plot.dot + theme(legend.text = element_text(size=16))
plot.dot <- plot.dot + facet_grid(Module ~ ., space = "free", scales = "free")
pdf("./TCells_RNA_Fig1G2.pdf", height = 8, width = 10)
plot.dot
dev.off()

##FIGURE 1H - CD8 GZMM
normalized_counts <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/Pseudobulk/Cluster6/DifferentialGeneCounts_All_Cluster6.txt", sep = '\t')
sig_genes <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/Pseudobulk/Combined/CD8_GZMM.txt", sep = '\t')

sig_genes_down <- c("Jun","Fos","Jund","Cxcr5","Ccr9","Irf4","Tcf7","Cd27","Lef1")
#sig_genes_down <- c("Ccr7","Ifi214","Klf2","Sell","Lef1","Foxp1","Bach2")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

sig_genes_up <- c("Gzmk","Gzmm","Gzmb","Prf1","Nkg7","S100a4","S100a6","S100a10","S100a13","Ccr2","Cxcr6","Ccr5","Ccl5","Il7r","Ctla2a","Cdk4","Fasl", "Eomes") 
#sig_genes_up <- c("Cxcl2","Cxcl3","Ccr2","Ccr1","Ccl5","Cxcr6","S100a4","S100a6","S100a11","Gzmk","Fosb","Itga1","Lmna","Ctla2a","Klrg1","Il7r","Casp1","Casp4","Lag3","Dusp1","Junb","Pdcd1","Anxa1","Cdk6")
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
ann_colors = list(group_id = c(Young = "lightgray", Old = "black"))

pdf("./TCells_RNA_Fig1H.pdf", height = 6, width = 6)
pheatmap(sig_norm, 
         color = heat_colors, 
         annotation_colors = ann_colors,
         cluster_rows = T, 
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         treeheight_row = 0, treeheight_col = 0,
         annotation_col = df, 
         fontsize_row = 13,
         scale = "row")      
dev.off()

##FIGURE 1I - Memory CD4 CD44
normalized_counts <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/Pseudobulk/Cluster9/DifferentialGeneCounts_All_Cluster9.txt", sep = '\t')
sig_genes <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/Pseudobulk/Combined/MemoryCD4_CD44.txt", sep = '\t')

#sig_genes_down <- c("Jun","Fos","Jund","Cxcr5","Ccr9","Irf4","Tcf7","Cd27","Lef1")
sig_genes_down <- c("Ccr7","Ifi214","Klf2","Sell","Lef1","Foxp1","Bach2")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

#sig_genes_up <- c("Gzmk","Gzmm","Gzmb","Prf1","Nkg7","S100a4","S100a6","S100a10","S100a13","Ccr2","Cxcr6","Ccr5","Ccl5","Il7r","Ctla2a","Cdk4","Fasl", "Eomes") 
sig_genes_up <- c("Cxcl2","Cxcl3","Ccr2","Ccr1","Ccl5","Cxcr6","S100a4","S100a6","S100a11","Gzmk","Fosb","Itga1","Lmna","Ctla2a","Klrg1","Il7r","Casp1","Casp4","Lag3","Dusp1","Junb","Pdcd1","Anxa1","Cdk6")
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
ann_colors = list(group_id = c(Young = "lightgray", Old = "black"))

pdf("./TCells_RNA_Fig1I.pdf", height = 7, width = 7)
pheatmap(sig_norm, 
         color = heat_colors, 
         annotation_colors = ann_colors[1],
         cluster_rows = T, 
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         treeheight_row = 0, treeheight_col = 0,
         annotation_col = df, 
         fontsize_row = 13,
         scale = "row")      
dev.off()



##FIGURE 1E
ImmMod <- read.csv("~/Desktop/Mice_BC/FinalFigurePanels/TCELLS_DG_Type/vp2008_hpea-list.csv")
WP <- read.csv("~/Desktop/Mice_BC/FinalFigurePanels/TCELLS_DG_Type/wp_hpea-list.csv")
ImmMod$Module <- "ImmuneModules"
WP$Module <- "Wikipathways"

ImmMod <- ImmMod[grep("Inflammation|Cytotoxic cells", ImmMod$module.name), ]
WP <- WP[grep("WP437|WP2704|WP382|WP69|WP585|WP49|WP205|WP3404|WP366|WP2112|WP75",WP$module.name),]
ImmMod <- ImmMod[(ImmMod$adj.p <= 0.1),]
WP <- WP[(WP$adj.p <= 0.1),]
HPEA.table.filtered <- rbind(ImmMod, WP)
df.plot <- HPEA.table.filtered


filter.pathways <- TRUE
fdr.cutoff <- 0.1

if(filter.pathways){
  if (sum(df.plot$adj.p < fdr.cutoff) == 0){
    stop("You can't filter because there are no pathways to be displayed!")
  }
  df.plot <- subset(df.plot, adj.p < fdr.cutoff)
}

df.plot$Status <- factor(df.plot$Status, levels=c("Up", "Down"))
#df.plot$contrast <- factor(df.plot$contrast, levels=c("NaiveCD4_LEF1","CD4_CCR7_JUN","MemoryCD4_CD44","CD8_CCR7","CD8_CCR9","CD8_GZMK","CD8_GZMM","IL7R_ICOS","NK"))
df.plot$contrast <- factor(df.plot$contrast, levels=c("Naive CD4 LEF1+","CD4 CCR7+ JUN+","Memory CD4 CD44+","CD8 CCR7+ ","CD8 CCR9+","CD8 GZMK+","CD8 GZMM +","CD8 CD69+ ISG15+","IL7R+ ICOS+","NK Cells"))

# create ggplot
color_values <- color_values[c(4,7)]
plot.dot <- ggplot2::ggplot(df.plot,
                            ggplot2::aes(x = contrast,
                                         y = module.name,
                                         size = ifelse(adj.p < fdr.cutoff, -log(adj.p), NA),
                                         color = Status))

plot.dot <- plot.dot + ggplot2::geom_point()

plot.dot <- plot.dot + facet_grid(scales='free', space = "free")
plot.dot <- plot.dot + ggplot2::labs(x = "Contrast",
                                     y = "Pathways",
                                     color = "Sign",
                                     size = "-log10(adj.p)",
                                     caption = paste0("FDR < ", fdr.cutoff))
plot.dot <- plot.dot + ggplot2::scale_color_manual(values = color_values)
plot.dot <- plot.dot + theme(axis.text =element_text(size=13),
                             axis.title=element_text(size=15,face="bold"))
plot.dot <- plot.dot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot.dot <- plot.dot + theme(legend.text = element_text(size=16))
plot.dot <- plot.dot + facet_grid(Module ~ ., space = "free", scales = "free")
pdf("./TCells_RNA_Enrichemnt_Type.pdf", height = 10, width = 12)
plot.dot
dev.off()






RNA <- TCells_RNA@meta.data
ATAC <- TCells_ATAC@meta.data
RNA <- RNA[,c(4,29)]
ATAC <- ATAC[,c(21,36)]
RNA <- RNA %>% group_by_all() %>% dplyr::summarise(COUNT = n())
ATAC <- ATAC %>% group_by_all() %>% dplyr::summarise(COUNT = n())
A1 <- cast(ATAC, predicted.id ~ dataset, sum)
A2 <- cast(RNA, NewAnnotations ~ replicate, sum)
write.table(A1, "TCell_ATAC_Cellcounts.txt", sep = '\t', row.names = F)
write.table(A2, "TCell_RNA_Cellcounts.txt", sep = '\t', row.names = F)


