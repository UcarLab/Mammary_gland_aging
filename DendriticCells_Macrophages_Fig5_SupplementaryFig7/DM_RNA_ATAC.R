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

DM_RNA <- readRDS("./DM_RNA.rds")
DM_RNA@meta.data$Annotation <- gsub("Ma Macrophages", "M2 Macrophages", DM_RNA@meta.data$Annotation)
DM_RNA@meta.data$Annotation <- gsub("Mb Macrophages", "M1 Macrophages", DM_RNA@meta.data$Annotation)
for (i in 1:nrow(DM_RNA@meta.data)) {
  if (DM_RNA@meta.data$Annotation[i] == "IL4I1+ TRAF1+ FLT3+ CCR7+ CCL22+") {
    DM_RNA@meta.data$Annotation[i] = "mregDC"
  }
}
Idents(DM_RNA) <- "Annotation"
my_levels <- c("Activated Monocytes","M1 Macrophages","M2 Macrophages", "mregDC", "FCGR1+ NAPSA+", "DC", "cDC1")
my_levels <- rev(my_levels)
DM_RNA@meta.data$Annotation <- factor(x = DM_RNA@meta.data$Annotation, levels = my_levels)

DM_ATAC <- readRDS("./DM_ATAC.rds")
DM_ATAC@meta.data$predicted.id <- gsub("Ma Macrophages", "M2 Macrophages", DM_ATAC@meta.data$predicted.id)
DM_ATAC@meta.data$predicted.id <- gsub("Mb Macrophages", "M1 Macrophages", DM_ATAC@meta.data$predicted.id)
for (i in 1:nrow(DM_ATAC@meta.data)) {
  if (DM_ATAC@meta.data$predicted.id[i] == "IL4I1+ TRAF1+ FLT3+ CCR7+ CCL22+") {
    DM_ATAC@meta.data$predicted.id[i] = "mregDC"
  }
}
Idents(DM_ATAC) <- "predicted.id"
my_levels <- c("Activated Monocytes","M1 Macrophages","M2 Macrophages", "mregDC", "FCGR1+ NAPSA+", "DC", "cDC1")
DM_ATAC@meta.data$predicted.id <- factor(x = DM_ATAC@meta.data$predicted.id, levels = my_levels)

##FIGURE 1A - RNA UMAP
colr <- c("orangered1","dodgerblue3","purple","firebrick4","goldenrod2","lightblue2","darkblue")
colr <- rev(colr)
p2 <- DimPlot(DM_RNA, reduction = "humap7", group.by = "Annotation", label = F, label.size = 8, cols = colr) 
p2
pdf("./Figure.pdf", height = 6, width = 8)
p1
dev.off()

##FIGURE 1B - RNA AGE UMAP
Idents(DM_RNA) <- "orig.ident"
OldDM_RNA <- subset(DM_RNA, idents = "mm10_18mths")
YoungDM_RNA <- subset(DM_RNA, idents = "mm10_3mths")
p1 <- DimPlot(OldDM_RNA, reduction = "humap7", group.by = "Annotation", cols = colr,label = F, label.size = 8)
p2 <- DimPlot(YoungDM_RNA, reduction = "humap7", group.by = "Annotation", cols = colr,label = F, label.size = 8)
pdf("./DM_RNA_ATAC_Figure1B.pdf", height = 4, width = 6)
p1
p2
dev.off()

##FIGURE 1C - ATAC UMAP
colr <- c("orangered1","dodgerblue3","purple","firebrick4","goldenrod2","lightblue2","darkblue")
p1 <- DimPlot(DM_ATAC, reduction = "umap5", group.by = "predicted.id", label = F, label.size = 8, cols = colr) 
pdf("./DM_RNA_ATAC_Figure1C.pdf", height = 6, width = 8)
p1
dev.off()

##FIGURE 1D - ATAC AGE UMAP
Idents(DM_ATAC) <- "age"
OldDM_ATAC <- subset(DM_ATAC, idents = "18M")
YoungDM_ATAC <- subset(DM_ATAC, idents = "3M")
p1 <- DimPlot(OldDM_ATAC, reduction = "umap5", group.by = "predicted.id", cols = colr,label = F, label.size = 8)
p2 <- DimPlot(YoungDM_ATAC, reduction = "umap5", group.by = "predicted.id", cols = colr,label = F, label.size = 8) 
pdf("./DM_RNA_ATAC_Figure1D.pdf", height = 4, width = 6)
p1
p2
dev.off()

##FIGURE 1E
cd_genes <- list('Monocytes'=c("Cd14","S100a8","S100a9"),
                 'Macrophages'=c("Itgax" ,"Csf1r", "Fcgr3", "Adgre1", "Ms4a7", "Cd68", "H2-Ab1", "Fcgr1", "C1qa"),
                 'M1 Macrophages' = c("Mmp12", "Mmp13", "Spic","Cd80", "Cd86","Cd38"),
                 'M2 Macrophages' = c("Mrc1", "Cd209f", "Cd163","Cd83","Egr2","Cd9","Cd74","Bcl2","Arg1"),
                 'mregDC' = c("Ccr7", "Fscn1","Socs2","Relb", "Il4i1", "Cd40", "Cd274"),
                 'Dendritic Cells' = c("Cd24a", "Cd209a", "Traf1", "Napsa", "Flt3", "Fcgr2b", "Itgam"),
                 'cDC1' = c("Cadm1", "Clec9a","Xcr1","Naaa"))
                 #'cDC2' = c("Fcer1a"),
                 #'preDC' = c("Cd5"),
                 #'Others' = c("Thbd", "Il3ra", "Siglec1"))
#Markers <- c("Il4i1", "Cd274", "Ido1", "Cd14", "Aif1", "Itgax" ,"Csf1r", "Fcgr3", "Adgre1", "Ms4a7", "Cd68", "H2-Ab1", "Fcgr1", "C1qa", "C1qb", "C1qc", "S100a8", "S100a9", "Mrc1", "Cd209f", "Cd163", "Mmp12", "Mmp13", "Spic", "Cd24a", "Cd209a", "Traf1", "Napsa", "Flt3", "Fcgr2b", "Itgam", "Cadm1", "Clec9a","Xcr1", "Fcer1a", "Cd5", "Thbd", "Il3ra", "Siglec1","Cd80", "Cd86", "Ccr7", "Fscn1","Socs2","Relb","Naaa")
#DM_RNA@meta.data$Annotation <- factor(x = DM_RNA@meta.data$Annotation, levels = rev(my_levels))
p1 <- DotPlot(DM_RNA,features = cd_genes, cols = c("lightgrey","red"), group.by = "Annotation", dot.scale = 4, cluster.idents = F) + RotatedAxis() + theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
p1
pdf("./DM_RNA_ATAC_Fig1E_Filtered.pdf", height = 8, width = 18)
p1
dev.off()

pdf("./PosterFeaturePlots_DM.pdf", height = 12, width = 15)
FeaturePlot(DM_RNA, features = c("Mmp12","Pdcd1","Mmp13", "Cxcl2", "Ccl2","Cd163"), reduction = "humap7", cols = c("lightgrey", "red3"),label = F, ncol = 3, keep.scale="all", slot = "scale.data")
FeaturePlot(DM_RNA, features = c("Itgax" ,"Csf1r", "Fcgr3", "Adgre1", "Ms4a7", "Cd68", "H2-Ab1", "Fcgr1", "C1qa"), reduction = "humap7", cols = c("lightgrey", "red3"),label = F, ncol = 3, keep.scale="all", slot = "scale.data")
FeaturePlot(DM_RNA, features = c("Mmp12", "Mmp13", "Spic","Cd80", "Cd86","Cd38"), reduction = "humap7", cols = c("lightgrey", "red3"),label = F, ncol = 3, keep.scale="all", slot = "scale.data")
FeaturePlot(DM_RNA, features = c("Mrc1", "Cd209f", "Cd163","Cd83","Egr2"), reduction = "humap7", cols = c("lightgrey", "red3"),label = F, ncol = 3, keep.scale="all", slot = "scale.data")
FeaturePlot(DM_RNA, features = c("Ccr7", "Fscn1","Socs2","Relb", "Il4i1", "Cd40", "Cd274"), reduction = "humap7", cols = c("lightgrey", "red3"),label = F, ncol = 3, keep.scale="all", slot = "scale.data")
FeaturePlot(DM_RNA, features = c("Cd24a", "Cd209a", "Traf1", "Napsa", "Flt3", "Fcgr2b", "Itgam"), reduction = "humap7", cols = c("lightgrey", "red3"),label = F, ncol = 3, keep.scale="all", slot = "scale.data")
FeaturePlot(DM_RNA, features = c("Cadm1", "Clec9a","Xcr1","Naaa"), reduction = "humap7", cols = c("lightgrey", "red3"),label = F, ncol = 3, keep.scale="all", slot = "scale.data")
FeaturePlot(DM_RNA, features = c("Ccr2", "Csf1r", "Marco", "Pdl2", "Cd40", "Ccl2", "Csf1", "Cd16" , "Pdgfb"), reduction = "humap7", cols = c("lightgrey", "red3"),label = F, ncol = 3, keep.scale="all", slot = "scale.data")
dev.off()


pdf("./FigureFeaturePlots_DM.pdf", height = 10, width = 15)
FeaturePlot(DM_RNA, features = c("Mrc1", "Xcr1", "Fscn1", "Mmp12", "Itgax", "Cd163"), reduction = "humap7", cols = c("lightgrey", "red3"),label = F, ncol = 3, keep.scale="all", slot = "scale.data")
dev.off()



A <- FindAllMarkers(DM_RNA, assay = "RNA")
A <- A[A$avg_log2FC >= 0,]
B <- A %>%
  dplyr::group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = F) %>%
  as.data.frame() 
#write.table(A,"AllMarkers_DM.txt", sep = '\t')
A <- read.csv("../AllMarkers_DM.txt", sep = '\t')
pdf("./PosterFeaturePlots_DM_Select.pdf", height = 5, width = 15)
FeaturePlot(DM_RNA, features = c("Mrc1","Xcr1","Fscn1"), reduction = "humap7", cols = c("lightgrey", "red3"),label = F, ncol = 3, keep.scale="all", slot = "scale.data")
dev.off()

  ##FIGURE 1F
#DM_RNA@meta.data$Annotation <- factor(x = DM_RNA@meta.data$Annotation, levels = rev(my_levels))
Meta <- DM_RNA@meta.data
Meta <- Meta[,c(1,4,31)]
Meta <- Meta %>% group_by_all() %>% summarise(COUNT = n())
Meta$replicate <- gsub("M3","M18",Meta$replicate)
NewMeta = cast(Meta, replicate + Annotation ~ orig.ident, sum) 
NewMeta$Ratio <- log2(NewMeta$mm10_18mths/NewMeta$mm10_3mths)
NewMeta <- NewMeta[,c(1,2,5)]
NewMeta$Annotation <- as.character(NewMeta$Annotation)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep1","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep2","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep3","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep4","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep5","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep6","Base",0)
NewMeta$Annotation <- as.factor(NewMeta$Annotation)
NewMeta$Ratio <- as.numeric(NewMeta$Ratio)
colr <- c("black","orangered1","dodgerblue3","purple","firebrick4","goldenrod2","lightblue2","darkblue")

p1 <- ggplot(NewMeta,aes(x = reorder(Annotation, -Ratio), y = Ratio, fill = Annotation)) + 
  geom_boxplot()+
  scale_fill_manual(values=(colr)) +
  xlab("CellTypes") +
  ylab("log2 FC (Old/Young)") +
  ggtitle("Cell Type Proportions") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  stat_compare_means(method = "anova", label.y = 10)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = "Base", hide.ns = T, label.y = 8)   

stat1 <- compare_means(Ratio ~ Annotation,  data = NewMeta, ref.group = "Base", method = "t.test")
write.table(stat1, "./RNA_FoldChange_Stats_DM.txt", row.names = F, sep = '\t')

p1                                                               
pdf("./DM_RNA_ATAC_Figure1F.pdf", height = 6, width = 8)
p1
dev.off()

##FIGURE 1G
Meta <- DM_ATAC@meta.data
Meta <- Meta[,c(22,23,36)]
Meta <- Meta %>% group_by_all() %>% summarise(COUNT = n())
#Meta$replicate <- gsub("M3","M18",Meta$replicate)
NewMeta = cast(Meta, rep + predicted.id ~ age, sum) 
NewMeta$Ratio <- log2(NewMeta$`18M`/NewMeta$`3M`)
NewMeta <- NewMeta[,c(1,2,5)]
NewMeta$predicted.id <- as.character(NewMeta$predicted.id)
NewMeta[nrow(NewMeta) + 1,] <- c("rep4","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("rep5","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("rep6","Base",0)
NewMeta$predicted.id <- as.factor(NewMeta$predicted.id)
NewMeta$Ratio <- as.numeric(NewMeta$Ratio)
colr <- c("black","orangered1","dodgerblue3","purple","firebrick4","goldenrod2","lightblue2","darkblue")

p1 <- ggplot(NewMeta,aes(x = reorder(predicted.id, -Ratio), y = Ratio, fill = predicted.id)) + 
  geom_boxplot()+
  scale_fill_manual(values=colr) +
  xlab("CellTypes") +
  ylab("log2 FC (Old/Young)") +
  ggtitle("Cell Type Proportions") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  stat_compare_means(method = "anova", label.y = 10)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = "Base", hide.ns = T, label.y = 8)   

stat1 <- compare_means(Ratio ~ predicted.id,  data = NewMeta, ref.group = "Base", method = "t.test")
write.table(stat1, "./ATAC_FoldChange_Stats_DM.txt", row.names = F, sep = '\t')

p1                                                               
pdf("./DM_RNA_ATAC_Figure1G.pdf", height = 6, width = 8)
p1
dev.off()

##FIGURE 1E
# p1 <- FeaturePlot(TCells, features = c("Gzmk","Pdcd1","Ctla4", "Lag3","Cd44","Icos","Il17a", "Ccl5","Ccr5","Eomes"), reduction = "humap33", cols = c("lightgrey", "red3"),label = F, ncol = 3, keep.scale="all", slot = "scale.data")
# pdf("./TCells_RNA_Fig1E.pdf", height = 12, width = 12)
# p1
# dev.off()

##FIGURE 1F
A <- read_xlsx("~/Desktop/Mice_BC/Subclustering_DM/ATAC/DiffPeaks_Type/Counts2.xlsx")
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
ImmMod <- read.csv("~/Desktop/Mice_BC/FinalFigurePanels/DC_DG_Type/vp2008_hpea-list.csv")
WP <- read.csv("~/Desktop/Mice_BC/FinalFigurePanels/DC_DG_Type/wp_hpea-list.csv")
ImmMod$Module <- "ImmuneModules"
WP$Module <- "Wikipathways"

ImmMod <- ImmMod[grep("Inflammation I", ImmMod$module.name), ]
ImmMod <- ImmMod[!(ImmMod$module.name == "Inflammation II"),]
WP <- WP[grep("WP619|WP75|WP615|WP382|WP1839|WP585|WP1836|WP1984|WP205|WP2703|WP530",WP$module.name),]
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
#df.plot$contrast <- gsub("IL4I1_TRAF1_FLT3_CCR7_CCL22", "mregDC", df.plot$contrast)
#df.plot$contrast <- factor(df.plot$contrast, levels=c("CD4_CCR7_JUN", "NaiveCD4_LEF1","MemoryCD4_CD44","CD8_CCR7","CD8_CCR9","CD8_GZMK","CD8_GZMM","CD8_CD69_ISG15","IL7R_ICOS","NK"))
df.plot$contrast <- factor(df.plot$contrast, levels=c("ActivatedMonocytes", "M1_Macrophages","M2_Macrophages","mregDC","FCGR1_NAPSA","DC","cDC1"))

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
pdf("./DM_RNA_ATAC_Figure1H_RNA.pdf", height = 6, width = 10)
plot.dot
dev.off()

##FIGURE 1I - Memory CD4 CD44
normalized_counts <- read.csv("~/Desktop/Mice_BC/Subclustering_DM/Pseudobulk/Cluster3/DifferentialGeneCounts_All_Cluster3.txt", sep = '\t')
sig_genes <- read.csv("~/Desktop/Mice_BC/Subclustering_DM/Pseudobulk/Combined/DC.txt", sep = '\t')
sig_genes$Direction <- "Opening"
for (i in 1:nrow(sig_genes)) {
  if (sig_genes$log2FoldChange[i] < 0) {
    sig_genes$Direction[i] <- "Closing"
  }
}
sig_genes <- sig_genes[,c(1,3)]
colnames(sig_genes) <- c("genes", "Direction")
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

pdf("./DM_RNA_ATAC_Figure1I.pdf", height = 9, width = 8)
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


##Figure 1J - MOtifs
Motifs <- read.csv("~/Desktop/Mice_BC/Subclustering_DM/ATAC/DiffPeaks_Type/Chromvar/TotalMotifs.txt", sep = '\t', header = F)
colnames(Motifs) <- c("ID", "p_val", "log_FC_Change", "pct.1", "pct.2", "adj_p_val", "CellType")
Motifs$CellType <- gsub(".txt","",Motifs$CellType)
Motifs$log_FC_Change <- as.numeric(Motifs$log_FC_Change)
MotifNames <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/ATAC/DiffPeaks/Motifs/Final_Total_Motifs1.txt", sep = '\t', header = F)
MotifNames <- MotifNames[-1,]
MotifNames <- MotifNames[,c(2,3)]
colnames(MotifNames) <- c("ID", "MotifName")
MotifNames <- unique(MotifNames)
Motifs <- merge(MotifNames, Motifs)
Motifs2 <- Motifs
Motifs <- Motifs[(Motifs$log_FC_Change >= 7.5 | Motifs$log_FC_Change <= -7.5),]
Motifs <- Motifs[(Motifs$adj_p_val <= 0.05),]

Motifs1 <- Motifs %>% group_by(CellType) %>% top_n(25, abs(log_FC_Change))
MotifsList <- as.data.frame(unique(Motifs1$MotifName))
colnames(MotifsList) <- "MotifName"

FinalMotifs <- merge(Motifs2, MotifsList)
FinalMotifs <- FinalMotifs[,c(1,4,7,8)]
FinalMotifs$Significance <- "*"

for (i in 1:nrow(FinalMotifs)) {
  if (FinalMotifs$adj_p_val[i] > 0.05) {
    FinalMotifs$Significance[i] = ""
  }
}
Clustering <- FinalMotifs[,c(1,2,4)]
Clustering1 <- cast(Clustering, MotifName ~ CellType, mean, value = 'log_FC_Change')
rownames(Clustering1) <- Clustering1$MotifName
Clustering1 <- Clustering1[,-1]
data <- scale(t(Clustering1))

ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
ord
FinalMotifs$CellType <- factor(FinalMotifs$CellType, levels = c("ActMono","FCGR1_NAPSA","M2_Macrophages","FCGR1_NAPSA","DC","cDC1","mregDC"))
RS <- FinalMotifs[,c(1,5)]
RS1 <- RS %>% group_by_all() %>% summarise(COUNT = n())
RS1 <- RS1[!(RS1$Significance=="" & RS1$COUNT==7),]
RS1 <- as.data.frame(unique(RS1$MotifName))
colnames(RS1) <- "MotifName"

FinalMotifs <- merge(RS1, FinalMotifs)

pdf("./DM_RNA_ATAC_Figure1J.pdf", height = 12, width = 12)
ggplot(FinalMotifs, aes(CellType, MotifName, fill= log_FC_Change)) +
  geom_tile() +
  scale_fill_distiller(palette = "PuOr", ) +
  geom_text(aes(label=Significance), color="black", size=4) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##FIGURE 1K
A <- read_xlsx("~/Desktop/Mice_BC/Subclustering_DM/ATAC/DiffPeaks_Type/Chipseeker_Annots/Counts2.xlsx")
colnames(A) <- c("CellType", "ATAC_Closing", "ATAC_Opening")

Closing <- A[,c(1,2)]
Opening <- A[,c(1,3)]
colnames(Opening) <- c("CellType","Counts")
colnames(Closing) <- c("CellType","Counts")
Opening$Direction <- "Opening"
Closing$Direction <- "Closing"
Closing$Counts <- as.numeric(Closing$Counts) * -1
FinalCounts <- rbind(Opening, Closing)
long_DF <- FinalCounts %>% gather("Type", "Count", Counts)
long_DF$CellType <- factor(long_DF$CellType, levels=c("Activated Monocytes","M1 Macrophages","M2 Macrophages", "mregDC", "FCGR1+ NAPSA+", "DC", "cDC1"))
long_DF$Count <- as.numeric(long_DF$Count)

pdf("./DM_RNA_ATAC_Figure1K.pdf", height = 10, width = 10)
ggplot(long_DF[order(long_DF$Type, decreasing = F),], aes(y=Count, x=CellType, fill = Direction)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("lightgrey","black")) +
  xlab("Gene Counts") +
  theme(axis.text =element_text(size=16),axis.title=element_text(size=16,face="bold")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Cell Types")+
  guides(fill=guide_legend(title="Direction"))
dev.off()

