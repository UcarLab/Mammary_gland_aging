setwd("/Users/sharms/Desktop/Mice_BC/FinalFigurePanels//")
library(Seurat)
library(harmony)
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
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(magrittr)
library(Matrix)
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

# TCells <- readRDS("../Subclustering_BCells/BCells_RNA.rds")
# TCells <- FindVariableFeatures(TCells, selection.method = "vst")
# TCells <- RunPCA(TCells, features = VariableFeatures(object = TCells), reduction.name="pca_subset", npcs = 50)
# #pdf("Elbow.pdf", height = 10, width = 10)
# ElbowPlot(TCells, ndims=50, reduction = "pca")
# #dev.off()
# TCells <- RunHarmony(TCells, group.by.vars="expt", reduction = "pca_subset", dims.use = 1:50, reduction.save="har50")
# TCells <- FindNeighbors(TCells, dims = 1:30, reduction="harmony", graph.name = "hknn5")
# TCells <- FindClusters(TCells, resolution = 0.5, graph.name = "hknn5")
# TCells <- RunUMAP(TCells, dims = 1:30, reduction = "harmony", reduction.name="humap5")
# #pdf("HUMAP_PCA50_res0.5_NewPlots.pdf", height = 7, width = 7)
# DimPlot(TCells, reduction = "humap5", label = T)
# DimPlot(TCells, reduction = "humap5", group.by = "expt")
# #dev.off()
# saveRDS(TCells, "./BCells_RNA.rds")

Markers <- FindAllMarkers(TCells, only.pos = T, logfc.threshold = 0.25)

Mouse_BCell_markers <- list(
  'B1'=c("CD138","CD43","CD5"),
  'Follicular'=c('CD45R','CD38','CD19','CD21','CD22','CD23','PAX5'),
  'Marginal'=c('IGM','CD1D','CD9','CD35','EBF','E2A','OCT2'),
  'Transitional'=c('B220','CD93','CD24','BR3','TACI'),
  'Germinal'=c("GL7",'CD95','PNA',"BCL6"),
  'Plasma'=c('BCMA','CD126','CD184','CD320','BLIMP1',"IRF4","XBP1"),
  'Memory'=c('CD80', 'CD73', 'CD273', 'CD84', 'CD86', 'PBF1')
)


# MarkerGenes_BCells <- list('Nonswitched_BCells'=c('CD1A','CD48','PTPRC','CD44','B2M','FCGR2B','CD27','ICAM1','CD84'),
#                            'Transitional_BCells'=c('CD24','MME','CD79B','CD9','CD38','CD72','CD55','BTLA'), #,'TREML','TNFRSF16'
#                            'Switched_BCells'=c('ATP1B3','ITGB1','CD82','CD53','CD99','TNFRSF14','CD180','CR1','MRC1'),
#                            'ABC_DN2_BCells'=c('ITGAX','TBX21','FGR','FCRL5','TFEC','ZEB2',
#                                               'ZBTB32','CD86','AICDA','EBI3','SIGLEC6',
#                                               'PRDM1'), #,'CXCR3'
#                            'MBC_BCells'=c('IGHA1','IGHG1','IGHG3'),
#                            'ISG_BCells'=c('IFI44L','ISG15','IFI6','XAF1','MX1'),
#                            'naive_BCells'=c('CD1D','CXCR5','CCR7','FCER2','CD200'),
#                            'PCs_BCells'=c("MZB1","TNFRSF17"))

rownames(TCells@assays[["RNA"]]@scale.data) <- toupper(rownames(TCells@assays[["RNA"]]@scale.data))
rownames(TCells@assays[["RNA"]]@data) <- toupper(rownames(TCells@assays[["RNA"]]@data))

# BCell_markers <- list('B Lineage'=c("MS4A1", "BANK1","CD79A"),
#                       'Follicular'=c('B220','CD38','CD19','CD21','CD22','CD23','PAX5'),
#                       'Naive'=c('SELL','CCR7','IL4R','CD1D','FCER2','CD200'),
#                       'Memory'=c('CD27','CD24','CD80', 'CD73', 'CD273', 'CD84', 'CD86', 'PBF1'),
#                       'Transitional'=c('TCL1A','CD9','MME','FCRL1','SOX4','CD1A','CD48','PTPRC','CD44','B2M','FCGR2B','ICAM1'),
#                       'ZEB2'=c('ZEB2'),
#                       'ABC/DN2'=c('ITGAX', 'TBX21', 'FGR','TFEC', 'ZBTB32', 'CXCR5', 'TLR7', 'CR2'),
#                       'ISG'=c('IFI27', 'ISG15', 'IFI44', 'STAT1', 'MX1','IFI44L','IFI6','XAF1'),
#                       'Igs'=c('IGHG1', 'IGHD', 'IGHM', 'IGHA1'),
#                       'DCs'=c('CD1C',"CST3","FCER1G"),
#                       "pDCs"=c("IRF7", "TCF4"),
#                       'PCs'=c("MZB1","TNFRSF17",'JCHAIN',"CD126", "CD184", "CD320"))

BCell_markers <- list('B Lineage'=c("MS4A1", "BANK1","CD79A"),
                      'Follicular'=c('CD38','CD19','CR2','CD22','PAX5'),
                      'Naive'=c('SELL','CCR7',"IL4RA",'FCER2A','CD200'),
                      'Memory'=c('CD27','CD80', 'NT5E', 'PDCD1LG2', 'CD84','CD86'),
                      'Transitional'=c('CD9','MME','FCRL1','SOX4','CD48','PTPRC','CD44','B2M','FCGR2B','ICAM1'),
                      'Traficking' = c("CXCR4","CXCR5"),
                      'ABC/DN2'=c('TFEC', 'ZBTB32',"ZEB2", 'TLR7'),
                      'ISG'=c('IFI27', 'ISG15', 'IFI44', 'STAT1', 'MX1','XAF1'),
                      'Igs'=c('IGHG1', 'IGHD', 'IGHM'),
                      'DCs'=c("CST3","FCER1G"),
                      "pDCs"=c("IRF7", "TCF4"),
                      'PCs'=c("MZB1","TNFRSF17",'JCHAIN', "CD320","XBP1","PRDM1","SDC1"))

TCells@meta.data$FinalAnnotation <- as.factor(TCells@meta.data$Annotation)
TCells@meta.data$FinalAnnotation <- factor(TCells@meta.data$FinalAnnotation, levels = c("B_C6","B_C5","B_C4","B_C3","B_C2","B_C1","B_C0"))
pdf("Dotplot_PCA30_res0.5_NewPlots1.pdf", height = 12, width = 24)
DotPlot(TCells, features = BCell_markers, cols = c("lightgrey","red"), group.by = "FinalAnnotation", dot.scale = 4, cluster.idents = F) + RotatedAxis()
dev.off()


BCell_markers_F <- list('B Lineage'=c("MS4A1", "BANK1","CD79A"),
                      'Follicular'=c('CD38','CD19','CR2','CD22','PAX5'),
                      'Naive'=c('SELL','CCR7',"IL4RA",'FCER2A','CD200'),
                      'Memory'=c('PDCD1LG2', 'CD84','CD86'),
                      'Transitional'=c('FCRL1','SOX4','CD48','PTPRC','CD44','B2M','FCGR2B','ICAM1'),
                      'Traficking' = c("CXCR4","CXCR5"),
                      'ABC/DN2'=c("ZEB2", 'TLR7'),
                      'ISG'=c('IFI27', 'ISG15', 'IFI44', 'STAT1', 'MX1','XAF1'),
                      'Igs'=c('IGHG1', 'IGHD', 'IGHM'),
                      'DCs'=c("CST3","FCER1G"),
                      "pDCs"=c("IRF7", "TCF4"),
                      'PCs'=c("MZB1","XBP1"))

pdf("Dotplot_PCA30_res0.5_NewPlots1_F.pdf", height = 12, width = 20)
DotPlot(TCells, features = BCell_markers_F, cols = c("lightgrey","red"), group.by = "FinalAnnotation", dot.scale = 4, cluster.idents = F) + RotatedAxis()
dev.off()



TCells <- readRDS("./BCells_RNA.rds")
#colr <- c("orangered1","dodgerblue3","purple","firebrick4","goldenrod2","lightblue2","darkblue")
colr <- c("orangered1","firebrick4","purple","lightblue2","dodgerblue3","goldenrod2","darkblue")
pdf("HUMAP_PCA50_res0.5_NewPlots_COLOR.pdf", height = 7, width = 7)
DimPlot(TCells, reduction = "humap5", label = T, label.size = 6, cols = colr, group.by = "Annotation")
dev.off()



TCells@meta.data$Annotation <- paste("B_C", TCells@meta.data$seurat_clusters, sep = "")
Meta <- TCells@meta.data
Meta <- Meta[,c(1,4,12)]
Meta <- Meta %>% group_by_all() %>% summarise(COUNT = n())
Meta$replicate <- gsub("M3","M18",Meta$replicate)
NewMeta = cast(Meta, replicate + Annotation ~ orig.ident, sum) 
NewMeta$Ratio <- log2(NewMeta$mm10_18mths/NewMeta$mm10_3mths)
NewMeta <- NewMeta[,c(1,2,5)]
NewMeta$Annotation <- as.character(NewMeta$Annotation)
# NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep1","Base",0)
# NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep2","Base",0)
# NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep3","Base",0)
# NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep4","Base",0)
# NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep5","Base",0)
# NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep6","Base",0)
# NewMeta$Annotation <- as.factor(NewMeta$Annotation)
NewMeta$Ratio <- as.numeric(NewMeta$Ratio)
colr <- c("orangered1","firebrick4","purple","lightblue2","dodgerblue3","goldenrod2","darkblue")

p1 <- ggplot(NewMeta,aes(x = reorder(Annotation, -Ratio), y = Ratio, fill = Annotation)) + 
  geom_boxplot()+
  scale_fill_manual(values=(colr)) +
  xlab("CellTypes") +
  ylab("log2 FC (Old/Young)") +
  ggtitle("Cell Type Proportions") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  ylim(-4,2)
  #stat_compare_means(method = "anova", label.y = 10)+        # Add global annova p-value
  #stat_compare_means(label = "p.signif", method = "t.test",ref.group = "Base", hide.ns = T, label.y = 8)   

#stat1 <- compare_means(Ratio ~ Annotation,  data = NewMeta, ref.group = "Base", method = "t.test")
#write.table(stat1, "./RNA_FoldChange_Stats_BCells.txt", row.names = F, sep = '\t')

p1                                                               
pdf("./BCells_Cellration_FigureB1.pdf", height = 6, width = 8)
p1
dev.off()

BCell_markers <- c("Tyrobp", "Clec12a", "Cd300a", "Cd7", "Chdh", "Mycl")
pdf("./BCells_Cellration_Figure.pdf", height = 6, width = 8)
DotPlot(TCells, features = BCell_markers, cluster.idents = T) + RotatedAxis()
dev.off()
