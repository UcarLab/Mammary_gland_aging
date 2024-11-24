setwd("~/path/to/output")
library(Seurat)
library(Signac)
library(ggplot2)

DM_RNA <- readRDS("./DM_RNA.rds")
Idents(TCells) <- "SubclusterAnnotations"
pdf("~/path/to/MFig_5f.pdf", height = 15, width = 8)
FeaturePlot(TCells, features = c("Mrc1","Cd163","Fscn1","Mmp12","Xcr1", "Itgax"), reduction = "humap7", cols = c("lightgrey", "red3"),label = F, ncol = 2, keep.scale="all", slot = "scale.data")
dev.off()
