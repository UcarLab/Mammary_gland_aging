setwd("~/path/to/output")
library(Seurat)
library(Signac)
library(ggplot2)

TCells <- readRDS("./TCells_RNA.rds")
Idents(TCells) <- "SubclusterAnnotations"
pdf("~/path/to/SFig6a.pdf", height = 10, width = 20)
FeaturePlot(TCells, features = c("Gzmk","Il17a","Itgae","Pdcd1","Gzmm", "Ccl5", "Ctla4", "Lag3"), reduction = "humap33", cols = c("lightgrey", "red3"),label = F, ncol = 4, keep.scale="all", slot = "scale.data")
dev.off()
