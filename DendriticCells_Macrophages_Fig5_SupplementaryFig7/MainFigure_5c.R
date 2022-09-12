setwd("~/path/to/output")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(tidyverse)

DM_ATAC <- readRDS("./DM_ATAC.rds")
Idents(DM_ATAC) <- "SubclusterAnnotations"

my_levels <- c("C1","C2","C3", "C4", "C5", "C6", "C7")
DM_ATAC@meta.data$SubclusterAnnotations <- factor(x = DM_ATAC@meta.data$SubclusterAnnotations, levels = my_levels)
colr <- c("orangered1","dodgerblue3","purple","firebrick4","goldenrod2","lightblue2","darkblue")
pdf("~/path/to/MFig_5c.pdf", height = 6, width = 8)
DimPlot(DM_ATAC, reduction = "humap7", group.by = "Annotation", label = F, label.size = 8, cols = colr) 
dev.off()

