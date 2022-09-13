setwd("~/path/to/output")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(tidyverse)

BCells_RNA <- readRDS("./BCells_RNA.rds")
Idents(BCells_RNA) <- "SubclusterAnnotations"

my_levels <- c("C1","C2","C3", "C4", "C5", "C6", "C7")
BCells_RNA@meta.data$SubclusterAnnotations <- factor(x = BCells_RNA@meta.data$SubclusterAnnotations, levels = my_levels)
colr <- c("orangered1","firebrick4","purple","lightblue2","dodgerblue3","goldenrod2","darkblue")
pdf("~/path/to/SFig_7b.pdf", height = 6, width = 8)
DimPlot(BCells_RNA, reduction = "humap7", group.by = "Annotation", label = F, label.size = 8, cols = colr) 
dev.off()

