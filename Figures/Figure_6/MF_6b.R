setwd("~/path/to/output")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(tidyverse)

DM_RNA <- readRDS("path/to/object")
Idents(DM_RNA) <- "SubclusterAnnotations"
my_levels <- c("C1","C2","C3", "C4", "C5", "C6", "C7")
DM_RNA@meta.data$SubclusterAnnotations <- factor(x = DM_RNA@meta.data$SubclusterAnnotations, levels = rev(my_levels))

cd_genes <- list('Monocytes'=c("S100a8","S100a9"),
                 'M1 like' = c("Cd80", "Cd86","Cd38"),
                 'M2 like' = c("Mrc1", "Cd209f", "Cd163"),
                 'mregDC' = c("Ccr7", "Fscn1","Socs2","Relb", "Il4i1"),
                 'cDC1' = c("Clec9a","Xcr1"),
                 'cDC2' = c("Clec10a","Fcer1a"))

pdf("~/path/to/pdf", height = 8, width = 18)
DotPlot(DM_RNA,features = cd_genes, cols = c("lightgrey","red"), group.by = "SubclusterAnnotations", dot.scale = 4, cluster.idents = F, keep.scale="all", slot = "scale.data") + RotatedAxis() + theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
dev.off()

