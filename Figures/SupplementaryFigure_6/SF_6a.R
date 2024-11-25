setwd("~/path/to/output")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(tidyverse)

TCells <- readRDS("path/to/object")
##Confirm cell identities are stored in the right metadata column
Idents(TCells) <- "NewAnnotations"
my_levels <- c("C1","C2","C3", "C4", "C5", "C6", "C7")
TCells@meta.data$SubclusterAnnotations <- factor(x = TCells@meta.data$NewAnnotations, levels = rev(my_levels))

cd_genes <- list('TCell'=c("Cd3e","Cd4","Cd8a"),
                 'Naive and Memory' = c("Sell", "Il7r","Lef1","Ccr7","Ccr9"),
                 'Cytotoxic' = c("Gzma", "Gzmk", "Gzmm","Ccl5"),
                 'Checkpoint' = c("Pdcd1","Ctla4"),
                 'Interferon' = c("Isg15","Irf1","Irf7"),
                 'GD T' = c("Il17ra","Trdc","Trgv2"),
                 'NK' = c("Nrc1","Nkg7","Klrc1"),
                 'Other' = c("Eomes","Icos","Itgae"))

pdf("~/path/to/pdf", height = 8, width = 18)
DotPlot(TCells,features = cd_genes, cols = c("lightgrey","red"), group.by = "SubclusterAnnotations", dot.scale = 4, cluster.idents = F, keep.scale="all", slot = "scale.data") + RotatedAxis() + theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
dev.off()

