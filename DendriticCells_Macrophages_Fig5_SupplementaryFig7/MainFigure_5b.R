setwd("~/path/to/output")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(tidyverse)

DM_RNA <- readRDS("./DM_RNA.rds")
Idents(DM_RNA) <- "SubclusterAnnotations"
my_levels <- c("C1","C2","C3", "C4", "C5", "C6", "C7")
DM_RNA@meta.data$SubclusterAnnotations <- factor(x = DM_RNA@meta.data$SubclusterAnnotations, levels = my_levels)

cd_genes <- list('Monocytes'=c("Cd14","S100a8","S100a9"),
                 'Macrophages'=c("Itgax" ,"Csf1r", "Fcgr3", "Adgre1", "Ms4a7", "Cd68", "H2-Ab1", "Fcgr1", "C1qa"),
                 'M1 Macrophages' = c("Mmp12", "Mmp13", "Spic","Cd80", "Cd86","Cd38"),
                 'M2 Macrophages' = c("Mrc1", "Cd209f", "Cd163","Cd83","Egr2","Cd9","Cd74","Bcl2","Arg1"),
                 'mregDC' = c("Ccr7", "Fscn1","Socs2","Relb", "Il4i1", "Cd40", "Cd274"),
                 'Dendritic Cells' = c("Cd24a", "Cd209a", "Traf1", "Napsa", "Flt3", "Fcgr2b", "Itgam"),
                 'cDC1' = c("Cadm1", "Clec9a","Xcr1","Naaa"))

pdf("~/path/to/MFig_5b.pdf", height = 8, width = 18)
DotPlot(DM_RNA,features = cd_genes, cols = c("lightgrey","red"), group.by = "SubclusterAnnotations", dot.scale = 4, cluster.idents = F, keep.scale="all", slot = "scale.data") + RotatedAxis() + theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
dev.off()
