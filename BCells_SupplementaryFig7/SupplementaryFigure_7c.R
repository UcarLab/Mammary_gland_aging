setwd("~/path/to/output")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(tidyverse)

BCells_RNA <- readRDS("./BCells_RNA.rds")
Idents(BCells_RNA) <- "SubclusterAnnotations"

my_levels <- c("C1","C2","C3", "C4", "C5", "C6", "C7")
BCells_RNA@meta.data$SubclusterAnnotations <- factor(x = BCells_RNA@meta.data$SubclusterAnnotations, levels = rev(my_levels))

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

pdf("~/path/to/SFig_7c.pdf", height = 12, width = 24)
DotPlot(TCells, features = BCell_markers, cols = c("lightgrey","red"), group.by = "SubclusterAnnotations", dot.scale = 4, cluster.idents = F) + RotatedAxis()
dev.off()
