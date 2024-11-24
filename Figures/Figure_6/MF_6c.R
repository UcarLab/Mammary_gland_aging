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
p1 <- DimPlot(DM_ATAC, reduction = "humap7", group.by = "Annotation", label = F, label.size = 8, cols = colr) 

Meta <- DM_ATAC@meta.data
Meta <- Meta[,c(22,23,36)]
Meta <- Meta %>% group_by_all() %>% summarise(COUNT = n())
NewMeta = cast(Meta, rep + SubclusterAnnotations ~ age, sum) 
NewMeta$Ratio <- log2(NewMeta$`18M`/NewMeta$`3M`)
NewMeta <- NewMeta[,c(1,2,5)]
NewMeta$SubclusterAnnotations <- as.character(NewMeta$SubclusterAnnotations)
NewMeta$Ratio <- as.numeric(NewMeta$Ratio)

p2 <- ggplot(NewMeta,aes(x = reorder(SubclusterAnnotations, -Ratio), y = Ratio, fill = SubclusterAnnotations)) + 
  geom_boxplot()+
  scale_fill_manual(values=colr) +
  xlab("CellTypes") +
  ylab("log2 FC (Old/Young)") +
  ggtitle("Cell Type Proportions") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  geom_hline(yintercept=0, linetype="dashed", color = "black")
  
## For Calculating P Values, we added a zero base ratio to refer to. Then we used the compare_means function from ggpubr to calculate the p values. We added the p values to the plots manually later on.
NewMeta[nrow(NewMeta) + 1,] <- c("rep4","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("rep5","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("rep6","Base",0)
stat <- compare_means(Ratio ~ SubclusterAnnotations,  data = NewMeta, ref.group = "Base", method = "t.test")
write.table(stat, "./PValues_MFig_5c.txt", row.names = F, sep = '\t')

pdf("~/path/to/MFig_5c.pdf", height = 6, width = 8)
p1
p2
dev.off()


