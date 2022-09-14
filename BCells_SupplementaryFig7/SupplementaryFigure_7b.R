setwd("~/path/to/output")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)

BCells_RNA <- readRDS("./BCells_RNA.rds")
Idents(BCells_RNA) <- "SubclusterAnnotations"

my_levels <- c("C1","C2","C3", "C4", "C5", "C6", "C7")
BCells_RNA@meta.data$SubclusterAnnotations <- factor(x = BCells_RNA@meta.data$SubclusterAnnotations, levels = my_levels)
colr <- c("orangered1","firebrick4","purple","lightblue2","dodgerblue3","goldenrod2","darkblue")
p1 <- DimPlot(BCells_RNA, reduction = "humap7", group.by = "SubclusterAnnotations", label = F, label.size = 8, cols = colr) 

Meta <- BCells_RNA@meta.data
Meta <- Meta[,c(1,4,12)]
Meta <- Meta %>% group_by_all() %>% summarise(COUNT = n())
Meta$replicate <- gsub("M3","M18",Meta$replicate)
NewMeta = cast(Meta, replicate + SubclusterAnnotations ~ orig.ident, sum) 
NewMeta$Ratio <- log2(NewMeta$mm10_18mths/NewMeta$mm10_3mths)
NewMeta <- NewMeta[,c(1,2,5)]
NewMeta$Annotation <- as.character(NewMeta$SubclusterAnnotations)
NewMeta$Ratio <- as.numeric(NewMeta$Ratio)

p2 <- ggplot(NewMeta,aes(x = reorder(SubclusterAnnotations, -Ratio), y = Ratio, fill = SubclusterAnnotations)) + 
  geom_boxplot()+
  scale_fill_manual(values=(colr)) +
  xlab("CellTypes") +
  ylab("log2 FC (Old/Young)") +
  ggtitle("Cell Type Proportions") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  ylim(-4,2)

## For Calculating P Values, we added a zero base ratio to refer to. Then we used the compare_means function from ggpubr to calculate the p values. We added the p values to the plots manually later on.
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep1","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep2","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep3","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep4","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep5","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep6","Base",0)
NewMeta$Annotation <- as.factor(NewMeta$SubclusterAnnotations)
stat <- compare_means(Ratio ~ SubclusterAnnotations,  data = NewMeta, ref.group = "Base", method = "t.test")
write.table(stat, "./PValues_SFig_7b.txt", row.names = F, sep = '\t')

pdf("~/path/to/SFig_7b.pdf", height = 6, width = 8)
p1
p2
dev.off()

