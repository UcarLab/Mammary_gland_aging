setwd("~/path/to/output")
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape)
library(readxl)
library(purrr)

FinalMotifs <- read.csv("~/path/to/Motifs_SFig_6j.txt", sep = "\t")
Clustering <- FinalMotifs[,c(1,2,4)]
Clustering1 <- cast(Clustering, MotifName ~ CellType, mean, value = 'log_FC_Change')
rownames(Clustering1) <- Clustering1$MotifName
Clustering1 <- Clustering1[,-1]
data <- scale(t(Clustering1))

ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
ord

FinalMotifs$CellType <- factor(FinalMotifs$CellType, levels = c("Memory CD4","Naive CD4-1","Naive CD4-2","γδ T","NK"," Naive CD8-1","Naive CD8-2", "CD8 GZMK", "CD8 GZMM"))
RS <- FinalMotifs[,c(1,5)]
RS1 <- RS %>% group_by_all() %>% summarise(COUNT = n())
RS1 <- RS1[!(RS1$Significance=="" & RS1$COUNT==3),]
RS1 <- as.data.frame(unique(RS1$MotifName))
colnames(RS1) <- "MotifName"

FinalMotifs <- merge(RS1, FinalMotifs)

pdf("~/path/to/SFig_6j.pdf", height = 12, width = 12)
ggplot(FinalMotifs, aes(MotifName, CellType, fill= log_FC_Change)) +
  geom_tile() +
  scale_fill_distiller(palette = "PuOr", ) +
  geom_text(aes(label=Significance), color="black", size=4) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
