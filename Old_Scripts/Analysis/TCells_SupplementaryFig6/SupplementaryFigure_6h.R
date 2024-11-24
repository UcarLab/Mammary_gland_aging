##FIGURE 1F
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)


Counts <- read_xlsx("~/path/to/SupplementaryFigure_6h.xlsx")
Closing <- Counts[,c(1,2)]
Opening <- Counts[,c(1,3)]
colnames(Opening) <- c("CellType","Counts")
colnames(Closing) <- c("CellType","Counts")
Opening$Direction <- "Upregulated"
Closing$Direction <- "Downregulated"
Closing$Counts <- as.numeric(Closing$Counts) * -1
FinalCounts <- rbind(Opening, Closing)
long_DF <- FinalCounts %>% gather("Type", "Count", Counts)
long_DF$CellType <- factor(long_DF$CellType, levels=c("Naive CD4-1", "Naive CD4-2","Memory CD4","Naive CD8-1","Naive CD8-2","CD8 GZMK","CD8 GZMM","γδ T","NK"))
long_DF$Count <- as.numeric(long_DF$Count)

## The γδ T label in the figure sometime does not show up (prints as .. T). We had to manually edit it in final figures. If that happens, please do the same
pdf("~/path/to/SFig_6h.pdf", height = 10, width = 14)
ggplot(long_DF[order(long_DF$Type, decreasing = F),], aes(y=Count, x=CellType, fill = Direction)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("blue","red")) +
  xlab("Gene Counts") +
  theme(axis.text =element_text(size=16),axis.title=element_text(size=16,face="bold")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Cell Types")+
  guides(fill=guide_legend(title="Direction"))
dev.off()
