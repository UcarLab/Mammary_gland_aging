library(readxl)  
library(plyr)
library(tibble)
library(ggplot2)
library(tidyr)
##Needs an excel file with each cell type as rows and upregulated and downregulated differential gene/peak counts as individual columns
Counts <- read_xlsx("path/to/diffgenecounts")
Counts <- as.data.frame(Counts)

## Assuming upregulated counts are in 2nd column and downregulated counts are in column 3.
Opening <- Counts[,c(1,2)]
Closing <- Counts[,c(1,3)]
colnames(Opening) <- c("CellType","Counts")
colnames(Closing) <- c("CellType","Counts")
Opening$Direction <- "Upregulated"
Closing$Direction <- "Downregulated"
Closing$Counts <- (Closing$Counts) * -1

FinalCounts <- rbind(Opening, Closing)
long_DF <- FinalCounts %>% gather("Type", "Count", Counts)
long_DF$CellType <- factor(long_DF$CellType, levels=c("Epi-C1","Epi-C2","Epi-C3","Epi-C4","Epi-C5","Epi-C6","Epi-C7","Epi-C8","Epi-C9","Epi-C10","Epi-C11",))

pdf("path/to/pdf", height = 10, width = 14)
ggplot(long_DF[order(long_DF$Type, decreasing = F),], aes(y=Count, x=CellType, fill = Direction)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("deepskyblue","deepskyblue","deepskyblue","springgreen3","springgreen3","springgreen3","springgreen3","purple","purple","purple","purple")) + 
  xlab("Gene Counts") + 
  ylab("Cell Types")+
  guides(fill=guide_legend(title="Direction"))
dev.off()
