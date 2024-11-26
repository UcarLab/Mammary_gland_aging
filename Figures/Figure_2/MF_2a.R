library(Seurat)
scRNA_seq_expt_comb_harmony_umap_2 <- readRDS("~/path/to/object")
# Visualization for Figure 2A
new.cluster.ids <- c("Myoepithelial", "Luminal-AV", "Myoepithelial", "Luminal-AV", "Luminal-HS", "Luminal-AV", "Luminal-HS", "Myoepithelial", "Luminal-AV", "Luminal-HS", "Myoepithelial")
seurat.cluster.ids <- scRNA_seq_expt_comb_harmony_umap_2$seurat_clusters
names(new.cluster.ids) <- levels(scRNA_seq_expt_comb_harmony_umap_2)
scRNA_seq_expt_comb_harmony_umap_2 <- RenameIdents(scRNA_seq_expt_comb_harmony_umap_2,new.cluster.ids)
pdf("/apth/to/pdf")
DimPlot(scRNA_seq_expt_comb_harmony_umap_2, cols = c("springgreen3","deepskyblue","purple"), reduction = "umap", label = TRUE, pt.size = .1, label.size = 9) + theme(text = element_text(size=20, face = "bold")) 
dev.off()

## For bar plot
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
long_DF$CellType <- factor(long_DF$CellType, levels=c("Luminal-AV", "Luminal-HS", "Myoepithelial"))

pdf("path/to/pdf", height = 10, width = 14)
ggplot(long_DF[order(long_DF$Type, decreasing = F),], aes(y=Count, x=CellType, fill = Direction)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("springgreen3","deepskyblue","purple")) + 
  xlab("Gene Counts") + 
  ylab("Cell Types")+
  guides(fill=guide_legend(title="Direction"))
dev.off()
