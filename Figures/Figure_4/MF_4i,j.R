library(Seurat)

SeuratObject <- readRDS("~/path/to/SeuratObject")

##Set identity to age. 
Idents(SeuratObject) <- "orig.ident"

## 4i
pdf("~/path/to/pdf")
VlnPlot(object = SeuratObject, features = "Cdkn1a", group.by = 'orig.ident', log=TRUE)+ stat_compare_means(method = "t.test") 
dev.off()

##4j
## Split by subcluster identities. Please confirm the name of metadata column matches up
pdf("~/path/to/pdf")
VlnPlot(object = SeuratObject, features = "Cdkn1a", split.by= 'orig.ident', log=TRUE)+ stat_compare_means(method = "t.test") 
dev.off()
