library(Seurat)

Seurat_Object <- readRDS("~/path/to/object")
genes <- #(c(list of Genes for feature plots))

pdf("~/path/to/pdf")
FeaturePlot(Seurat_Object, features = genes, reduction = "umap", cols = c("lightgrey", "red3"),label = F, ncol = 2, keep.scale="feature", slot = "data")
dev.off()
