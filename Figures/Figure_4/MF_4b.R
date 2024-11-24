library(Seurat)

scATAC_seq_expt_comb_harmony_umap_3 <- readRDS("path/to/seuratobject")
Idents(scATAC_seq_expt_comb_harmony_umap_3) <- scATAC_seq_expt_comb_harmony_umap_3$predicted.id
pdf("~/path/to/pdf", height=10, width=10)
DimPlot(scATAC_seq_expt_comb_harmony_umap_3, cols = c("#22BDC1","#8BB6E1","#292562"), reduction = "umap", label = TRUE, pt.size = .1, label.size = 9) + theme(text = element_text(size=20)) 
dev.off()
