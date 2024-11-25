library(Seurat)
library(Signac)
scATAC_seq_expt_comb_harmony_umap_2 <- readRDS("~/path/to/object")
Idents(scATAC_seq_expt_comb_harmony_umap_2) <- scATAC_seq_expt_comb_harmony_umap_2$predicted.id
# Visualization for Figure 2A
pdf("/apth/to/pdf")
DimPlot(scATAC_seq_expt_comb_harmony_umap_2, cols = c("springgreen3","deepskyblue","purple"), reduction = "umap", label = TRUE, pt.size = .1, label.size = 9) + theme(text = element_text(size=20, face = "bold")) 
dev.off()