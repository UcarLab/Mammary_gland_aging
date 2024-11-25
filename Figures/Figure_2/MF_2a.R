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