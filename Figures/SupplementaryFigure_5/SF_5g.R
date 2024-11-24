library(Seurat)
library(ggplot)
library(dplyr)
scRNA_seq_expt_comb_harmony_umap_3 <- readRDS("path/to/object")
pdf("path/to/pdf", height=20, width=20)
DotPlot(scRNA_seq_expt_comb_harmony_umap_3, features = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Pdgfra", 
                                                         "Acta2", "Myl9", "Mylk", "Myh11",
                                                         "Rgs5", "Des", "Notch3", "Ndufa4l2", "Higd1b",
                                                         "Pax7", "Bmp4",
                                                         "Pecam1", "Cdh5", "Eng", "Sox17", "Sele",
                                                         "Prox1", "Lyve1", "Itga2b"),
        dot.scale = 5) + coord_flip() + theme(text = element_text(size=12))
dev.off()
