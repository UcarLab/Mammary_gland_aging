###Plotting Single Cell Heatmap
library(Seurat)
library(ggplot)
library(dplyr)
scRNA_seq_expt_comb_harmony_umap_3 <- readRDS("path/to/object")
my_levels <- c("0", "6", "1", "4", "3", "7",  "5", "10", "2", "8", "9")
scRNA_seq_expt_comb_harmony_umap_3@active.ident <- factor(x = scRNA_seq_expt_comb_harmony_umap_3@active.ident, levels = my_levels)
scRNA_seq_expt_comb_markers_3 <- FindAllMarkers(scRNA_seq_expt_comb_harmony_umap_3, only.pos = FALSE, min.pct = 0.25, test.use = 'LR', latent.vars = 'expt', logfc.threshold = 0.25)
top5 <- scRNA_seq_expt_comb_markers_3 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

pdf("path/to/pdf", height=9, width=7)
DoHeatmap(scRNA_seq_expt_comb_harmony_umap_3, features = top5$gene, raster = F, lines.width =5, 
          group.colors = c("springgreen4", "deepskyblue", "goldenrod2", "purple", "turquoise", "lightblue2", "firebrick4", 
                           "violetred", "darkblue", "lightpink", "darkorange"))+ theme(text = element_text(size=5))
dev.off()
