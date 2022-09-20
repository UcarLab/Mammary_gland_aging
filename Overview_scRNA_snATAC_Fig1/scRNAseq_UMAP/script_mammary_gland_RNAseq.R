#Analyze single cell RNA-seq data for mammary gland aging
#author: Neerja Katiyar
##Aging and Cancer

### Load libraries
library(dplyr)
library(Seurat)
library(cowplot)
library(harmony)
library(pheatmap)

dir.create("./results/Plots")
dir.create("./results/Tables")
dir.create("./results/RDS_files")

scRNA_seq_expt1 <- readRDS("./data/mm_3_18_comb_before_filter_batch1.rds")
scRNA_seq_expt2 <- readRDS("./data/mm_3_18_comb_before_filter_batch2.rds")

scRNA_seq_expt1$expt <- sample(c("expt1"), size = ncol(scRNA_seq_expt1), replace = T)
scRNA_seq_expt2$expt <- sample(c("expt2"), size = ncol(scRNA_seq_expt2), replace = T)

scRNA_seq_expt_comb <- merge(scRNA_seq_expt1, y = c(scRNA_seq_expt2), add.cell.ids=c("expt1", "expt2"), project = "Mouse aging cancer")

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
scRNA_seq_expt_comb[["percent.mt"]] <- PercentageFeatureSet(scRNA_seq_expt_comb, pattern = "^mt-")

# Visualize QC metrics as a violin plot
pdf("./results/Plots/Violin_plot_scRNA_seq_RNA_feature.pdf")
VlnPlot(scRNA_seq_expt_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
saveRDS(scRNA_seq_expt_comb, "./results/RDS_files/scRNA_seq_expt_comb_before_filter_May26_2021.rds")

scRNA_seq_expt_comb <- readRDS("./results/RDS_files/scRNA_seq_expt_comb_before_filter_May26_2021.rds")
pdf("./results/Plots/Plots_qc_orig.pdf", height=10, width=10)
plot1 <- FeatureScatter(scRNA_seq_expt_comb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA_seq_expt_comb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

scRNA_seq_expt_comb <- subset(scRNA_seq_expt_comb, subset = nFeature_RNA > 500 & percent.mt < 10)

#Normalizing the data
scRNA_seq_expt <- NormalizeData(scRNA_seq_expt_comb, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features (feature selection)
scRNA_seq_expt <- FindVariableFeatures(scRNA_seq_expt, selection.method = "vst", nfeatures = length(rownames(scRNA_seq_expt)))

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNA_seq_expt), 10)
top5 <- head(VariableFeatures(scRNA_seq_expt), 5)

# plot variable features with and without labels
pdf("./results/Plots/Plot_variable_feature.pdf", height=10, width=10)
plot1 <- VariableFeaturePlot(scRNA_seq_expt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

saveRDS(scRNA_seq_expt, file="./results/RDS_files/scRNA_seq_expt_orig_May26_2021.rds")

scRNA_seq_expt <- readRDS("./results/RDS_files/scRNA_seq_expt_orig_May26_2021.rds")
all.genes <- rownames(scRNA_seq_expt)
scRNA_seq_expt <- ScaleData(scRNA_seq_expt, features = all.genes)

#Perform linear dimensional reduction
scRNA_seq_expt <- RunPCA(scRNA_seq_expt, features = VariableFeatures(object = scRNA_seq_expt))
print(scRNA_seq_expt[["pca"]], dims = 1:5, nfeatures = 5)

pdf("./results/Plots/Plot_Viz_loading.pdf")
VizDimLoadings(scRNA_seq_expt, dims = 1:2, reduction = "pca")
dev.off()

pdf("./results/Plots/Plot_DimPlot.pdf")
DimPlot(scRNA_seq_expt, reduction = "pca")
dev.off()

pdf("./results/Plots/Heatmap.pdf")
DimHeatmap(scRNA_seq_expt, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

pdf("./results/Plots/Heatmap_dim18.pdf")
DimHeatmap(scRNA_seq_expt, dims = 1:18, cells = 500, balanced = TRUE)
dev.off()

pdf("./results/Plots/Plot_ElbowPlot.pdf")
ElbowPlot(scRNA_seq_expt)
dev.off()

saveRDS(scRNA_seq_expt, file = "./results/RDS_files/scRNA_seq_expt_before_clustering_May26_2021.rds")

#Cluster the cells
scRNA_seq_expt <- FindNeighbors(scRNA_seq_expt, dims = 1:10)
scRNA_seq_expt <- FindClusters(scRNA_seq_expt, resolution = 0.5)

#Run non-linear dimensional reduction (UMAP/tSNE)
scRNA_seq_expt <- RunUMAP(scRNA_seq_expt, dims = 1:10)
pdf("./results/Plots/Plot_UMAP.pdf")
DimPlot(scRNA_seq_expt, reduction = "umap", label=TRUE)
dev.off()

saveRDS(scRNA_seq_expt, file = "./results/RDS_files/scRNA_seq_expt_out_May26_2021.rds")

scRNA_seq_expt <- readRDS(".results//RDS_files/scRNA_seq_expt_out_May26_2021.rds")
# Visualization
p1 <- DimPlot(scRNA_seq_expt, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(scRNA_seq_expt, reduction = "umap", label = TRUE)

pdf("./results/Plots/Dimplot_by_time_and_orig_cluster.pdf", height=10, width=20)
plot_grid(p1,p2)
dev.off()

# Visualization
p3 <- DimPlot(scRNA_seq_expt, reduction = "umap", group.by = "expt", label = TRUE)
p4 <- DimPlot(scRNA_seq_expt, reduction = "umap", label = TRUE)

pdf("./results/Plots/Dimplot_by_expt_and_orig_cluster.pdf", height=10, width=20)
plot_grid(p3,p4)
dev.off()

scRNA_seq_expt <- scRNA_seq_expt %>% RunHarmony("expt", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(scRNA_seq_expt, 'harmony')
harmony_embeddings[1:5, 1:5]
options(repr.plot.height = 5, repr.plot.width = 12)
pdf("./results/Plots/DimPlot_after_harmony.pdf")
p1 <- DimPlot(object = scRNA_seq_expt, reduction = "harmony", pt.size = .1, group.by = "expt")
p2 <- VlnPlot(object = scRNA_seq_expt, features = "harmony_1", group.by = "expt", pt.size = .1)
plot_grid(p1,p2)
dev.off()

scRNA_seq_expt <- scRNA_seq_expt %>% 
    RunUMAP(reduction = "harmony", dims = 1:10) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

options(repr.plot.height = 4, repr.plot.width = 6)
pdf("./results/Plots/Dimplot_umap_after_harmony.pdf")
DimPlot(scRNA_seq_expt, reduction = "umap", label = TRUE, pt.size = .1, label.size=5)
dev.off()

# Visualization
p3 <- DimPlot(scRNA_seq_expt, reduction = "umap", group.by = "expt", label = TRUE, label.size=5)
p4 <- DimPlot(scRNA_seq_expt, reduction = "umap", label = TRUE, label.size=5)

pdf("./results/Plots/Dimplot_by_expt_and_orig_cluster_after_harmony.pdf", height=10, width=20)
plot_grid(p3,p4)
dev.off()

p5 <- DimPlot(scRNA_seq_expt, reduction = "umap", split.by = "expt", label = TRUE, label.size=5)
p6 <- DimPlot(scRNA_seq_expt, reduction = "umap", label = TRUE, label.size=5)

pdf("./results/Plots/Dimplot_split_expt_and_orig_cluster_after_harmony.pdf", height=10, width=20)
plot_grid(p5,p6)
dev.off()

p7 <- DimPlot(scRNA_seq_expt, reduction = "umap", split.by = "orig.ident", label = TRUE, label.size=5)
p8 <- DimPlot(scRNA_seq_expt, reduction = "umap", label = TRUE, label.size=5)

pdf("./results/Plots/Dimplot_split_age_and_orig_cluster_after_harmony.pdf", height=10, width=20)
plot_grid(p7,p8)
dev.off()

saveRDS(scRNA_seq_expt, file = "./results/RDS_files/seurat_black6_scRNA_May26_2021.rds")

################Find all markers #####################
# find markers for every cluster compared to all remaining cells, report only the positive ones
scRNA_seq_expt.markers <- FindAllMarkers(scRNA_seq_expt, only.pos = FALSE, min.pct = 0.25, test.use = 'LR', latent.vars = 'expt', logfc.threshold = 0.25)
scRNA_seq_expt.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.table(scRNA_seq_expt.markers, "./results/Tables/Table_markers_all_prelim_pos_neg.txt", sep="\t", quote=FALSE)

top10 <- scRNA_seq_expt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, "./results/Tables/Table_markers_top10.txt", sep="\t", quote=FALSE)

top5 <- scRNA_seq_expt.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.table(top5, "./results/Tables/Table_markers_top5.txt", sep="\t", quote=FALSE)

pdf("./results/Plots/Heatmap_top_expressed.pdf", height=20, width=30)
DoHeatmap(scRNA_seq_expt, features = top10$gene) + NoLegend()
dev.off()

table_expt_rep <- table(scRNA_seq_expt$expt, scRNA_seq_expt$replicate)
write.table(table_expt_rep, "./results/Tables/Table_expt_rep.txt", sep="\t", quote=F)

table_expt_age <- table(scRNA_seq_expt$expt, scRNA_seq_expt$orig.ident)
write.table(table_expt_age, "./results/Tables/Table_expt_age.txt", sep="\t", quote=F)

saveRDS(scRNA_seq_expt, file = "./results/RDS_files/seurat_black6_scRNA_all_markers.rds")

#################################################
scRNA_seq_expt <- readRDS("../RDS_files/seurat_black6_scRNA_May26_2021.rds")

pdf("./results/Plots/FeaturePlot_Immune_Ptprc.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ptprc"))
dev.off()

pdf("./results/Plots/FeaturePlot_Epithelial_Epcam.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Epcam"))
dev.off()

pdf("./results/Plots/FeaturePlot_Dendritic_Itgax.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Itgax"))
dev.off()

pdf("./results/Plots/FeaturePlot_Bcells_Blnk.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Blnk"))
dev.off()

pdf("./results/Plots/FeaturePlots_Bcells_Cd79b.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd79b"))
dev.off()

pdf("./results/Plots/FeaturePlot_Tcells_naive_S100a4.pdf")
FeaturePlot(scRNA_seq_expt, features = c("S100a4"))
dev.off()

pdf("./results/Plots/FeaturePlot_Tcells_naive_Cd44.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd44"))
dev.off()

pdf("./results/Plots/FeaturePlot_Tcells_naive_Sell.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Sell"))
dev.off()

pdf("./results/Plots/FeaturePlot_Tcells_naive_Il7r.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Il7r"))
dev.off()

pdf("./results/Plots/FeaturePlot_Tcells_Cd3d.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd3d"))
dev.off()

pdf("./results/Plots/FeaturePlot_CD4_CD8_pos_Cd4.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd4"))
dev.off()

pdf("./results/Plots/FeaturePlot_CD4_CD8_pos_Cd8a.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd8a"))
dev.off()

pdf("./results/Plots/FeaturePlot_CD4_CD8_pos_Cd8b1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd8b1"))
dev.off()

pdf("./results/Plots/FeaturePlot_plasma_Bcells_marker_Jchain.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Jchain"))
dev.off()

pdf("./results/Plots/FeaturePlot_S100A4.pdf")
FeaturePlot(scRNA_seq_expt, features = c("S100a4"))
dev.off()

pdf("./results/Plots/FeaturePlot_Cxcr5.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cxcr5"))

pdf("./results/Plots/FeaturePlot_naive_Sell.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Sell"))
dev.off()

pdf("./results/Plots/FeaturePlot_naive_Lef1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Lef1"))
dev.off()

pdf("./results/Plots/FeaturePlot_naive_Ccr7.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ccr7"))
dev.off()

pdf("./results/Plots/FeaturePlot_some_immune_cells_Pdgfra.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Pdgfra"))
dev.off()

pdf("./results/Plots/FeaturePlot_luminal_Krt18.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Krt18"))
dev.off()

pdf("./results/Plots/FeaturePlot_luminal-HS_Cited1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cited1"))
dev.off()

pdf("./results/Plots/FeaturePlot_luminal-AV_Mfge8.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Mfge8"))
dev.off()

pdf("./results/Plots/FeaturePlot_Myoepithelial_Krt14.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Krt14"))
dev.off()

pdf("./results/Plots/FeaturePlot_Myoepithelial_Acta2.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Acta2"))
dev.off()

pdf("./results/Plots/FeaturePlot_Fibroblasts-ECM_producing_Fn1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Fn1"))
dev.off()

pdf("./results/Plots/FeaturePlot_Fibroblasts-ECM_producing_Col1a1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Col1a1"))
dev.off()

pdf("./results/Plots/FeaturePlot_Fibroblasts-Contractile_Acta2.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Acta2"))
dev.off()

pdf("./results/Plots/FeaturePlot_Fibroblasts-Contractile_Myl9.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Myl9"))
dev.off()

pdf("./results/Plots/FeaturePlot_Vascular_Pecam1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Pecam1"))
dev.off()

pdf("./results/Plots/FeaturePlot_Pericytes_Des.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Des"))
dev.off()

pdf("./results/Plots/FeaturePlot_Immune.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ptprc"))
dev.off()

pdf("./results/Plots/FeaturePlot_Macrophages_Aif1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Aif1"))
dev.off()

pdf("./results/Plots/FeaturePlot_Macrophages_Itgax.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Itgax"))
dev.off()

pdf("./results/Plots/FeaturePlot_Dendritic_Cd209a.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd209a"))
dev.off()

pdf("./results/Plots/FeaturePlot_NaturalKiller_Nkg7.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Nkg7"))
dev.off()

pdf("./results/Plots/FeaturePlot_Cytotoxic_Gzmb.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Gzmb"))
dev.off()

#Rename_idents, perform QC.
new.cluster.ids = c(`0`="Tcells_naive",`1`="Tcells_mem", `2`="Bcells", `3`="Bcells", `4`="Tcells_naive", `5`="Fibroblasts", `6`="Tcells_naive", `7`="Myoepithelial", `8`="Luminal-AV", `9`="Dendritic/Macrophages",`10`="Luminal-HS", `11`="Dendritic/Macrophages", `12`="Vascular", `13` = "Bcells", `14`="Jchain", `15`="Tcells_mem", `16`="Pericytes", `17`="Plasma")

scRNA_seq_expt <- RenameIdents(scRNA_seq_expt, new.cluster.ids)
pdf("Annotated_clusters_UMAP.pdf")
DimPlot(scRNA_seq_expt, cols = c("#FFCC33", "#FF9933", "#746CB1", "#1BBDC1", "#09713A", "#3DB54A", "#c00000", "#A0D082", "#8AB6E1", "#ff99cc", "#262262", "#BD77B2"), reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

scRNA_seq_expt_new <- subset(scRNA_seq_expt, idents = c("Doublet"), invert = TRUE)
pdf("Annotated_clusters_doublet_rmvd_UMAP_new.pdf")
DimPlot(scRNA_seq_expt_new, cols = c("#FFCC33", "#FF9933", "#746CB1", "#1BBDC1", "#09713A", "#3DB54A", "#c00000", "#A0D082", "#8AB6E1", "#262262", "#BD77B2"), reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

Idents(scRNA_seq_expt_new) <- scRNA_seq_expt@meta.data$orig.ident

scRNA_seq_expt_new$orig.ident <- factor(x = scRNA_seq_expt_new$orig.ident, levels = c("mm10_3mths", "mm10_18mths"))
p9 <- DimPlot(scRNA_seq_expt_new, cols = c("blue", "blue"), reduction = "umap", split.by = "orig.ident", label = TRUE)
p10 <- DimPlot(scRNA_seq_expt_new, cols = c("blue", "blue"), reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size=2)

pdf("Dimplot_split_age_and_orig_cluster_after_harmony_new2.pdf", height=10, width=20)
plot_grid(p9)
dev.off()

pdf("Dimplot_split_age_and_orig_cluster_after_harmony_ptsize2.pdf", height=10, width=20)
plot_grid(p10)
dev.off()

# Save a single object to a file
saveRDS(scRNA_seq_expt, "seuratObject_rename_ident_May26_2021.rds")

#Find_markers_avg_exprn
mydata_new = readRDS("seuratObject_rename_ident.rds")
mydata_new_rmv_doublet <- subset(mydata_new, idents = c("Doublet"), invert = TRUE)
scRNA_seq_expt.markers_new <- FindAllMarkers(mydata_new_rmv_doublet, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
scRNA_seq_expt.markers_new %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.table(scRNA_seq_expt.markers_new, "Table_markers_all_pos_neg_rename_clusters.txt", sep="\t", quote=FALSE)
top10_new <- scRNA_seq_expt.markers_new %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

pdf("Heatmap_top_expressed_rename_clusters_new.pdf", height=15, width=20)
DoHeatmap(mydata_new, features = top10_new$gene, size=3)
dev.off()

Clust_Celltype_clusters_orig <- subset(cluster.averages_tp_new_orig$RNA@scale.data, rownames(cluster.averages_tp_new_orig$RNA@scale.data) %in% marker_list)

cluster.averages_tp_new <- AverageExpression(mydata_new, return.seurat = TRUE)

Clust_Celltype_all <- subset(cluster.averages_tp_new$RNA@scale.data, rownames(cluster.averages_tp_new$RNA@scale.data) %in% top10_new$gene)

pdf("Heatmap_all_markers_cluster.pdf")
pheatmap(Clust_Celltype_all, annotation_legend = TRUE, cluster_cols=FALSE, cellwidth = 20,treeheight_row = 0, treeheight_col = 0, fontsize_row = 4, fontsize_col = 8)
dev.off()

