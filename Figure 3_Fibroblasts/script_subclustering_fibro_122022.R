##Subclustering Fibroblasts
##04112022
##Based off Neerja's Clustering script 
##This is based off of the Seurat Tutorial - https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
##Harmony Vignette: https://portals.broadinstitute.org/harmony/SeuratV3.html

##Revisting and adding pericytes/vascular

##Loading the libraries
library(dplyr)
library(Seurat)
library(cowplot)
library(harmony)
library(ggplot2)
library(reshape2)
library(ggpubr)

###Uploading the RDS with named clusters
BA_obj <- readRDS("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/seuratObject_rename_ident_May26_2021.rds")

###Survey the metadata of the uploaded file
dim(BA_obj@meta.data)   #48180 cells
table(Idents(BA_obj))   #3130 Fibroblasts
                        #417 Pericytes
                        #686 Vascular

###Generating column that contains the idents so that it is easier to subset
BA_obj$celltype <- paste(Idents(BA_obj))
head(BA_obj@meta.data)

###Looking at the data by cell type - verify this is the most recent clustering
DimPlot(BA_obj, 
        group.by = c("celltype"), label = T)

###Creating metadata object and verifying it looks as expected
Meta <- BA_obj@meta.data
dim(Meta)
head(Meta)

###Subsetting the fibroblasts and generating a list of the cell identifiers
Meta %>% dplyr::filter(Idents(BA_obj) == 'Fibroblasts'|Idents(BA_obj) == 'Pericytes'|Idents(BA_obj) == 'Vascular') -> Stroma
dim(Stroma) #4233 cells
cell_to_keep <- row.names(Stroma)
length(cell_to_keep)
#write.csv(cell_to_keep, "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/stroma_subset_04152022.csv")
rm(BA_obj)

###Loading the unfiltered/unlabeled RDS file
scRNA_seq_expt_comb_before_filter <- readRDS("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_before_filter_May26_2021.rds")

###Subsetting the unfiltered RDS file for the fibroblasts
subset_stroma_combined = subset(scRNA_seq_expt_comb_before_filter, cells = cell_to_keep)
dim(subset_stroma_combined@meta.data) #4233 cells
head(subset_stroma_combined@meta.data)
rm(scRNA_seq_expt_comb_before_filter)

#saveRDS(subset_stroma_combined, file = "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/stroma_subset_04152022_BA.rds")
#subset_stroma_combined <- readRDS("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/stroma_subset_04152022_BA.rds")

##Re-run QC just in case
scRNA_seq_expt_comb_filtered <- subset(subset_stroma_combined, subset = nFeature_RNA > 500 & percent.mt < 10)
rm(subset_stroma_combined)

#Normalizing the data
# "LogNormalize" normalizes the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000), and log-transforms the result
# Normalized values are stored in scRNA_seq_expt_comb_norm[["RNA"]]@data.
scRNA_seq_expt_comb_norm <- NormalizeData(scRNA_seq_expt_comb_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
rm(scRNA_seq_expt_comb_filtered)

#Identification of highly variable features (feature selection)
#Calculates features that exhibit high cell-to-cell variation
#We include all genes as features
scRNA_seq_expt_comb_hvf <- FindVariableFeatures(scRNA_seq_expt_comb_norm, selection.method = "vst", nfeatures = length(rownames(scRNA_seq_expt_comb_norm)))
rm(scRNA_seq_expt_comb_norm)

#Scaling the data
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data
all.genes <- rownames(scRNA_seq_expt_comb_hvf)
scRNA_seq_expt_comb_scal <- ScaleData(scRNA_seq_expt_comb_hvf, features = all.genes)
rm(scRNA_seq_expt_comb_hvf)

#Perform linear dimensional reduction
#Mostly using the elbow plot to determine number of principle components to use
#Seurat advises users to err on the higher side when choosing this parameter. 
scRNA_seq_expt_comb_pca <- RunPCA(scRNA_seq_expt_comb_scal, features = VariableFeatures(object = scRNA_seq_expt_comb_scal))
DimHeatmap(scRNA_seq_expt_comb_pca, dims = 1:24, cells = 500, balanced = TRUE)
ElbowPlot(scRNA_seq_expt_comb_pca, ndims = 100)
rm(scRNA_seq_expt_comb_scal)

#Cluster the cells
#This step constructs a KNN graph based on the euclidean distance in PCA space, 
#and refines the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity).
#resolution parameter that sets the 'granularity' of the downstream clustering
scRNA_seq_expt_comb_neighbors <- FindNeighbors(scRNA_seq_expt_comb_pca, dims = 1:18)
rm(scRNA_seq_expt_comb_pca)
scRNA_seq_expt_comb_clusters <- FindClusters(scRNA_seq_expt_comb_neighbors, resolution = 0.8)
rm(scRNA_seq_expt_comb_neighbors)

#Run non-linear dimensional reduction (UMAP/tSNE)
scRNA_seq_expt_comb_umap <- RunUMAP(scRNA_seq_expt_comb_clusters, dims = 1:20)
DimPlot(scRNA_seq_expt_comb_umap, reduction = "umap", label=TRUE)
rm(scRNA_seq_expt_comb_clusters)

#Passing Seurat object to Harmony and integrating out experiment
#Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.RNA.harmony; see ?make.names for more details on syntax validity
#Took 6/10 iterations
scRNA_seq_expt_comb_harmony <- scRNA_seq_expt_comb_umap %>% RunHarmony("expt", plot_convergence = TRUE)
rm(scRNA_seq_expt_comb_umap)

#To access embeddings
harmony_embeddings <- Embeddings(scRNA_seq_expt_comb_harmony, 'harmony')

scRNA_seq_expt_comb_harmony_umap <- scRNA_seq_expt_comb_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:18) %>%
  FindNeighbors(reduction = "harmony", dims = 1:18) %>%
  FindClusters(resolution = 0.8) %>%
  identity()

rm(scRNA_seq_expt_comb_harmony)

#Reordering idents to put young first in visualization
#Found here:https://github.com/satijalab/seurat/issues/2471
scRNA_seq_expt_comb_harmony_umap$orig.ident <- factor(x = scRNA_seq_expt_comb_harmony_umap$orig.ident, levels = c("mm10_3mths", "mm10_18mths"))

# Visualization
DimPlot(scRNA_seq_expt_comb_harmony_umap, reduction = "umap", label = TRUE, pt.size = .1, label.size = 9) + theme(text = element_text(size=20, face = "bold")) 
DimPlot(scRNA_seq_expt_comb_harmony_umap, reduction = "umap", group.by = "expt", label = TRUE)
DimPlot(scRNA_seq_expt_comb_harmony_umap, reduction = "umap", split.by = "replicate", label = TRUE)
DimPlot(scRNA_seq_expt_comb_harmony_umap, reduction = "umap", split.by = "orig.ident", label = TRUE, label.size = 9) + theme(text = element_text(size=20, face = "bold")) 

#Save the rds file
#saveRDS(scRNA_seq_expt_comb_harmony_umap, file = "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_harmony_umap_initial_stroma_subset_04152022_BA.rds")
scRNA_seq_expt_comb_harmony_umap <- readRDS("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_harmony_umap_initial_stroma_subset_04152022_BA.rds")

##Find markers that define clusters
scRNA_seq_expt_comb_markers <- FindAllMarkers(scRNA_seq_expt_comb_harmony_umap, only.pos = FALSE, min.pct = 0.25, test.use = 'LR', latent.vars = 'expt', logfc.threshold = 0.25)
#write.csv(scRNA_seq_expt_comb_markers, "C:/Users/angarb/Box/Aging single cells/Subclustering/Single_Cell_June2021/Table_markers_all_prelim_pos_neg.csv")

top10 <- scRNA_seq_expt_comb_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#write.csv(top10, "Fibro_050322_Table_markers_top10.txt", sep="\t", quote=FALSE)

##Identifying immune and lymphatic clusters
FeaturePlot(scRNA_seq_expt_comb_harmony_umap, features = c("Cd3d", "Ptprc", "Ms4a1"))
FeaturePlot(scRNA_seq_expt_comb_harmony_umap, features = c("Pecam1", "Pdpn", "Prox1", "Lyve1", "Itga2b"))
#, "Sox17", "Mmrn1", "Prox1", "Flt4", "Ccl21a"
#https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000704

##########################################
###Removing muscle contamination and immune doublets
###Creating metadata object and verifying it looks as expected
dim(scRNA_seq_expt_comb_harmony_umap@meta.data)   #4233 cells
table(Idents(scRNA_seq_expt_comb_harmony_umap))

Meta_Second <- scRNA_seq_expt_comb_harmony_umap@meta.data
Meta_Second %>% dplyr::filter(Idents(scRNA_seq_expt_comb_harmony_umap) == "16"| Idents(scRNA_seq_expt_comb_harmony_umap) == "10") -> Doublets_Muscle
dim(Doublets_Muscle) #154 cells
cell_to_remove <- row.names(Doublets_Muscle)
length(cell_to_remove)

final_cells_to_keep <- setdiff(cell_to_keep, cell_to_remove)
length(final_cells_to_keep) #4079 cells

subset_fibro_combined_2 = subset(scRNA_seq_expt_comb_before_filter, cells = final_cells_to_keep)
table(subset_fibro_combined_2$orig.ident)

#Normalizing the data
scRNA_seq_expt_comb_norm_2 <- NormalizeData(subset_fibro_combined_2, normalization.method = "LogNormalize", scale.factor = 10000)
rm(subset_fibro_combined_2)

#Identification of highly variable features (feature selection)
scRNA_seq_expt_comb_hvf_2 <- FindVariableFeatures(scRNA_seq_expt_comb_norm_2, selection.method = "vst", nfeatures = length(rownames(scRNA_seq_expt_comb_norm_2)))
rm(scRNA_seq_expt_comb_norm_2)

all.genes.2 <- rownames(scRNA_seq_expt_comb_hvf_2)
scRNA_seq_expt_comb_scal_2 <- ScaleData(scRNA_seq_expt_comb_hvf_2, features = all.genes.2)
rm(scRNA_seq_expt_comb_hvf_2)

#Perform linear dimensional reduction
scRNA_seq_expt_comb_pca_2 <- RunPCA(scRNA_seq_expt_comb_scal_2, features = VariableFeatures(object = scRNA_seq_expt_comb_scal_2))
rm(scRNA_seq_expt_comb_scal_2)

ElbowPlot(scRNA_seq_expt_comb_pca_2)

#Cluster the cells
scRNA_seq_expt_comb_neighbors_2 <- FindNeighbors(scRNA_seq_expt_comb_pca_2, dims = 1:14)
rm(scRNA_seq_expt_comb_pca_2)

scRNA_seq_expt_comb_clusters_2 <- FindClusters(scRNA_seq_expt_comb_neighbors_2, resolution = 0.5)
rm(scRNA_seq_expt_comb_neighbors_2)

#Run non-linear dimensional reduction (UMAP/tSNE)
scRNA_seq_expt_comb_umap_2 <- RunUMAP(scRNA_seq_expt_comb_clusters_2, dims = 1:14)
rm(scRNA_seq_expt_comb_clusters_2)

DimPlot(scRNA_seq_expt_comb_umap_2, reduction = "umap", label=TRUE)

scRNA_seq_expt_comb_harmony_2 <- scRNA_seq_expt_comb_umap_2 %>% RunHarmony("expt", plot_convergence = TRUE)
rm(scRNA_seq_expt_comb_umap_2)

harmony_embeddings_2 <- Embeddings(scRNA_seq_expt_comb_harmony_2, 'harmony')

scRNA_seq_expt_comb_harmony_umap_2 <- scRNA_seq_expt_comb_harmony_2 %>%
  RunUMAP(reduction = "harmony", dims = 1:14) %>%
  FindNeighbors(reduction = "harmony", dims = 1:14) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

DimPlot(scRNA_seq_expt_comb_harmony_umap_2, reduction = "umap", label = TRUE, pt.size = .1, label.size = 9) + theme(text = element_text(size=20, face = "bold")) 
DimPlot(scRNA_seq_expt_comb_harmony_umap_2, reduction = "umap", split.by = "orig.ident", label = TRUE, label.size = 9) + theme(text = element_text(size=20, face = "bold")) 

##Identifying immune and lymphatic clusters
FeaturePlot(scRNA_seq_expt_comb_harmony_umap_2, features = c("Cd3d", "Ptprc", "Ms4a1"))
FeaturePlot(scRNA_seq_expt_comb_harmony_umap_2, features = c("Pecam1", "Pdpn", "Prox1", "Lyve1", "Itga2b"))

VlnPlot(scRNA_seq_expt_comb_harmony_umap_2, features = "nFeature_RNA") 

##Find markers that define clusters
scRNA_seq_expt_comb_markers <- FindAllMarkers(scRNA_seq_expt_comb_harmony_umap_2, only.pos = FALSE, min.pct = 0.25, test.use = 'LR', latent.vars = 'expt', logfc.threshold = 0.25)
#write.csv(scRNA_seq_expt_comb_markers, "C:/Users/angarb/Box/Aging single cells/Subclustering/Single_Cell_June2021/Table_markers_all_prelim_pos_neg.csv")

top10 <- scRNA_seq_expt_comb_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#write.csv(top10, "Fibro_Round2_050622_Table_markers_top10.txt", sep="\t", quote=FALSE)

DoHeatmap(scRNA_seq_expt_comb_harmony_umap_2, features = top10$gene) + NoLegend() + theme(text = element_text(size=10)) 

#Save the rds file
#saveRDS(scRNA_seq_expt_comb_harmony_umap_2, file = "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_harmony_umap_second_stroma_subset_05092022_BA.rds")
scRNA_seq_expt_comb_harmony_umap_2 <- readRDS("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_harmony_umap_second_stroma_subset_05092022_BA.rds")

scRNA_seq_expt_comb_filtered_2 <- subset(scRNA_seq_expt_comb_harmony_umap_2, subset = nFeature_RNA > 900 & percent.mt < 10)
dim(scRNA_seq_expt_comb_filtered_2@meta.data)#3905

##########################################
###Removing low quality cells
###Creating metadata object and verifying it looks as expected
dim(scRNA_seq_expt_comb_harmony_umap_2@meta.data)   #4079 cells
table(Idents(scRNA_seq_expt_comb_harmony_umap_2))

Meta_Third <- scRNA_seq_expt_comb_harmony_umap_2@meta.data
Meta_Third %>% dplyr::filter(Idents(scRNA_seq_expt_comb_harmony_umap_2) == "6"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "10") -> Low_Quality
dim(Low_Quality) #207 cells
cell_to_remove <- row.names(Low_Quality)
length(cell_to_remove)

df <- subset(scRNA_seq_expt_comb_harmony_umap_2, subset = Cd3d>0 & Ptprc>0 | Ms4a1>0 & Ptprc>0)
Meta_Fourth <- df@meta.data
dim(Meta_Fourth) #41 cells
cell_to_remove.2 <- row.names(Meta_Fourth)
length(cell_to_remove.2)

final_cells_to_keep_3 <- setdiff(final_cells_to_keep, cell_to_remove)
final_cells_to_keep_4 <- setdiff(final_cells_to_keep_3, cell_to_remove.2)

length(final_cells_to_keep_4) #3832 cells

subset_fibro_combined_3 = subset(scRNA_seq_expt_comb_before_filter, cells = final_cells_to_keep_4)
table(subset_fibro_combined_3$orig.ident)

#Normalizing the data
scRNA_seq_expt_comb_norm_3 <- NormalizeData(subset_fibro_combined_3, normalization.method = "LogNormalize", scale.factor = 10000)
rm(subset_fibro_combined_3)

#Identification of highly variable features (feature selection)
scRNA_seq_expt_comb_hvf_3 <- FindVariableFeatures(scRNA_seq_expt_comb_norm_3, selection.method = "vst", nfeatures = length(rownames(scRNA_seq_expt_comb_norm_3)))
rm(scRNA_seq_expt_comb_norm_3)

all.genes.3 <- rownames(scRNA_seq_expt_comb_hvf_3)
scRNA_seq_expt_comb_scal_3 <- ScaleData(scRNA_seq_expt_comb_hvf_3, features = all.genes.3)
rm(scRNA_seq_expt_comb_hvf_3)

#Perform linear dimensional reduction
scRNA_seq_expt_comb_pca_3 <- RunPCA(scRNA_seq_expt_comb_scal_3, features = VariableFeatures(object = scRNA_seq_expt_comb_scal_3))
rm(scRNA_seq_expt_comb_scal_3)

ElbowPlot(scRNA_seq_expt_comb_pca_3)

#Cluster the cells
scRNA_seq_expt_comb_neighbors_3 <- FindNeighbors(scRNA_seq_expt_comb_pca_3, dims = 1:14)
rm(scRNA_seq_expt_comb_pca_3)

scRNA_seq_expt_comb_clusters_3 <- FindClusters(scRNA_seq_expt_comb_neighbors_3, resolution = 0.4)
rm(scRNA_seq_expt_comb_neighbors_3)

#Run non-linear dimensional reduction (UMAP/tSNE)
scRNA_seq_expt_comb_umap_3 <- RunUMAP(scRNA_seq_expt_comb_clusters_3, dims = 1:14)
rm(scRNA_seq_expt_comb_clusters_3)

DimPlot(scRNA_seq_expt_comb_umap_3, reduction = "umap", label=TRUE)

scRNA_seq_expt_comb_harmony_3 <- scRNA_seq_expt_comb_umap_3 %>% RunHarmony("expt", plot_convergence = TRUE)
rm(scRNA_seq_expt_comb_umap_3)

harmony_embeddings_3 <- Embeddings(scRNA_seq_expt_comb_harmony_3, 'harmony')

scRNA_seq_expt_comb_harmony_umap_3 <- scRNA_seq_expt_comb_harmony_3 %>%
  RunUMAP(reduction = "harmony", dims = 1:14) %>%
  FindNeighbors(reduction = "harmony", dims = 1:14) %>%
  FindClusters(resolution = 0.4) %>%
  identity()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/uMAPS/Fibroblast_uMAP_Subcluster.pdf", height=10, width=10)
DimPlot(scRNA_seq_expt_comb_harmony_umap_3, reduction = "umap", cols = c("springgreen4","goldenrod2","darkblue","turquoise","purple", "firebrick4", "deepskyblue", "lightblue2", "lightpink", "darkorange", "violetred"),
        label = TRUE, pt.size = .1, label.size = 9) + theme(text = element_text(size=20, face = "bold"))
dev.off()

#Reordering idents to put young first in visualization
#Found here:https://github.com/satijalab/seurat/issues/2471
scRNA_seq_expt_comb_harmony_umap_3$orig.ident <- factor(x = scRNA_seq_expt_comb_harmony_umap_3$orig.ident, levels = c("mm10_3mths", "mm10_18mths"))

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/uMAPS/Stroma_Aging_Split.pdf", height=10, width=20)
scRNA_seq_expt_comb_harmony_umap_3$orig.ident <- factor(x = scRNA_seq_expt_comb_harmony_umap_3$orig.ident, levels = c("mm10_3mths", "mm10_18mths"))
DimPlot(scRNA_seq_expt_comb_harmony_umap_3, reduction = "umap", cols = c("springgreen4","goldenrod2","darkblue","turquoise","purple", "firebrick4", "deepskyblue", "lightblue2", "lightpink", "darkorange", "violetred"), split.by = "orig.ident", label = TRUE, label.size = 9) + theme(text = element_text(size=20)) 
dev.off()

DimPlot(scRNA_seq_expt_comb_harmony_umap_3, reduction = "umap", split.by = "replicate", label = TRUE, label.size = 9) + theme(text = element_text(size=20, face = "bold")) 

VlnPlot(scRNA_seq_expt_comb_harmony_umap_3, features = "nFeature_RNA") 

#saveRDS(scRNA_seq_expt_comb_harmony_umap_3, file = "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_harmony_umap_third_stroma_subset_05192022_BA.rds")
scRNA_seq_expt_comb_harmony_umap_3 <- readRDS("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_harmony_umap_third_stroma_subset_05192022_BA.rds")

Meta_Fibro <- scRNA_seq_expt_comb_harmony_umap_3@meta.data
#write.csv(Meta_Fibro , "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/Metadata_Fibroblasts.csv")

### Visualization 
#new.cluster.ids <- c("Fibroblasts", "Fibroblasts", "Vascular", "Fibroblasts", "Fibroblasts", "Pericytes", "Fibroblasts", 
#                     "Fibroblasts", "Vascular", "Vascular", "Fibroblasts")
#seurat.cluster.ids <- scRNA_seq_expt_comb_harmony_umap_3$seurat_clusters
#names(new.cluster.ids) <- levels(scRNA_seq_expt_comb_harmony_umap_3)
#scRNA_seq_expt_comb_harmony_umap_3 <- RenameIdents(scRNA_seq_expt_comb_harmony_umap_3,new.cluster.ids)

#pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/uMAPS/Fibroblast_uMAP.pdf", height=10, width=10)
#DimPlot(scRNA_seq_expt_comb_harmony_umap_3, cols = c("#22BDC1","#8BB6E1","#292562"), reduction = "umap", label = TRUE, pt.size = .1, label.size = 9) + theme(text = element_text(size=20)) 
#dev.off()

dim(scRNA_seq_expt_comb_harmony_umap_3@meta.data)   #3832 cells
table(Idents(scRNA_seq_expt_comb_harmony_umap_3))
#0    6    1    4    3   10    7    5    2    8    9 
#1068  143  798  380  429   22  141  226  530   52   43 
table(scRNA_seq_expt_comb_harmony_umap_3$orig.ident)
#mm10_3mths mm10_18mths 
#1090        2742

#scRNA_seq_expt_comb_harmony_umap_3$agewithcluster <- paste(scRNA_seq_expt_comb_harmony_umap_3$seurat_clusters, scRNA_seq_expt_comb_harmony_umap_3$orig.ident, sep = "_")
#table <- as.data.frame(table(scRNA_seq_expt_comb_harmony_umap_3$agewithcluster))
#write.csv(table, "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/Cell_counts_per_cluster_fibro.csv")
#scRNA_seq_expt_comb_harmony_umap_3$agewithcluster <- paste(scRNA_seq_expt_comb_harmony_umap_3$seurat_clusters, scRNA_seq_expt_comb_harmony_umap_3$orig.ident, sep = "_")

#scRNA_seq_expt_comb_harmony_umap_3$agewithreplicate <- paste(scRNA_seq_expt_comb_harmony_umap_3$seurat_clusters, scRNA_seq_expt_comb_harmony_umap_3$replicate, sep = "_")
#table <- as.data.frame(table(scRNA_seq_expt_comb_harmony_umap_3$agewithreplicate))
#write.csv(table, "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/Cell_counts_per_replicate_fibro.csv")

##Identifying immune and lymphatic clusters
FeaturePlot(scRNA_seq_expt_comb_harmony_umap_3, features = c("Cd3d", "Ptprc", "Ms4a1"))
FeaturePlot(scRNA_seq_expt_comb_harmony_umap_3, features = c("Pecam1", "Pdpn", "Prox1", "Lyve1", "Itga2b"))

#################################################
##Plotting Single Cell Heatmap
my_levels <- c("0", "6", "1", "4", "3", "7",  "5", "10", "2", "8", "9")
scRNA_seq_expt_comb_harmony_umap_3@active.ident <- factor(x = scRNA_seq_expt_comb_harmony_umap_3@active.ident, levels = my_levels)

scRNA_seq_expt_comb_markers_3 <- FindAllMarkers(scRNA_seq_expt_comb_harmony_umap_3, only.pos = FALSE, min.pct = 0.25, test.use = 'LR', latent.vars = 'expt', logfc.threshold = 0.25)
write.csv(scRNA_seq_expt_comb_markers_3, "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/Table_markers_all_prelim_pos_neg_fibro3_01092023.csv")

scRNA_seq_expt_comb_markers_3 <- read.csv("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/Table_markers_all_prelim_pos_neg_fibro3_01092023.csv")

top5 <- scRNA_seq_expt_comb_markers_3 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/Heatmap_top_expressedwithlegend.pdf", height=9, width=7)
DoHeatmap(scRNA_seq_expt_comb_harmony_umap_3, features = top5$gene, raster = F, lines.width =5, 
          group.colors = c("springgreen4", "deepskyblue", "goldenrod2", "purple", "turquoise", "lightblue2", "firebrick4", 
                           "violetred", "darkblue", "lightpink", "darkorange"))+ theme(text = element_text(size=5))
dev.off()

########################
####Cell Proportion Plots####
##Using Sid's code to plot log2FC os cell pop changes

Meta <- scRNA_seq_expt_comb_harmony_umap_3@meta.data
#Meta$celltype.age <- paste(Idents(scRNA_seq_expt_comb_harmony_umap_2), scRNA_seq_expt_comb_harmony_umap_2$orig.ident, sep = "_")
#Meta$celltype.age.rep <- paste(Meta$celltype.age, Meta$replicate, sep = "_")
Meta.2 <- Meta[,c(1,4,8)]
Meta.2 <- Meta.2 %>% group_by_all() %>% summarise(COUNT = n())
#write.csv(Meta.2,"C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/Fibro_Cell_Count_122022.csv")
Meta.2 <- read.csv("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/Fibro_Cell_Count_122022_BA.csv")
Meta.2$replicate <- gsub("M3","M18",Meta.2$replicate)
NewMeta = dcast(Meta.2, replicate + RNA_snn_res.0.4 ~ orig.ident, sum, value.var = "COUNT_New") 
NewMeta$Ratio <- log2(NewMeta$mm10_18mths/NewMeta$mm10_3mths)
NewMeta <- NewMeta[,c(1,2,5)]
NewMeta$RNA_snn_res.0.4 <- as.character(NewMeta$RNA_snn_res.0.4)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep1","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep2","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep3","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep4","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep5","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep6","Base",0)
NewMeta$RNA_snn_res.0.4 <- as.factor(NewMeta$RNA_snn_res.0.4)
NewMeta$Ratio <- as.numeric(NewMeta$Ratio)

colr <- c("springgreen4","goldenrod2", "violetred", "darkblue","turquoise","purple", "firebrick4", "deepskyblue", "lightblue2", "lightpink", "darkorange", "black")

p1 <- ggplot(NewMeta,aes(x = reorder(RNA_snn_res.0.4, -Ratio), y = Ratio, fill = RNA_snn_res.0.4)) + 
  geom_boxplot()+
  scale_fill_manual(values=(colr)) +
  xlab("CellTypes") +
  ylab("log2 FC (Old/Young)") + 
  ggtitle("Cell Type Proportions") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.15, binwidth = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  #ylim(c(-2,5)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  stat_compare_means(method = "anova", label.y = 10) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = "Base", hide.ns = T, label.y = 8)   
p1                                                                            
pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/Fibroblasts_RNA_Pvals.pdf", height = 7, width = 8)
p1
dev.off()

##########################################
Meta %>% dplyr::filter(Idents(scRNA_seq_expt_comb_harmony_umap_3) == "0"| Idents(scRNA_seq_expt_comb_harmony_umap_3) == "6"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_3) == "1"| Idents(scRNA_seq_expt_comb_harmony_umap_3) == "4"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_3) == "3"| Idents(scRNA_seq_expt_comb_harmony_umap_3) == "10"|
                         Idents(scRNA_seq_expt_comb_harmony_umap_3) == "7") -> Fibro_Cells
dim(Fibro_Cells)
cell_to_keep_Fib <- row.names(Fibro_Cells)
length(cell_to_keep_Fib)

subset_Fibro = subset(scRNA_seq_expt_comb_harmony_umap_3, cells = cell_to_keep_Fib)
table(Idents(subset_Fibro))

###########
#Based on Djamel's code to make cell proportion plots
###Create a proportions table
y=prop.table(x=table(subset_Fibro$orig.ident, as.character(subset_Fibro$seurat_clusters)), margin=2)

###Adding proportions for all of the cells
nbr_row <- length(unique(subset_Fibro$orig.ident))
ALL <- vector("numeric", nbr_row) # prepare a container

for (i in 1:nbr_row) {
  res <- table(subset_Fibro$orig.ident)[[i]]/dim(subset_Fibro)[2]
  ALL[i] <- res         # change to assignment
}
ALL

y1 <- cbind(y, All=ALL)
y1

#Change to long version of the table and rename the columns
df1 <- melt(y1)
head(df1)
length(unique(subset_Fibro$seurat_clusters))
colnames(df1) <- c("orig.ident", "Clusters", "Percentage")
head(df1)

#To order the clusters
ordered_clusters <- c('All', '0', '6', '1', '4', '3','10','7')

#To color the ages
col_orig.ident = c('mm10_3mths'= "#9e9ac8",
                   'mm10_18mths'= "#66c2a4")

##Plot with ggplot
p_Age <- ggplot(data=df1, aes(x=as.character(Clusters), y=Percentage, 
                              fill=factor(orig.ident, levels=c('mm10_3mths',"mm10_18mths"))))+ 
  geom_bar(stat="identity", color="black", position =position_fill(reverse = TRUE)) +
  scale_fill_manual(values=col_orig.ident) + 
  scale_x_discrete(limits=ordered_clusters) +
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y =  element_text(size = 20))+
  theme(axis.title = element_text(size = 20))+
  labs(fill='orig.ident') + 
  ylab("Proportion of cells") +
  xlab("Clusters") +
  ggtitle("orig.ident")+
  geom_hline(yintercept = 0.2683663,linetype="dotted", size=1) 

print(p_Age)

ggsave("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/Fibroblast_Prop.pdf", 
       p_Age , width=3, height=1.7,  units="in", scale=3)

########################
pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/Dotplot_Major_Markers.pdf", height=20, width=20)
DotPlot(scRNA_seq_expt_comb_harmony_umap_3, features = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Pdgfra", 
                                                         "Acta2", "Myl9", "Mylk", "Myh11",
                                                         "Rgs5", "Des", "Notch3", "Ndufa4l2", "Higd1b",
                                                         "Pax7", "Bmp4",
                                                         "Pecam1", "Cdh5", "Eng", "Sox17", "Sele",
                                                         "Prox1", "Lyve1", "Itga2b"),
        dot.scale = 5) + coord_flip() + theme(text = element_text(size=12))
dev.off()

Meta_6 <- scRNA_seq_expt_comb_harmony_umap_3@meta.data
Meta_6 %>% dplyr::filter(Idents(scRNA_seq_expt_comb_harmony_umap_3) == "6"| Idents(scRNA_seq_expt_comb_harmony_umap_3) == "0"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_3) == "1"| Idents(scRNA_seq_expt_comb_harmony_umap_3) == "4"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_3) == "3"| Idents(scRNA_seq_expt_comb_harmony_umap_3) == "10"| 
                                                                                     Idents(scRNA_seq_expt_comb_harmony_umap_3) == "7") -> Fibroblast
dim(Fibroblast)
cell_to_keep <- row.names(Fibroblast)
subset_Fibro = subset(scRNA_seq_expt_comb_harmony_umap_3, cells = cell_to_keep)
my_levels <- c("0", "6", "1", "4", "3", "10", "7")
subset_Fibro@active.ident <- factor(x = subset_Fibro@active.ident, levels = my_levels)

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Manuscript/Figures/Figure 3_Fibroblasts/uMAPS/Fibroblast_uMAP_Subcluster_Fibroblast_Only.pdf", height=10, width=10)
DimPlot(subset_Fibro, reduction = "umap", cols = c("springgreen4", "deepskyblue", "goldenrod2", "purple", "turquoise", "violetred", "lightblue2"),
        label = TRUE, pt.size = .1, label.size = 9) + theme(text = element_text(size=20)) + xlim(c(-7,10))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/Dotplot_Revised.pdf", height=20, width=20)
DotPlot(subset_Fibro, features = c("Penk", "Fbln1","Tnc", "Fmo1", "Postn", "Mfap4", "Spon1","Crabp1","Col15a1", 
                                   "Col6a3", "Col5a3", "Col4a1", "Col4a2", "Lpl", "Fabp4", "Pparg", "Fap",
                                   "Gdf10", "F3", "Hmcn2", "Gpc3", "Dpt", "Ly6a",
                                   "Pi16", "Dpp4", "Ly6c1", "Tek", "Aldh1a3","Fn1", "Col14a1","Sema3c", "Esr1",
                                   "Ccl19", "Cxcl12", "Cd74"),
        dot.scale = 5) + coord_flip() + theme(text = element_text(size=12))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/uMAPS/Fibroblasts_Pi16.pdf", height = 5, width = 5)
FeaturePlot(subset_Fibro, features = c("Pi16"), cols = c("lightgrey", "red3"), label = F, slot = "scale.data") +xlim(c(-7,10))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/uMAPS/Fibroblasts_Col15a1.pdf", height = 5, width = 5)
FeaturePlot(subset_Fibro, features = c("Col15a1"), cols = c("lightgrey", "red3"), label = F, slot = "scale.data") +xlim(c(-7,10))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/uMAPS/Fibroblasts_Ccl19.pdf", height = 5, width = 5)
FeaturePlot(subset_Fibro, features = c("Ccl19"), cols = c("lightgrey", "red3"), label = F, slot = "scale.data") +xlim(c(-7,10))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/uMAPS/Fibroblasts_Cxcl12.pdf", height = 5, width = 5)
FeaturePlot(subset_Fibro, features = c("Cxcl12"), cols = c("lightgrey", "red3"), label = F, slot = "scale.data") +xlim(c(-7,10))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/uMAPS/Fibroblasts_Fn1.pdf", height = 5, width = 5)
FeaturePlot(subset_Fibro, features = c("Fn1"), cols = c("lightgrey", "red3"), label = F, slot = "scale.data") +xlim(c(-7,10))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/uMAPS/Fibroblasts_Col5a3.pdf", height = 5, width = 5)
FeaturePlot(subset_Fibro, features = c("Col5a3"), cols = c("lightgrey", "red3"), label = F, slot = "scale.data") +xlim(c(-7,10))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/uMAPS/Fibroblasts_Fap.pdf", height = 5, width = 5)
FeaturePlot(subset_Fibro, features = c("Fap"), cols = c("lightgrey", "red3"), label = F, slot = "scale.data") +xlim(c(-7,10))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/uMAPS/Fibroblasts_Fabp4.pdf", height = 5, width = 5)
FeaturePlot(subset_Fibro, features = c("Fabp4"), cols = c("lightgrey", "red3"), label = F, slot = "scale.data") +xlim(c(-7,10))
dev.off()

FeaturePlot(subset_Fibro, features = c("Sgip1"), cols = c("lightgrey", "red3"), label = F, slot = "scale.data") +xlim(c(-7,10))

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/uMAPS/Fibroblasts_Cdkn1a.pdf", height = 5, width = 10)
FeaturePlot(subset_Fibro, features = c("Cdkn1a"), split.by = "orig.ident", cols = c("lightgrey", "red3"), label = F, slot = "scale.data") 
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure 3_Fibroblasts/uMAPS/Fibroblasts_Cdkn2a.pdf", height = 5, width = 10)
FeaturePlot(subset_Fibro, features = c("Cdkn2a"), split.by = "orig.ident", cols = c("lightgrey", "red3"), label = F, slot = "scale.data") 
dev.off()

###########
df <- subset(subset_Fibro, subset = Cdkn1a > 0 & Cdkn2a > 0)
table(df@meta.data$seurat_clusters)
tbl <- table(df@meta.data$orig.ident)
pie(tbl, main = "Cdkn1a & Cdkn2a")             
tbl

df <- subset(subset_Fibro, subset = Cdkn1a > 0)
table(df@meta.data$seurat_clusters)
tbl <- table(df@meta.data$orig.ident)
pie(tbl, main = "Cdkn1a")             
tbl

df <- subset(subset_Fibro, subset = Cdkn2a > 0)
table(df@meta.data$seurat_clusters)
tbl <- table(df@meta.data$orig.ident)
pie(tbl, main = "Cdkn2a")             
tbl


