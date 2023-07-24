##Subclustering Fibroblasts
##04112022
##Based off Neerja's Original Clustering script 
##This is based off of the Seurat Tutorial - https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
##Harmony Vignette: https://portals.broadinstitute.org/harmony/SeuratV3.html

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
# "Sox17", "Mmrn1", "Prox1", "Flt4", "Ccl21a"
#https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000704

##########################################
###Removing muscle contamination and immune doublets (corresponding to clusters 16 and 10 from previous analysis)
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

##########################################
###Removing low quality cells (corresponds to clusters 6 and 10 from previous analysis)
###Removing T and B cells
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

#saveRDS(scRNA_seq_expt_comb_harmony_umap_3, file = "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_harmony_umap_third_stroma_subset_05192022_BA.rds")
scRNA_seq_expt_comb_harmony_umap_3 <- readRDS("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_harmony_umap_third_stroma_subset_05192022_BA.rds")

