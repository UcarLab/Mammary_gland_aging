##Subclustering Epithelial Cells
## Revised 12092022
##Based off Neerja's Clustering script 
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
dim(BA_obj@meta.data) #48180 cells
head(BA_obj@meta.data)
table(Idents(BA_obj)) #Luminal-AV 2843
                      #Luminal-HS 1494
                      #Myoepithelial #2971

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

###Subsetting the epithelial cells and generating a list of the cell identifiers
Meta %>% dplyr::filter(Idents(BA_obj) == 'Luminal-HS'| Idents(BA_obj) == 'Luminal-AV'| Idents(BA_obj) == 'Myoepithelial' ) -> Epi
dim(Epi) #7308
cell_to_keep <- row.names(Epi)
length(cell_to_keep)
#write.csv(cell_to_keep, "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/epithelial_cell_subset_06112021.csv")
rm(BA_obj)

###Loading the unfiltered/unlabeled RDS file
scRNA_seq_expt_comb_before_filter <- readRDS("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_before_filter_May26_2021.rds")

###Subsetting the unfiltered RDS file for the epithelial cells
subset_epi_combined = subset(scRNA_seq_expt_comb_before_filter, cells = cell_to_keep)
dim(subset_epi_combined@meta.data) 
head(subset_epi_combined@meta.data)
rm(scRNA_seq_expt_comb_before_filter)

#saveRDS(subset_epi_combined, file = "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/subset_epi_combined_June11_2021_BA.rds")
#subset_epi_combined <- readRDS("C:/Users/angarb/Box/AAnczukow_lab_NGS_raw_data/ging single cells/Subclustering/Single_Cell_June2021/subset_epi_combined_June11_2021_BA.rds")

##Check QC (we have prefiltered these for nFeature_RNA > 500 & percent.mt < 10)
plot1 <- FeatureScatter(subset_epi_combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(subset_epi_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

##Re-run QC just in case
scRNA_seq_expt_comb_filtered <- subset(subset_epi_combined, subset = nFeature_RNA > 500 & percent.mt < 10)
rm(subset_epi_combined)

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
scRNA_seq_expt_comb_neighbors <- FindNeighbors(scRNA_seq_expt_comb_pca, dims = 1:10)
rm(scRNA_seq_expt_comb_pca)
scRNA_seq_expt_comb_clusters <- FindClusters(scRNA_seq_expt_comb_neighbors, resolution = 0.8)
rm(scRNA_seq_expt_comb_neighbors)

#Run non-linear dimensional reduction (UMAP/tSNE)
scRNA_seq_expt_comb_umap <- RunUMAP(scRNA_seq_expt_comb_clusters, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
DimPlot(scRNA_seq_expt_comb_umap, reduction = "umap", label=TRUE)
rm(scRNA_seq_expt_comb_clusters)

#Passing Seurat object to Harmony and integrating out experiment
#Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.RNA.harmony; see ?make.names for more details on syntax validity
#Took 5/10 iterations
scRNA_seq_expt_comb_harmony <- scRNA_seq_expt_comb_umap %>% RunHarmony("expt", plot_convergence = TRUE)
rm(scRNA_seq_expt_comb_umap)

#To access embeddings
harmony_embeddings <- Embeddings(scRNA_seq_expt_comb_harmony, 'harmony')

scRNA_seq_expt_comb_harmony_umap <- scRNA_seq_expt_comb_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:10) %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>%
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
#saveRDS(scRNA_seq_expt_comb_harmony_umap, file = "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_harmony_umap_initial_epi_subset_06102022_BA.rds")
#scRNA_seq_expt_comb_harmony_umap <- readRDS("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_harmony_umap_initial_epi_subset_06102022_BA.rds")

##Identifying immune clusters
FeaturePlot(scRNA_seq_expt_comb_harmony_umap, features = c("Epcam", "Ptprc", "Cd52", "Blnk", "Itgax", "Cd3d", "Cd4", "Cd8a", "Ms4a1", "Cd79a", "Jchain"))
VlnPlot(scRNA_seq_expt_comb_harmony_umap, features = "nFeature_RNA") #Cluster 10 seems to be of low quality

##########################################
##########################################
###Removing immune doublets and low quality cells

###Creating metadata object and verifying it looks as expected
Meta_2 <- scRNA_seq_expt_comb_harmony_umap@meta.data
dim(Meta_2) #7308

Meta_2 %>% dplyr::filter(Idents(scRNA_seq_expt_comb_harmony_umap) == "9"| Idents(scRNA_seq_expt_comb_harmony_umap) == "13") -> Doublets
dim(Doublets)
cell_to_remove <- row.names(Doublets)
length(cell_to_remove) #221

final_cells_to_keep <- setdiff(cell_to_keep, cell_to_remove)
length(final_cells_to_keep) #7087

Meta_Q <- scRNA_seq_expt_comb_harmony_umap@meta.data
Meta_Q %>% dplyr::filter(scRNA_seq_expt_comb_harmony_umap$nFeature_RNA < 1000) -> RNA #134

###Loading the unfiltered/unlabeled RDS file
scRNA_seq_expt_comb_before_filter <- readRDS("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_before_filter_May26_2021.rds")

###Subsetting the unfiltered RDS file for the epithelial cells
subset_epi_combined_2 = subset(scRNA_seq_expt_comb_before_filter, cells = final_cells_to_keep)
table(subset_epi_combined_2$orig.ident)
#mm10_18mths  mm10_3mths 
#4011        3076

##Run QC (just in case)with increased filter
scRNA_seq_expt_comb_filtered_2 <- subset(subset_epi_combined_2, subset = nFeature_RNA > 1000 & percent.mt < 10)
rm(subset_epi_combined_2)

table(scRNA_seq_expt_comb_filtered_2$orig.ident)
#mm10_18mths  mm10_3mths 
#3947        3006 
#6953

#Normalizing the data
scRNA_seq_expt_comb_norm_2 <- NormalizeData(scRNA_seq_expt_comb_filtered_2, normalization.method = "LogNormalize", scale.factor = 10000)
rm(scRNA_seq_expt_comb_filtered_2)

#Identification of highly variable features (feature selection)
scRNA_seq_expt_comb_hvf_2 <- FindVariableFeatures(scRNA_seq_expt_comb_norm_2, selection.method = "vst", nfeatures = length(rownames(scRNA_seq_expt_comb_norm_2)))

all.genes_2 <- rownames(scRNA_seq_expt_comb_hvf_2)
scRNA_seq_expt_comb_scal_2 <- ScaleData(scRNA_seq_expt_comb_hvf_2, features = all.genes_2)
rm(scRNA_seq_expt_comb_hvf_2)

#Perform linear dimensional reduction
scRNA_seq_expt_comb_pca_2 <- RunPCA(scRNA_seq_expt_comb_scal_2, features = VariableFeatures(object = scRNA_seq_expt_comb_scal_2))
rm(scRNA_seq_expt_comb_scal_2)

ElbowPlot(scRNA_seq_expt_comb_pca_2)

#Cluster the cells
scRNA_seq_expt_comb_neighbors_2 <- FindNeighbors(scRNA_seq_expt_comb_pca_2, dims = 1:14)
rm(scRNA_seq_expt_comb_pca_2)

scRNA_seq_expt_comb_clusters_2 <- FindClusters(scRNA_seq_expt_comb_neighbors_2, resolution = 0.8)
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
  FindClusters(resolution = 0.8) %>%
  identity()

DimPlot(scRNA_seq_expt_comb_harmony_umap_2, reduction = "umap", label = TRUE, pt.size = .1, label.size = 9) + theme(text = element_text(size=20, face = "bold")) 
#DimPlot(scRNA_seq_expt_comb_harmony_umap_2, reduction = "umap", group.by = "expt", label = TRUE)
DimPlot(scRNA_seq_expt_comb_harmony_umap_2, reduction = "umap", group.by = "orig.ident", label = TRUE)

new.cluster.ids <- c("C8", "C4", "C9", "C5", "C1", "C6", "C2", "C11", "C7", "C3", "C10")

names(new.cluster.ids) <- levels(scRNA_seq_expt_comb_harmony_umap_2)
scRNA_seq_expt_comb_harmony_umap_2 <- RenameIdents(scRNA_seq_expt_comb_harmony_umap_2,new.cluster.ids)

scRNA_seq_expt_comb_harmony_umap_2$new_cluster_ids <- paste(Idents(scRNA_seq_expt_comb_harmony_umap_2))
head(scRNA_seq_expt_comb_harmony_umap_2@meta.data)

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/Age_uMAP/Epithelial_Aging_Recolor.pdf", height=10, width=20)
DimPlot(scRNA_seq_expt_comb_harmony_umap_2, reduction = "umap", cols = c("#BD9ED8", "#0C911F", "#F8B4C0", "#B2EDD1", 
                                                                           "#3ABCED", "#0AC48F", "#B3DFED","#BC46B6", 
                                                                           "#153827", "#2C73BA", "#46206B"), label = TRUE, pt.size = .1, label.size = 9) + theme(text = element_text(size=20, face = "bold")) 
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/Age_uMAP/Epithelial_Aging_Split_Recolor.pdf", height=10, width=20)
scRNA_seq_expt_comb_harmony_umap_2$orig.ident <- factor(x = scRNA_seq_expt_comb_harmony_umap_2$orig.ident, levels = c("mm10_3mths", "mm10_18mths"))
DimPlot(scRNA_seq_expt_comb_harmony_umap_2, reduction = "umap", split.by = "orig.ident", cols = c("#BD9ED8", "#0C911F", "#F8B4C0", "#B2EDD1", 
                                                                         "#3ABCED", "#0AC48F", "#B3DFED","#BC46B6", 
                                                                         "#153827", "#2C73BA", "#46206B"), label = TRUE, pt.size = .1, label.size = 9) + theme(text = element_text(size=20, face = "bold")) 
dev.off()

#saveRDS(scRNA_seq_expt_comb_harmony_umap_2, file = "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_harmony_umap_second_epi_subset_12092022_BA_PC14Res0.8.rds")
scRNA_seq_expt_comb_harmony_umap_2 <-readRDS("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_harmony_umap_second_epi_subset_12092022_BA_PC14Res0.8.rds")

# Visualization for Figure 2A
#new.cluster.ids <- c("Myoepithelial", "Luminal-AV", "Myoepithelial", "Luminal-AV", "Luminal-HS", "Luminal-AV", "Luminal-HS", "Myoepithelial", "Luminal-AV", "Luminal-HS", "Myoepithelial")
#seurat.cluster.ids <- scRNA_seq_expt_comb_harmony_umap_2$seurat_clusters
#names(new.cluster.ids) <- levels(scRNA_seq_expt_comb_harmony_umap_2)
#scRNA_seq_expt_comb_harmony_umap_2 <- RenameIdents(scRNA_seq_expt_comb_harmony_umap_2,new.cluster.ids)
#DimPlot(scRNA_seq_expt_comb_harmony_umap_2, cols = c("#3DB54A","#A0D082","#09713A"), reduction = "umap", label = TRUE, pt.size = .1, label.size = 9) + theme(text = element_text(size=20, face = "bold")) 
#pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure2 epithelial/Epithelial_Aging_Split.pdf", height=10, width=20)
#DimPlot(scRNA_seq_expt_comb_harmony_umap, cols = c("#3DB54A","#A0D082","#09713A"), reduction = "umap", split.by = "orig.ident", label = F)
#dev.off()

########################
####Cell Proportion Plots####
##Using Sid's code to plot log2FC os cell pop changes

Meta <- scRNA_seq_expt_comb_harmony_umap_2@meta.data
#Meta$celltype.age <- paste(Idents(scRNA_seq_expt_comb_harmony_umap_2), scRNA_seq_expt_comb_harmony_umap_2$orig.ident, sep = "_")
#Meta$celltype.age.rep <- paste(Meta$celltype.age, Meta$replicate, sep = "_")
Meta.2 <- Meta[,c(1,4,9)]
Meta.2 <- Meta.2 %>% group_by_all() %>% summarise(COUNT = n())
#write.csv(Meta.2,"C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/Epi_Cell_Count_122022.csv")
Meta.2 <- read.csv("C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/Epi_Cell_Count_122022_BA.csv")
Meta.2$replicate <- gsub("M3","M18",Meta.2$replicate)
NewMeta = dcast(Meta.2, replicate + new_cluster_ids ~ orig.ident, sum, value.var = "COUNT_New") 
NewMeta$Ratio <- log2(NewMeta$mm10_18mths/NewMeta$mm10_3mths)
NewMeta <- NewMeta[,c(1,2,5)]
NewMeta$new_cluster_ids <- as.character(NewMeta$new_cluster_ids)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep1","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep2","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep3","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep4","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep5","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep6","Base",0)
NewMeta$new_cluster_ids <- as.factor(NewMeta$new_cluster_ids)
NewMeta$Ratio <- as.numeric(NewMeta$Ratio)

#colr1 <- c("orangered1","deepskyblue","lightblue2","darkblue","purple","dodgerblue3", "goldenrod2", "firebrick4", "lightpink", "springgreen4")
colr <- c("springgreen4", "#3ABCED", "#46206B", "#BC46B6", "#B3DFED", "#2C73BA", "#0C911F", "#B2EDD1", "#0AC48F",  "#153827", "#BD9ED8", "#F8B4C0")
         
p1 <- ggplot(NewMeta,aes(x = reorder(new_cluster_ids, -Ratio), y = Ratio, fill = new_cluster_ids)) + 
  geom_boxplot()+
  scale_fill_manual(values=(colr)) +
  xlab("CellTypes") +
  ylab("log2 FC (Old/Young)") + 
  ggtitle("Cell Type Proportions") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, binwidth = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  stat_compare_means(method = "anova", label.y = 10) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = "Base", hide.ns = T, label.y = 8)   
p1                                                                            
pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/EpithelialCells_RNA_Pvals.pdf", height = 7, width = 8)
p1
dev.off()

#################################################
##Plotting Single Cell Heatmap

my_levels <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11")
scRNA_seq_expt_comb_harmony_umap_2@active.ident <- factor(x = scRNA_seq_expt_comb_harmony_umap_2@active.ident, levels = my_levels)

scRNA_seq_expt_comb_markers <- FindAllMarkers(scRNA_seq_expt_comb_harmony_umap_2, only.pos = FALSE, min.pct = 0.25, test.use = 'LR', latent.vars = 'expt', logfc.threshold = 0.25)
write.csv(scRNA_seq_expt_comb_markers, "C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/Table_markers_all_prelim_pos_neg_12142022.csv")

scRNA_seq_expt_comb_markers <- read.csv("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/Top5_Heatmap/Table_markers_all_prelim_pos_neg_12142022.csv")

top5 <- scRNA_seq_expt_comb_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/Top5_Heatmap/Heatmap_top_expressedwithlegend.pdf", height=9, width=4)
DoHeatmap(scRNA_seq_expt_comb_harmony_umap_2, features = top5$gene, raster = F, lines.width =5, 
          group.colors = c("#3ABCED", "#B3DFED", "#2C73BA", "#0C911F", "#B2EDD1", "#0AC48F",  "#153827", "#BD9ED8", "#F8B4C0", "#46206B", "#BC46B6"))+ theme(text = element_text(size=5))
dev.off()

#####################################
#####Luminal Markers

table(Idents(scRNA_seq_expt_comb_harmony_umap_2))
Meta <- scRNA_seq_expt_comb_harmony_umap_2@meta.data
Meta %>% dplyr::filter(Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C1"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C2"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C3"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C4"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C5"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C6"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C7") -> Luminal_Cells
dim(Luminal_Cells)
cell_to_keep <- row.names(Luminal_Cells)
length(cell_to_keep)

subset_Luminal = subset(scRNA_seq_expt_comb_harmony_umap_2, cells = cell_to_keep)
table(Idents(subset_Luminal))

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/Dotplot_Luminal.pdf", height=20, width=20)
DotPlot(subset_Luminal, features = c("Epcam", "Prlr", "Cited1", "Pgr", "Prom1", "Esr1", "Ptn", "Gpx3", "Wfdc2", "Ly6a",
                                                         "Rcan1", "Fxyd2", "Tph1", "Fam3c", "Capsl",
                                                         "Kit", "Aldh1a3", "Cd14", "Elf5", 
                                                         "Mfge8", "Trf", "Csn3", "Wfdc18", "Ltf", "Cst3", "Igfbp5",
                                                         "Hey1", "Txnip", "Crip1", "Alox15", "Alox12e", "Palmd", "Rspo1",
                                                         "Mcm2", "Mcm3", "Mcm5", "Cdk1", "Mki67", "Top2a", "E2f1", "Hmgb2", "Tubb5", "Stmn1"),
        dot.scale = 5) + coord_flip() + theme(text = element_text(size=12))
dev.off()


#####
Meta %>% dplyr::filter(Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C1"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C2"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C3") -> Luminal_HS_Cells
dim(Luminal_HS_Cells)
cell_to_keep_HS <- row.names(Luminal_HS_Cells)
length(cell_to_keep_HS)

subset_Luminal_HS = subset(scRNA_seq_expt_comb_harmony_umap_2, cells = cell_to_keep_HS)
table(Idents(subset_Luminal_HS))

##################
#Based on Djamel's code to make cell proportion plots
###Create a proportions table
y=prop.table(x=table(subset_Luminal_HS$orig.ident, as.character(subset_Luminal_HS$new_cluster_ids)), margin=2)

###Adding proportions for all of the cells
nbr_row <- length(unique(subset_Luminal_HS$orig.ident))
ALL <- vector("numeric", nbr_row) # prepare a container

for (i in 1:nbr_row) {
  res <- table(subset_Luminal_HS$orig.ident)[[i]]/dim(subset_Luminal_HS)[2]
  ALL[i] <- res         # change to assignment
}
ALL

y1 <- cbind(y, All=ALL)
y1

#Change to long version of the table and rename the columns
df1 <- melt(y1)
head(df1)
length(unique(subset_Luminal_HS$new_cluster_ids))
colnames(df1) <- c("orig.ident", "Clusters", "Percentage")
head(df1)

#To order the clusters
ordered_clusters <- c('All', 'C1', 'C2', 'C3')

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
  geom_hline(yintercept = 0.6525856,linetype="dotted", size=1) 

print(p_Age)

ggsave("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/LuminalHS_Prop.pdf", 
       p_Age , width=3, height=1.7,  units="in", scale=3)

###########################
###Plotting uMAP and Feature Plots

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPs/Luminal_HS/HS_UMAP.pdf", height=20, width=20)
DimPlot(subset_Luminal_HS, reduction = "umap", cols = c("#3ABCED", "#B3DFED", "#2C73BA"), 
        label = TRUE, label.size = 7) + theme(text = element_text(size=15)) + ylim(c(9.5,14.5)) +xlim(c(-6.5,-0.5))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPs/Luminal_HS/LuminalHS_Fxyd2.pdf", height=20, width=20)
FeaturePlot(subset_Luminal_HS, features = c("Fxyd2"), cols = c("lightgrey", "red3"), label = F, ncol = 1, slot = "scale.data") + ylim(c(9.5,14.5)) + xlim(c(-6.5,-0.5))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPs/Luminal_HS/LuminalHS_Tph1.pdf", height=20, width=20)
FeaturePlot(subset_Luminal_HS, features = c("Tph1"), cols = c("lightgrey", "red3"), label = F, ncol = 1, slot = "scale.data") + ylim(c(9.5,14.5)) + xlim(c(-6.5,-0.5))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPs/Luminal_HS/LuminalHS_Rcan1.pdf", height=20, width=20)
FeaturePlot(subset_Luminal_HS, features = c("Rcan1"), cols = c("lightgrey", "red3"), label = F, ncol = 1, slot = "scale.data") + ylim(c(9.5,14.5)) + xlim(c(-6.5,-0.5))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPs/Luminal_HS/LuminalHS_Kit.pdf", height=20, width=20)
FeaturePlot(subset_Luminal_HS, features = c("Kit"), cols = c("lightgrey", "red3"), label = F, ncol = 1, slot = "scale.data") + ylim(c(9.5,14.5)) + xlim(c(-6.5,-0.5))
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPs/Luminal_HS/LuminalHS_Wfdc18.pdf", height=20, width=20)
FeaturePlot(subset_Luminal_HS, features = c("Wfdc18"), cols = c("lightgrey", "red3"), label = F, ncol = 1, slot = "scale.data") + ylim(c(9.5,14.5)) + xlim(c(-6.5,-0.5))
dev.off()

#############
#####
Meta %>% dplyr::filter(Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C4"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C5"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C6"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C7") -> Luminal_AV_Cells
dim(Luminal_AV_Cells)
cell_to_keep_AV <- row.names(Luminal_AV_Cells)
length(cell_to_keep_AV)

subset_Luminal_AV = subset(scRNA_seq_expt_comb_harmony_umap_2, cells = cell_to_keep_AV)
table(Idents(subset_Luminal_AV))

##################
#Based on Djamel's code to make cell proportion plots
###Create a proportions table
y=prop.table(x=table(subset_Luminal_AV$orig.ident, as.character(subset_Luminal_AV$new_cluster_ids)), margin=2)

###Adding proportions for all of the cells
nbr_row <- length(unique(subset_Luminal_AV$orig.ident))
ALL <- vector("numeric", nbr_row) # prepare a container

for (i in 1:nbr_row) {
  res <- table(subset_Luminal_AV$orig.ident)[[i]]/dim(subset_Luminal_AV)[2]
  ALL[i] <- res         # change to assignment
}
ALL

y1 <- cbind(y, All=ALL)
y1

#Change to long version of the table and rename the columns
df1 <- melt(y1)
head(df1)
length(unique(subset_Luminal_AV$new_cluster_ids))
colnames(df1) <- c("orig.ident", "Clusters", "Percentage")
head(df1)

#To order the clusters
ordered_clusters <- c('All', 'C4', 'C5', 'C6', 'C7')

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
  geom_hline(yintercept = 0.5493166,linetype="dotted", size=1) 

print(p_Age)

ggsave("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/LuminalAV_Prop.pdf", 
       p_Age , width=3, height=1.7,  units="in", scale=3)

###########################
###Plotting uMAP and Feature Plots

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPS/Luminal_AV/AV_UMAP.pdf", height=20, width=20)
DimPlot(subset_Luminal_AV, reduction = "umap", cols = c("#0C911F", "#B2EDD1", "#0AC48F", "#153827"), 
        label = TRUE, label.size = 7) + theme(text = element_text(size=15)) + ylim(c(-11,-2)) + xlim(c(-8.5,-2.0)) 
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPS/Luminal_AV/LuminalAV_Alox12e.pdf", height=20, width=20)
FeaturePlot(subset_Luminal_AV, features = c("Alox12e"), cols = c("lightgrey", "red3"), label = F, ncol = 1, slot = "scale.data") + ylim(c(-11,-2)) + xlim(c(-8.5,-2.0)) 
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPS/Luminal_AV/LuminalAV_Rspo1.pdf", height=20, width=20)
FeaturePlot(subset_Luminal_AV, features = c("Rspo1"), cols = c("lightgrey", "red3"), label = F, ncol = 1, slot = "scale.data") + ylim(c(-11,-2)) + xlim(c(-8.5,-2.0)) 
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPS/Luminal_AV/LuminalAV_Csn3.pdf", height=20, width=20)
FeaturePlot(subset_Luminal_AV, features = c("Csn3"), cols = c("lightgrey", "red3"), label = F, ncol = 1, slot = "scale.data") + ylim(c(-11,-2)) + xlim(c(-8.5,-2.0)) 
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPS/Luminal_AV/LuminalAV_Mcm2.pdf", height=20, width=20)
FeaturePlot(subset_Luminal_AV, features = c("Mcm2"), cols = c("lightgrey", "red3"), label = F, ncol = 1, slot = "scale.data") + ylim(c(-11,-2)) + xlim(c(-8.5,-2.0)) 
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPS/Luminal_AV/LuminalAV_Stmn1.pdf", height=20, width=20)
FeaturePlot(subset_Luminal_AV, features = c("Stmn1"), cols = c("lightgrey", "red3"), label = F, ncol = 1, slot = "scale.data") + ylim(c(-11,-2)) + xlim(c(-8.5,-2.0))
dev.off()

#############
#####
Meta %>% dplyr::filter(Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C8"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C9"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C10"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C11") -> Myoep_Cells
dim(Myoep_Cells)
cell_to_keep_Myo <- row.names(Myoep_Cells)
length(cell_to_keep_Myo)

subset_Myo = subset(scRNA_seq_expt_comb_harmony_umap_2, cells = cell_to_keep_Myo)
table(Idents(subset_Myo))

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/Dotplot_Myoepithelial.pdf", height=20, width=20)
DotPlot(subset_Myo, features = c("Epcam", "Krt15", "Krt17", "Krt14", "Krt5", "Acta2", "Myl9", "Mylk", "Myh11",
                                     "Tagln", "Postn", "Actg2", 
                                     "Ccl2", "Cxcl10","Irgm1", "Stat1", "Tap1", "H2-Ab1", "H2-Eb1", "Cd74"),
        dot.scale = 5) + coord_flip() + theme(text = element_text(size=12))
dev.off()

##################
#Based on Djamel's code to make cell proportion plots
###Create a proportions table
y=prop.table(x=table(subset_Myo$orig.ident, as.character(subset_Myo$new_cluster_ids)), margin=2)

###Adding proportions for all of the cells
nbr_row <- length(unique(subset_Myo$orig.ident))
ALL <- vector("numeric", nbr_row) # prepare a container

for (i in 1:nbr_row) {
  res <- table(subset_Myo$orig.ident)[[i]]/dim(subset_Myo)[2]
  ALL[i] <- res         # change to assignment
}
ALL

y1 <- cbind(y, All=ALL)
y1

#Change to long version of the table and rename the columns
df1 <- melt(y1)
head(df1)
length(unique(subset_Myo$new_cluster_ids))
colnames(df1) <- c("orig.ident", "Clusters", "Percentage")
head(df1)

#To order the clusters
ordered_clusters <- c('All', 'C8', 'C9', 'C10', 'C11')

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
  geom_hline(yintercept = 0.2168465,linetype="dotted", size=1) 

print(p_Age)

ggsave("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/Myo_Prop.pdf", 
       p_Age , width=3, height=1.7,  units="in", scale=3)

###########################
###Plotting uMAP and Feature Plots

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPS/Myo/Myo_UMAP.pdf", height=20, width=20)
DimPlot(subset_Myo, reduction = "umap", cols = c("#BD9ED8", "#F8B4C0", "#46206B", "#BC46B6"), 
        label = TRUE, label.size = 7) + theme(text = element_text(size=15)) + ylim(c(-3,4.2)) + xlim(c(4,10.5))  
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPS/Myo/Myo_Krt15.pdf", height=20, width=20)
FeaturePlot(subset_Myo, features = c("Krt15"), cols = c("lightgrey", "red3"), label = F, ncol = 1, slot = "scale.data") + ylim(c(-3,4.2)) + xlim(c(4,10.5)) 
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPS/Myo/Myo_Actg2.pdf", height=20, width=20)
FeaturePlot(subset_Myo, features = c("Actg2"), cols = c("lightgrey", "red3"), label = F, ncol = 1, slot = "scale.data") + ylim(c(-3,4.2)) + xlim(c(4,10.5)) 
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/uMAPS/Myo/Myo_Cd74.pdf", height=20, width=20)
FeaturePlot(subset_Myo, features = c("Cd74"), cols = c("lightgrey", "red3"), label = F, ncol = 1, slot = "scale.data") + ylim(c(-3,4.2)) + xlim(c(4,10.5)) 
dev.off()

#############
#> sessionInfo()
#R version 4.1.2 (2021-11-01)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 19042)

#Matrix products: default

#locale:
#  [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                           LC_TIME=English_United States.1252    

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggpubr_0.4.0       reshape2_1.4.4     ggplot2_3.3.5      harmony_0.1.0      Rcpp_1.0.7         cowplot_1.1.1      SeuratObject_4.0.4 Seurat_4.0.5       dplyr_1.0.7       

#loaded via a namespace (and not attached):
#  [1] nlme_3.1-153          matrixStats_0.61.0    spatstat.sparse_2.0-0 RcppAnnoy_0.0.19      RColorBrewer_1.1-2    httr_1.4.2            backports_1.4.0       sctransform_0.3.2     tools_4.1.2           utf8_1.2.2           
#[11] R6_2.5.1              irlba_2.3.3           rpart_4.1-15          KernSmooth_2.23-20    uwot_0.1.10           mgcv_1.8-38           DBI_1.1.1             lazyeval_0.2.2        colorspace_2.0-2      withr_2.4.3          
#[21] tidyselect_1.1.1      gridExtra_2.3         compiler_4.1.2        plotly_4.10.0         labeling_0.4.2        scales_1.1.1          lmtest_0.9-39         spatstat.data_2.1-0   ggridges_0.5.3        pbapply_1.5-0        
#[31] goftest_1.2-3         stringr_1.4.0         digest_0.6.29         spatstat.utils_2.2-0  pkgconfig_2.0.3       htmltools_0.5.2       parallelly_1.29.0     fastmap_1.1.0         htmlwidgets_1.5.4     rlang_0.4.12         
#[41] rstudioapi_0.13       shiny_1.7.1           farver_2.1.0          generics_0.1.1        zoo_1.8-9             jsonlite_1.7.2        ica_1.0-2             car_3.0-12            magrittr_2.0.1        patchwork_1.1.1      
#[51] Matrix_1.3-4          munsell_0.5.0         fansi_0.5.0           abind_1.4-5           reticulate_1.22       lifecycle_1.0.1       stringi_1.7.6         carData_3.0-4         MASS_7.3-54           Rtsne_0.15           
#[61] plyr_1.8.6            grid_4.1.2            parallel_4.1.2        listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1         crayon_1.4.2          miniUI_0.1.1.1        deldir_1.0-6          lattice_0.20-45      
#[71] splines_4.1.2         tensor_1.5            pillar_1.6.4          igraph_1.2.9          spatstat.geom_2.3-0   ggsignif_0.6.3        future.apply_1.8.1    codetools_0.2-18      leiden_0.3.9          glue_1.5.1           
#[81] data.table_1.14.2     png_0.1-7             vctrs_0.3.8           httpuv_1.6.3          gtable_0.3.0          RANN_2.6.1            purrr_0.3.4           spatstat.core_2.3-2   polyclip_1.10-0       tidyr_1.1.4          
#[91] scattermore_0.7       future_1.23.0         assertthat_0.2.1      mime_0.12             broom_0.7.10          xtable_1.8-4          RSpectra_0.16-0       rstatix_0.7.0         later_1.3.0           survival_3.2-13      
#[101] viridisLite_0.4.0     tibble_3.1.6          cluster_2.1.2         globals_0.14.0        fitdistrplus_1.1-6    ellipsis_0.3.2        ROCR_1.0-11  
