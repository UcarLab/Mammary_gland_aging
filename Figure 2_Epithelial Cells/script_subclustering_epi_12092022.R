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

#saveRDS(scRNA_seq_expt_comb_harmony_umap_2, file = "C:/Users/angarb/Box/Anczukow_lab_NGS_raw_data/Aging single cells/Subclustering/Single_Cell_June2021/scRNA_seq_expt_comb_harmony_umap_second_epi_subset_12092022_BA_PC14Res0.8.rds")

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
