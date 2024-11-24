
#Script to run pseudobulk DE (differential expression) analysis.
# Load libraries
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

file_obj_orig <- readRDS("../RDS_files/seuratObject_rename_ident_May26_2021.rds")

file_obj <- subset(file_obj_orig, idents = c("Doublet"), invert = TRUE )
#file_obj <- subset(file_obj1, idents = c("Pericytes", "Plasma"), invert = TRUE )
# Bring in Seurat object
# Extract raw counts and metadata to create SingleCellExperiment object
counts <- file_obj@assays$RNA@counts 
metadata <- file_obj@meta.data

# Set up metadata as desired for aggregation and DE analysis
clust_tmp <- sapply(gsub("_", "-", file_obj@active.ident),`[`,1)
clust <- factor(clust_tmp)
metadata$cluster_id <- clust

sample_tmp <- sapply(gsub("_rep",":rep",metadata$replicate),`[`,1)
sample <- factor(sample_tmp)
metadata$sample_id <- sample

metadata$group_id <- factor(file_obj@meta.data$orig.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                           colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Explore the raw counts for the dataset
## Check the assays present
assays(sce)

## Explore the raw counts for the dataset
dim(counts(sce))
counts(sce)[1:6, 1:6]

## Explore the cellular metadata for the dataset
dim(colData(sce))
head(colData(sce))

# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids

# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(sce$sample_id))

# Total number of samples 
ns <- length(sids)
ns

# Generate sample level metadata
## Determine the number of cells per sample
table(sce$sample_id)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
                select(-"cluster_id")
ei

# Perform QC if not already performed
dim(sce)

# Calculate quality control (QC) metrics
#sce <- calculateQCMetrics(sce)

sce <-addPerCellQC(sce)

# Get cells w/ few/many detected genes
#sce$is_outlier <- isOutlier(
#        metric = sce$total_features_by_counts,
#        nmads = 2, type = "both", log = TRUE)

# Remove outlier cells
#sce <- sce[, !sce$is_outlier]
#dim(sce)

## Remove lowly expressed genes which have less than 10 cells with any counts
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

#write.table(sce, "sce_pseudobulk_filtered_low_cells.txt", sep="\t", quote=FALSE)

# Aggregate the counts per sample_id and cluster_id
# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)
dim(pb)
pb[1:6, 1:6]

# Not every cluster is present in all samples; create a vector that represents how to split samples
rownames(pb) <- sapply(gsub("_M",":M",rownames(pb)),`[`,1)
rownames(pb) <- sapply(gsub("_rep","-rep",rownames(pb)),`[`,1)

splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = ":",  
                                    n = 2), 
                 `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
        lapply(function(u) 
                set_colnames(t(u), 
                             stringr::str_extract(rownames(u), "(?<=:)[:graph:]+")))

class(pb)

# Explore the different components of list
str(pb)

# Print out the table of cells in each cluster-sample group
options(width = 100)
table(sce$cluster_id, sce$sample_id)

#Differential gene expression with DESeq2
# Get sample names for each of the cell type clusters
# prep. data.frame for plotting
get_sample_ids <- function(x){
        pb[[x]] %>%
                colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
        unlist()

# Get cluster IDs for each of the samples
samples_list <- map(1:length(kids), get_sample_ids)
get_cluster_ids <- function(x){
        rep(names(pb)[x], 
            each = length(samples_list[[x]]))
}
de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
        unlist()

# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id", "expt")]) 

metadata <- gg_df %>%
        dplyr::select(cluster_id, sample_id, group_id, expt) 
        
head(metadata)

clusters <- levels(factor(metadata$cluster_id))
clusters

################################3
dir.create("results")
dir.create("DESeq2")
dir.create("DESeq2/pairwise")

# Function to run DESeq2 and get results for all clusters
## x is index of cluster in clusters vector on which to run function
## A is the sample group to compare
## B is the sample group to compare against (base level)

#fibro_pseudobulk <- read.table("fibro_pseudobulk_DE_genes.txt", sep="\t", header=T)
#fibro_sc <- read.table("fibro_sc_DE_genes.txt", sep="\t", header=T)
#fibro_sc_pseudo_overlap <- read.table("fibro_sc_pseudo_overlap_DE_genes.txt", sep="\t", header=T)

#lumAV_pseudobulk <- read.table("LumAV_pseudobulk_DE_genes.txt", sep="\t", header=T)
#lumAV_sc <- read.table("LumAV_sc_DE_genes.txt", sep="\t", header=T)
#lumAV_sc_pseudo_overlap <- read.table("LumAV_sc_pseudobulk_overlap.txt", sep="\t", header=T)

library(ggplot2)
get_dds_resultsAvsB_new <- function(x, A, B){
        cluster_metadata <- metadata[which(metadata$cluster_id == clusters[x]), ]
        rownames(cluster_metadata) <- cluster_metadata$sample_id
        counts <- pb[[clusters[x]]]
        cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
        
        #all(rownames(cluster_metadata) == colnames(cluster_counts))        
        dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                      colData = cluster_metadata, 
                                      design = ~ expt + group_id)

        # Transform counts for data visualization
        rld <- rlog(dds, blind=TRUE)
        
	 # Plot PCA
        DESeq2::plotPCA(rld, intgroup = "group_id")
        clusters[x] = gsub("/","_", clusters[x])
        #ggsave(paste0(clusters[x], "_specific_PCAplot.png"))
        ggsave(paste0("./results/", clusters[x], "_specific_PCAplot.png"))

	 # Plot PCA with batch effect
        pca_batch_4color <- DESeq2::plotPCA(rld, intgroup = c("group_id", "expt"))
        #state1 <- (c(rep("", 3), rep("Y", 3)), 2))
	#plot_new_batch <- pca_batch_4color + geom_label(size = 2, aes(label = group))
        plot_new_batch <- pca_batch_4color + scale_colour_manual(values=c("violet","magenta","royalblue1","mediumpurple1"))
	clusters[x] = gsub("/","_or_", clusters[x])
        #ggsave(paste0(clusters[x], "_specific_PCAplot.png"))
        ggsave(filename = paste0("./results/", clusters[x], "_specific_PCAplot_batch.png"), plot = plot_new_batch, width=10, height=10)

        # Transform counts for data visualization
        assay(rld) <- limma::removeBatchEffect(assay(rld), rld$expt)

         # Plot PCA - batch effect corrected
        DESeq2::plotPCA(rld, intgroup = "group_id")
        clusters[x] = gsub("/","_or_", clusters[x])
        #ggsave(paste0(clusters[x], "_specific_PCAplot.png"))
        ggsave(paste0("./results/", clusters[x], "_specific_PCAplot_correct.png"), width=10, height=10)
	 
	# Plot PCA - batch effect corrected (4 colors)
        pca_4color <- DESeq2::plotPCA(rld, intgroup = c("group_id", "expt"))
	#plot_new_correct <- pca_4color + geom_label(size = 2, aes(label = group))
	plot_new_correct <- pca_4color + scale_colour_manual(values=c("violet","magenta","royalblue1","mediumpurple1"))
        clusters[x] = gsub("/","_or_", clusters[x])
        #ggsave(paste0(clusters[x], "_specific_PCAplot.png"))
        ggsave(filename = paste0("./results/", clusters[x], "_specific_PCAplot_correct_4colors.png"), plot = plot_new_correct, width=10, height=10)

        # Extract the rlog matrix from the object and compute pairwise correlation values
        rld_mat <- assay(rld)
        rld_cor <- cor(rld_mat)
        
        # Plot heatmap
        #png(paste0(clusters[x], "_specific_heatmap.png"))
        png(paste0("./results/", clusters[x], "_specific_heatmap.png"))
        pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])
        dev.off()
        
        # Run DESeq2 differential expression analysis
        dds <- DESeq(dds)
        
	#select <- fibro_pseudobulk$gene
	#select <- fibro_sc_pseudo_overlap$gene
	#df <- as.data.frame(colData(rld)[,c("group_id","expt")])
	#pdf("Fibro_heatmap_overlap_sc_pseudobulk_genes.pdf", height=30)
	#pheatmap(assay(rld)[rownames(rld) %in% select,], scale="row", cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
	#dev.off()
	
	#select <- lumAV_sc_pseudo_overlap$gene
	#df <- as.data.frame(colData(rld)[,c("group_id","expt")])
	#pdf("LumAV_heatmap_overlap_sc_pseudobulk_genes.pdf", height=30)
	#pheatmap(assay(rld)[rownames(rld) %in% select,], scale="row", cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
	#dev.off()

	#select <- fibro_sc$gene
	#df <- as.data.frame(colData(rld)[,c("group_id","expt")])
	#pdf("Fibro_heatmap_singlecell_genes.pdf", height=50)
	#pheatmap(assay(rld)[rownames(rld) %in% select,], scale="row", cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
	#dev.off()

        # Plot dispersion estimates
        #png(paste0(clusters[x], "_dispersion_plot.png"))
        png(paste0("results/", clusters[x], "_dispersion_plot.png"))
        plotDispEsts(dds)
        dev.off()

        # Output results of Wald test for contrast for A vs B
        contrast <- c("group_id", levels(cluster_metadata$group_id)[A], levels(cluster_metadata$group_id)[B])
        
        # resultsNames(dds)
        res <- results(dds, 
                       contrast = contrast,
                       alpha = 0.05)
        
        res <- lfcShrink(dds, 
                         contrast =  contrast,
                         res=res, type="normal")
        # Set thresholds
        padj_cutoff_small <- 0.1
        padj_cutoff <- 0.05
	log2FoldChange_cutoff <- 0.50
	log2FoldChange_strict_cutoff <- 1.0
        
        # Turn the results object into a tibble for use with tidyverse functions
        res_tbl <- res %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble()
        
        write.csv(res_tbl,
                  #paste0(clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_all_genes.csv"),
                  paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_all_genes.csv"),
                  quote = FALSE, 
                  row.names = FALSE)
        
	#select <- res_tbl$gene
	#df <- as.data.frame(colData(rld)[,c("group_id","expt")])
	#heatmap_filename_celltype <- paste0(clusters[x], "_heatmap_pseudobulk_genes.pdf")
	#pdf(heatmap_filename_celltype, height=30)
	#pheatmap(assay(rld)[rownames(rld) %in% select,], scale="row", cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
	#dev.off()
	
        # Subset the significant results
        sig_res_small <- dplyr::filter(res_tbl, padj < padj_cutoff_small) %>%
                dplyr::arrange(padj)
        
        write.csv(sig_res_small,
                  #paste0(clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_genes.csv"),
                  paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_genes_FDR0.1.csv"),
                  quote = FALSE, 
                  row.names = FALSE)

        # Subset the significant results
        sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
                dplyr::arrange(padj)
        
        write.csv(sig_res,
                  #paste0(clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_genes.csv"),
                  paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_genes.csv"),
                  quote = FALSE, 
                  row.names = FALSE)
        
	##########-------------------------------------------------#############
	
	# Subset the significant results
        sig_res_upreg <- dplyr::filter(res_tbl, padj < padj_cutoff & log2FoldChange > 0) %>%
                dplyr::arrange(padj)
        
        write.csv(sig_res_upreg,
                  #paste0(clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_genes.csv"),
                  paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_genes_upreg.csv"),
                  quote = FALSE, 
                  row.names = FALSE)
       
        # Subset the significant results
        sig_res_downreg <- dplyr::filter(res_tbl, padj < padj_cutoff & log2FoldChange < 0) %>%
                dplyr::arrange(padj)
        
        write.csv(sig_res_downreg,
                  #paste0(clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_genes.csv"),
                  paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_genes_downreg.csv"),
                  quote = FALSE, 
                  row.names = FALSE)
       
   # Subset the significant results
#        sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
 #               dplyr::arrange(padj)

	# Subset the significant results
	sig_res_fc <- dplyr::filter(sig_res, abs(log2FoldChange) > log2FoldChange_cutoff) %>%
        dplyr::arrange(log2FoldChange)

	# Check significant genes output
	#sig_res_fc

	# Write significant results to file
        write.csv(sig_res_fc,
                  #paste0(clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_genes.csv"),
                  paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_fc_fdr_genes.csv"), quote = FALSE, row.names = FALSE)
	
	sig_res_strict_fc <- dplyr::filter(sig_res, abs(log2FoldChange) > log2FoldChange_strict_cutoff) %>%
        dplyr::arrange(log2FoldChange)

	# Check significant genes output
	#sig_res_fc

	# Write significant results to file
        write.csv(sig_res_strict_fc,
                  #paste0(clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_genes.csv"),
                  paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_fc_fdr_strict_genes.csv"), quote = FALSE, row.names = FALSE)
	
        ## ggplot of top genes
        normalized_counts <- counts(dds, normalized = TRUE)
        
	#select <- res_tbl$gene
        #df <- as.data.frame(colData(rld)[,c("group_id","expt")])
        #heatmap_filename_celltype <- paste0(clusters[x], "_heatmap_pseudobulk_genes.pdf")
        #pdf(heatmap_filename_celltype, height=30)
        #pheatmap(assay(rld)[rownames(rld) %in% select,], scale="row", cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
        #dev.off()
	
	
	# Extract normalized counts for only the significant genes
	sig_norm <- data.frame(normalized_counts) %>%
        	rownames_to_column(var = "gene") %>%
        	dplyr::filter(gene %in% sig_res$gene)
        
	rownames(sig_norm) <- sig_norm[,1]
	df <- as.data.frame(colData(dds)[,c("group_id","expt")])
	rownames(df) <- names(sig_norm[,-1])
	
	print("DE pseudobulk only....")
	print(dim(sig_norm))
	
	# Set a color palette
        
	file_sig_norm <- paste0("./DESeq2/",clusters[x], "_sig_norm_data.txt")
	write.table(sig_norm, file_sig_norm, sep="\t", quote=FALSE)

	filename_heatmap_row <- paste0("../Heatmaps/Heatmap_", clusters[x], "_pseudobulk_DE_genes_cluster_rows.pdf")
	pdf(filename_heatmap_row, height=30)
        pheatmap(sig_norm[,2:length(colnames(sig_norm))], scale="row", cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
        dev.off()
	
	filename_heatmap_column <- paste0("../Heatmaps/Heatmap_", clusters[x], "_pseudobulk_DE_genes_cluster_cols.pdf")
	pdf(filename_heatmap_column, height=30)
        pheatmap(sig_norm[,2:length(colnames(sig_norm))], scale="row", cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df)
        dev.off()
	
	filename_heatmap_rows_column <- paste0("../Heatmaps/Heatmap_", clusters[x], "_pseudobulk_DE_genes_cluster_rows_and_cols.pdf")
	pdf(filename_heatmap_rows_column, height=30)
        pheatmap(sig_norm[,2:length(colnames(sig_norm))], scale="row", cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df)
        dev.off()

	##############--------------------------------------------###########################
	#select <- sig_res$gene
	#select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
	#df <- as.data.frame(colData(dds)[,c("group_id","expt")])
	#pdf("Heatmap_all_DE_genes_cluster_rows.pdf")
	#pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
	#dev.off()
        ## Order results by padj values
        top20_sig_genes <- sig_res %>%
                dplyr::arrange(padj) %>%
                dplyr::pull(gene) %>%
                head(n=20)
        
        top20_sig_norm <- data.frame(normalized_counts) %>%
                rownames_to_column(var = "gene") %>%
                dplyr::filter(gene %in% top20_sig_genes)
        
        gathered_top20_sig <- top20_sig_norm %>%
                gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")
     
        gathered_top20_sig$samplename <- gsub("\\.", ":", gathered_top20_sig$samplename) 
        gathered_top20_sig <- inner_join(ei[, c("sample_id", "group_id" )], gathered_top20_sig, by = c("sample_id" = "samplename"))
        
        ## plot using ggplot2
        ggplot(gathered_top20_sig) +
                geom_point(aes(x = gene, 
                               y = normalized_counts, 
                               color = group_id), 
                           position=position_jitter(w=0.1,h=0)) +
                scale_y_log10() +
                xlab("Genes") +
                ylab("log10 Normalized Counts") +
                ggtitle("Top 20 Significant DE Genes") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
                theme(plot.title = element_text(hjust = 0.5))
        ggsave(paste0("./results/",clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_top20_DE_genes.png"))

	#Generate heatmap.
	#select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
	#df <- as.data.frame(colData(dds)[,c("condition","type")])
	#heatmap_file <- paste0("Heatmap_", celltype, ".pdf")
	#pdf(heatmap_file)
	#pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
	#dev.off()
}

# Run the script on all clusters comparing stim condition relative to control condition
map(1:length(clusters), get_dds_resultsAvsB_new, A = 1, B = 2)

