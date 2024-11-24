library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(JASPAR2022)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)

##scRNAseq
pbmc.data <- Read10X(data.dir = "~/path/to/filtered_feature_bc_directory")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
##Subset data based on QC
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
## Before running the next steps, please remove the doublets identified from scrublet and DoubletDecon. Merge each object and proceed to normalization.

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

##For subset analyses, change all parameters accordingly and run from PCA step onwards.
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

##Batch correction step is only needed for RNAseq data. Not the ATAC. 
pbmc <- RunHarmony(pbmc,"stim", plot_convergence = TRUE)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "~/path/to/output")

##scATACseq
counts <- Read10X_h5(filename = "~/path/to/filtered_peak_bc_matrix.h5")
metadata <- read.csv(file = "~/path/to/singlecell.csv",header = TRUE,row.names = 1)
chrom_assay <- CreateChromatinAssay(counts = counts,sep = c(":", "-"),fragments = "~/path/to/fragments.tsv.gz",min.cells = 10,min.features = 200)
pbmc_atac <- CreateSeuratObject(counts = chrom_assay,assay = "peaks",meta.data = metadata)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
Annotation(pbmc_atac) <- annotations
pbmc <- NucleosomeSignal(object = pbmc_atac)
pbmc <- TSSEnrichment(object = pbmc_atac)
pbmc_atac$pct_reads_in_peaks <- pbmc_atac$peak_region_fragments / pbmc_atac$passed_filters * 100
pbmc_atac$blacklist_ratio <- FractionCountsInRegion(object = pbmc_atac, assay = 'peaks',regions = blacklist_mm10_regions)
pbmc_atac <- subset(x = pbmc_atac,subset = peak_region_fragments > 1000 & peak_region_fragments < 100000 & pct_reads_in_peaks > 40 & blacklist_ratio < 0.01 & nucleosome_signal < 4 & TSS.enrichment > 2)
## Also exclude all doublets identified using Amulet.
## For merging the ATAC objects, please refer: https://stuartlab.org/signac/articles/merging. Do merging and proceed to TFIDF step. 
pbmc_atac <- RunTFIDF(pbmc_atac)
pbmc_atac <- FindTopFeatures(pbmc_atac, min.cutoff = 'q0')
##For subset analyses, change all parameters accordingly and run from SVD step onwards.
pbmc_atac <- RunSVD(object = pbmc_atac)
pbmc_atac <- RunUMAP(object = pbmc_atac,reduction = 'lsi',dims = 2:30)
pbmc_atac <- FindNeighbors(object = pbmc_atac,reduction = 'lsi',dims = 2:30)
pbmc_atac <- FindClusters(object = pbmc_atac,algorithm = 3,resolution = 1.2,verbose = FALSE)
gene.activities <- GeneActivity(pbmc_atac)
pbmc_atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc_atac <- NormalizeData(object = pbmc_atac,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = median(pbmc_atac$nCount_RNA))
pbmc <- readRDS("allen_pbmc_atac.rds")
pbmc <- UpdateSeuratObject(pbmc)
pbmc <- FindVariableFeatures(object = pbmc,nfeatures = 5000)
transfer.anchors <- FindTransferAnchors(reference = pbmc,query = pbmc_atac,reduction = 'cca',dims = 1:30)
predicted.labels <- TransferData(anchorset = transfer.anchors,refdata = pbmc$subclass,weight.reduction = pbmc_atac[['lsi']],dims = 2:30)
pbmc_atac <- AddMetaData(object = pbmc_atac, metadata = predicted.labels)
pfm <- getMatrixSet(x = JASPAR2020,opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
pbmc_atac <- AddMotifs(object = pbmc_atac,genome = BSgenome.Hsapiens.UCSC.hg38,pfm = pfm)
pbmc_atac <- RunChromVAR(object = pbmc_atac,genome = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(pbmc_atac, file = "~/path/to/output")


##For single cell differential analyses
## Accordingly change assay for pbmc_atac to run differential analyses between motifs and peaks. Keep harmony assay for pbmc. 
## Change test type, and percentages of cells based on methods. This is simply a template. To replicate age wise and celltype wise differentials, please change identity and groups accordingly.
Idents(pbmc) <- pbmc$Celltype
DifferentialData <- FindMarkers(pbmc, ident.1 = "Old", ident.2 = "Young", group.by = "Age", subset.ident = "Celltype")
write.table(DifferentialData,"~/path/to/output")



##For pseudobulk RNAseq differential analyses: https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html


