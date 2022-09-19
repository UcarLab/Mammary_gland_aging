#Script for analysis of scATAC-seq merged with QC.
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
set.seed(1234)
library(GenomicRanges)
library(future)
library(cowplot)

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

peaks.mm_3m_rep4 <- read.table(
  file = "/OA20012/cellranger/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.mm_3m_rep5 <- read.table(
  file = "/OA20016/cellranger/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.mm_3m_rep6 <- read.table(
  file = "/OA20020/cellranger/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.mm_18m_rep4 <- read.table(
  file = "/OA20014/cellranger/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.mm_18m_rep5 <- read.table(
  file = "/OA20018/cellranger/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.mm_18m_rep6 <- read.table(
  file = "/OA20022/cellranger/peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.mm_3m_rep4 <- makeGRangesFromDataFrame(peaks.mm_3m_rep4)
gr.mm_3m_rep5 <- makeGRangesFromDataFrame(peaks.mm_3m_rep5)
gr.mm_3m_rep6 <- makeGRangesFromDataFrame(peaks.mm_3m_rep6)
gr.mm_18m_rep4 <- makeGRangesFromDataFrame(peaks.mm_18m_rep4)
gr.mm_18m_rep5 <- makeGRangesFromDataFrame(peaks.mm_18m_rep5)
gr.mm_18m_rep6 <- makeGRangesFromDataFrame(peaks.mm_18m_rep6)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.mm_3m_rep4, gr.mm_3m_rep5, gr.mm_3m_rep6, gr.mm_18m_rep4, gr.mm_18m_rep5, gr.mm_18m_rep6))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

#-------------------------------------------------------#
# load metadata
#Create Fragment Objects.
md.mm_3m_rep4 <- read.table(
  file = "/OA20012/cellranger/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
md.mm_3m_rep5 <- read.table(
  file = "/OA20016/cellranger/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.mm_3m_rep6 <- read.table(
  file = "/OA20020/cellranger/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.mm_18m_rep4 <- read.table(
  file = "/OA20014/cellranger/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.mm_18m_rep5 <- read.table(
  file = "/OA20018/cellranger/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.mm_18m_rep6 <- read.table(
  file = "/OA20022/cellranger/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

# perform an initial filtering of low count cells
md.mm_3m_rep4 <- md.mm_3m_rep4[md.mm_3m_rep4$passed_filters > 500, ]
md.mm_3m_rep5 <- md.mm_3m_rep5[md.mm_3m_rep5$passed_filters > 500, ]
md.mm_3m_rep6 <- md.mm_3m_rep6[md.mm_3m_rep6$passed_filters > 500, ]
md.mm_18m_rep4 <- md.mm_18m_rep4[md.mm_18m_rep4$passed_filters > 500, ]
md.mm_18m_rep5 <- md.mm_18m_rep5[md.mm_18m_rep5$passed_filters > 500, ]
md.mm_18m_rep6 <- md.mm_18m_rep6[md.mm_18m_rep6$passed_filters > 500, ] # sequenced deeper so set higher cutoff

# create fragment objects
frags.3m_rep4 <- CreateFragmentObject(
  path = "/OA20012/cellranger/fragments.tsv.gz",
  cells = rownames(md.mm_3m_rep4)
)
frags.3m_rep5 <- CreateFragmentObject(
  path = "/OA20016/cellranger/fragments.tsv.gz",
  cells = rownames(md.mm_3m_rep5)
)
## Computing hash
frags.3m_rep6 <- CreateFragmentObject(
  path = "/OA20020/cellranger/fragments.tsv.gz",
  cells = rownames(md.mm_3m_rep6)
)
## Computing hash
frags.18m_rep4 <- CreateFragmentObject(
  path = "/OA20014/cellranger/fragments.tsv.gz",
  cells = rownames(md.mm_18m_rep4)
)
frags.18m_rep5 <- CreateFragmentObject(
  path = "/OA20018/cellranger/fragments.tsv.gz",
  cells = rownames(md.mm_18m_rep5)
)
## Computing hash
frags.18m_rep6 <- CreateFragmentObject(
  path = "/OA20022/cellranger/fragments.tsv.gz",
  cells = rownames(md.mm_18m_rep6)
)

#Quantify peaks in each dataset
#We can now create a matrix of peaks x cell for each sample using the FeatureMatrix function. This function is parallelized using the future package. See the parallelization vignette for more information about using future.

mm10_3m_rep4.counts <- FeatureMatrix(
  fragments = frags.3m_rep4,
  features = combined.peaks,
  cells = rownames(md.mm_3m_rep4)
)
mm10_3m_rep5.counts <- FeatureMatrix(
  fragments = frags.3m_rep5,
  features = combined.peaks,
  cells = rownames(md.mm_3m_rep5)
)

mm10_3m_rep6.counts <- FeatureMatrix(
  fragments = frags.3m_rep6,
  features = combined.peaks,
  cells = rownames(md.mm_3m_rep6)
)

mm10_18m_rep4.counts <- FeatureMatrix(
  fragments = frags.18m_rep4,
  features = combined.peaks,
  cells = rownames(md.mm_18m_rep4)
)
mm10_18m_rep5.counts <- FeatureMatrix(
  fragments = frags.18m_rep5,
  features = combined.peaks,
  cells = rownames(md.mm_18m_rep5)
)

mm10_18m_rep6.counts <- FeatureMatrix(
  fragments = frags.18m_rep6,
  features = combined.peaks,
  cells = rownames(md.mm_18m_rep6)
)

##----------Metadata.............###################333
metadata_3m_rep4 <- read.csv(
  file = "/OA20012/cellranger/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_3m_rep5 <- read.csv(
  file = "/OA20016/cellranger/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_3m_rep6 <- read.csv(
  file = "/OA20020/cellranger/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_18m_rep4 <- read.csv(
  file = "/OA20014/cellranger/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_18m_rep5 <- read.csv(
  file = "/OA20018/cellranger/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_18m_rep6 <- read.csv(
  file = "/OA20022/cellranger/singlecell.csv",
  header = TRUE,
  row.names = 1
)

#Create the objects
#We will now use the quantified matrices to create a Seurat object for each dataset, storing the Fragment object for each dataset in the assay.
#mm10_3m_rep4_assay <- CreateChromatinAssay(mm10_3m_rep4.counts, sep=c(":", "-"), fragments=frags.3m_rep4, genome='mm10', min.cells=10, min.features=200)
mm10_3m_rep4_assay <- CreateChromatinAssay(mm10_3m_rep4.counts, sep=c(":", "-"), fragments=frags.3m_rep4, genome='mm10')
mm10_3m_rep4 <- CreateSeuratObject(mm10_3m_rep4_assay, assay = "ATAC", meta.data = metadata_3m_rep4)

#mm10_3m_rep5_assay <- CreateChromatinAssay(mm10_3m_rep5.counts, sep=c(":", "-"), fragments=frags.3m_rep5, genome='mm10', min.cells=10, min.features=200)
mm10_3m_rep5_assay <- CreateChromatinAssay(mm10_3m_rep5.counts, sep=c(":", "-"), fragments=frags.3m_rep5, genome='mm10')
mm10_3m_rep5 <- CreateSeuratObject(mm10_3m_rep5_assay, assay = "ATAC", meta.data = metadata_3m_rep5)

#mm10_3m_rep6_assay <- CreateChromatinAssay(mm10_3m_rep6.counts, sep=c(":", "-"), fragments=frags.3m_rep6, genome='mm10', min.cells=10, min.features=200)
mm10_3m_rep6_assay <- CreateChromatinAssay(mm10_3m_rep6.counts, sep=c(":", "-"), fragments=frags.3m_rep6, genome='mm10')
mm10_3m_rep6 <- CreateSeuratObject(mm10_3m_rep6_assay, assay = "ATAC", meta.data = metadata_3m_rep6)

#mm10_18m_rep4_assay <- CreateChromatinAssay(mm10_18m_rep4.counts, sep=c(":", "-"), fragments=frags.18m_rep4, genome='mm10', min.cells=10, min.features=200)
mm10_18m_rep4_assay <- CreateChromatinAssay(mm10_18m_rep4.counts, sep=c(":", "-"), fragments=frags.18m_rep4, genome='mm10')
mm10_18m_rep4 <- CreateSeuratObject(mm10_18m_rep4_assay, assay = "ATAC", meta.data = metadata_18m_rep4)

#mm10_18m_rep5_assay <- CreateChromatinAssay(mm10_18m_rep5.counts, sep=c(":", "-"), fragments=frags.18m_rep5, genome='mm10', min.cells=10, min.features=200)
mm10_18m_rep5_assay <- CreateChromatinAssay(mm10_18m_rep5.counts, sep=c(":", "-"), fragments=frags.18m_rep5, genome='mm10')
mm10_18m_rep5 <- CreateSeuratObject(mm10_18m_rep5_assay, assay = "ATAC", meta.data = metadata_18m_rep5)

#mm10_18m_rep6_assay <- CreateChromatinAssay(mm10_18m_rep6.counts, sep = c(":", "-"), fragments = frags.18m_rep6, genome='mm10', min.cells=10, min.features=200)
mm10_18m_rep6_assay <- CreateChromatinAssay(mm10_18m_rep6.counts, sep = c(":", "-"), fragments = frags.18m_rep6, genome='mm10')
mm10_18m_rep6 <- CreateSeuratObject(mm10_18m_rep6_assay, assay = "ATAC", meta.data = metadata_18m_rep6)

#Merge objects
#Now that the objects each contain an assay with the same set of features, we can use the standard merge function to merge the objects. This will also merge all the fragment objects so that we retain the fragment information for each cell in the final merged object.
# add information to identify dataset of origin
mm10_3m_rep4$dataset <- 'mm10_3m_rep4'
mm10_3m_rep5$dataset <- 'mm10_3m_rep5'
mm10_3m_rep6$dataset <- 'mm10_3m_rep6'
mm10_18m_rep4$dataset <- 'mm10_18m_rep4'
mm10_18m_rep5$dataset <- 'mm10_18m_rep5'
mm10_18m_rep6$dataset <- 'mm10_18m_rep6'

mm10_3m_rep4$age <- '3M'
mm10_3m_rep5$age <- '3M'
mm10_3m_rep6$age <- '3M'
mm10_18m_rep4$age <- '18M'
mm10_18m_rep5$age <- '18M'
mm10_18m_rep6$age <- '18M'

mm10_3m_rep4$rep <- 'rep4'
mm10_3m_rep5$rep <- 'rep5'
mm10_3m_rep6$rep <- 'rep6'
mm10_18m_rep4$rep <- 'rep4'
mm10_18m_rep5$rep <- 'rep5'
mm10_18m_rep6$rep <- 'rep6'

############################################################
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to mm10
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(mm10_3m_rep4) <- annotations
Annotation(mm10_3m_rep5) <- annotations
Annotation(mm10_3m_rep6) <- annotations
Annotation(mm10_18m_rep4) <- annotations
Annotation(mm10_18m_rep5) <- annotations
Annotation(mm10_18m_rep6) <- annotations

#Computing QC Metrics
# compute nucleosome signal score per cell
mm10_3m_rep4 <- NucleosomeSignal(object = mm10_3m_rep4)
mm10_3m_rep5 <- NucleosomeSignal(object = mm10_3m_rep5)
mm10_3m_rep6 <- NucleosomeSignal(object = mm10_3m_rep6)
mm10_18m_rep4 <- NucleosomeSignal(object = mm10_18m_rep4)
mm10_18m_rep5 <- NucleosomeSignal(object = mm10_18m_rep5)
mm10_18m_rep6 <- NucleosomeSignal(object = mm10_18m_rep6)

mm10_3m_rep4$nucleosome_group <- ifelse(mm10_3m_rep4$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
mm10_3m_rep5$nucleosome_group <- ifelse(mm10_3m_rep5$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
mm10_3m_rep6$nucleosome_group <- ifelse(mm10_3m_rep6$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
mm10_18m_rep4$nucleosome_group <- ifelse(mm10_18m_rep4$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
mm10_18m_rep5$nucleosome_group <- ifelse(mm10_18m_rep5$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
mm10_18m_rep6$nucleosome_group <- ifelse(mm10_18m_rep6$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# compute TSS enrichment score per cell
mm10_3m_rep4 <- TSSEnrichment(object = mm10_3m_rep4, fast = FALSE)
mm10_3m_rep5 <- TSSEnrichment(object = mm10_3m_rep5, fast = FALSE)
mm10_3m_rep6 <- TSSEnrichment(object = mm10_3m_rep6, fast = FALSE)
mm10_18m_rep4 <- TSSEnrichment(object = mm10_18m_rep4, fast = FALSE)
mm10_18m_rep5 <- TSSEnrichment(object = mm10_18m_rep5, fast = FALSE)
mm10_18m_rep6 <- TSSEnrichment(object = mm10_18m_rep6, fast = FALSE)

saveRDS(mm10_3m_rep4, "3m_rep4_TSS_enrich_orig.rds")
saveRDS(mm10_3m_rep5, "3m_rep5_TSS_enrich_orig.rds")
saveRDS(mm10_3m_rep6, "3m_rep6_TSS_enrich_orig.rds")
saveRDS(mm10_18m_rep4, "18m_rep4_TSS_enrich_orig.rds")
saveRDS(mm10_18m_rep5, "18m_rep5_TSS_enrich_orig.rds")
saveRDS(mm10_18m_rep6, "18m_rep6_TSS_enrich_orig.rds")

####################################################################

# add blacklist ratio and fraction of reads in peaks
mm10_3m_rep4$pct_reads_in_peaks <- mm10_3m_rep4$peak_region_fragments / mm10_3m_rep4$passed_filters * 100
mm10_3m_rep4$blacklist_ratio <- mm10_3m_rep4$blacklist_region_fragments / mm10_3m_rep4$peak_region_fragments
mm10_3m_rep4$high.tss <- ifelse(mm10_3m_rep4$TSS.enrichment > 2, 'High', 'Low')

mm10_3m_rep5$pct_reads_in_peaks <- mm10_3m_rep5$peak_region_fragments / mm10_3m_rep5$passed_filters * 100
mm10_3m_rep5$blacklist_ratio <- mm10_3m_rep5$blacklist_region_fragments / mm10_3m_rep5$peak_region_fragments
mm10_3m_rep5$high.tss <- ifelse(mm10_3m_rep5$TSS.enrichment > 2, 'High', 'Low')

mm10_3m_rep6$pct_reads_in_peaks <- mm10_3m_rep6$peak_region_fragments / mm10_3m_rep6$passed_filters * 100
mm10_3m_rep6$blacklist_ratio <- mm10_3m_rep6$blacklist_region_fragments / mm10_3m_rep6$peak_region_fragments
mm10_3m_rep6$high.tss <- ifelse(mm10_3m_rep6$TSS.enrichment > 2, 'High', 'Low')

mm10_18m_rep4$pct_reads_in_peaks <- mm10_18m_rep4$peak_region_fragments / mm10_18m_rep4$passed_filters * 100
mm10_18m_rep4$blacklist_ratio <- mm10_18m_rep4$blacklist_region_fragments / mm10_18m_rep4$peak_region_fragments
mm10_18m_rep4$high.tss <- ifelse(mm10_18m_rep4$TSS.enrichment > 2, 'High', 'Low')

mm10_18m_rep5$pct_reads_in_peaks <- mm10_18m_rep5$peak_region_fragments / mm10_18m_rep5$passed_filters * 100
mm10_18m_rep5$blacklist_ratio <- mm10_18m_rep5$blacklist_region_fragments / mm10_18m_rep5$peak_region_fragments
mm10_18m_rep5$high.tss <- ifelse(mm10_18m_rep5$TSS.enrichment > 2, 'High', 'Low')

mm10_18m_rep6$pct_reads_in_peaks <- mm10_18m_rep6$peak_region_fragments / mm10_18m_rep6$passed_filters * 100
mm10_18m_rep6$blacklist_ratio <- mm10_18m_rep6$blacklist_region_fragments / mm10_18m_rep6$peak_region_fragments
mm10_18m_rep6$high.tss <- ifelse(mm10_18m_rep6$TSS.enrichment > 2, 'High', 'Low')

dir.create("./Plots")
pdf("./Plots/TSSPlot_3m_rep4_orig.pdf")
TSSPlot(mm10_3m_rep4, group.by = 'high.tss') + NoLegend()
dev.off()
pdf("./Plots/TSSPlot_3m_rep5_orig.pdf")
TSSPlot(mm10_3m_rep5, group.by = 'high.tss') + NoLegend()
dev.off()
pdf("./Plots/TSSPlot_3m_rep6_orig.pdf")
TSSPlot(mm10_3m_rep6, group.by = 'high.tss') + NoLegend()
dev.off()

pdf("./Plots/TSSPlot_18m_rep4_orig.pdf")
TSSPlot(mm10_18m_rep4, group.by = 'high.tss') + NoLegend()
dev.off()
pdf("./Plots/TSSPlot_18m_rep5_orig.pdf")
TSSPlot(mm10_18m_rep5, group.by = 'high.tss') + NoLegend()
dev.off()
pdf("./Plots/TSSPlot_18m_rep6_orig.pdf")
TSSPlot(mm10_18m_rep6, group.by = 'high.tss') + NoLegend()
dev.off()

mm10_3m_rep4$nucleosome_group <- ifelse(mm10_3m_rep4$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
mm10_3m_rep5$nucleosome_group <- ifelse(mm10_3m_rep5$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
mm10_3m_rep6$nucleosome_group <- ifelse(mm10_3m_rep6$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
mm10_18m_rep4$nucleosome_group <- ifelse(mm10_18m_rep4$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
mm10_18m_rep5$nucleosome_group <- ifelse(mm10_18m_rep5$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
mm10_18m_rep6$nucleosome_group <- ifelse(mm10_18m_rep6$nucleosome_signal > 10, 'NS > 10', 'NS < 10')

pdf("./Plots/Fragment_hist_3m_rep4_orig.pdf")
FragmentHistogram(object = mm10_3m_rep4, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
dev.off()

pdf("./Plots/Fragment_hist_3m_rep5_orig.pdf")
FragmentHistogram(object = mm10_3m_rep5, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
dev.off()

pdf("./Plots/Fragment_hist_3m_rep6_orig.pdf")
FragmentHistogram(object = mm10_3m_rep6, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
dev.off()

pdf("./Plots/Fragment_hist_18m_rep4_orig.pdf")
FragmentHistogram(object = mm10_18m_rep4, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
dev.off()

pdf("./Plots/Fragment_hist_18m_rep5_orig.pdf")
FragmentHistogram(object = mm10_18m_rep5, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
dev.off()

pdf("./Plots/Fragment_hist_18m_rep6_orig.pdf")
FragmentHistogram(object = mm10_18m_rep6, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
dev.off()

pdf("./Plots/ViolinPlot_mm10_3M_rep4_orig.pdf", width=20)
VlnPlot(
  object = mm10_3m_rep4,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("./Plots/ViolinPlot_mm10_3M_rep5_orig.pdf", width=20)
VlnPlot(
  object = mm10_3m_rep5,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("./Plots/ViolinPlot_mm10_3M_rep6_orig.pdf", width=20)
VlnPlot(
  object = mm10_3m_rep6,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("./Plots/ViolinPlot_mm10_18M_rep4_orig.pdf", width=20)
VlnPlot(
  object = mm10_18m_rep4,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("./Plots/ViolinPlot_mm10_18M_rep5_orig.pdf", width=20)
VlnPlot(
  object = mm10_18m_rep5,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("./Plots/ViolinPlot_mm10_18M_rep6_orig.pdf", width=20)
VlnPlot(
  object = mm10_18m_rep6,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

########################################################################
mm10_3m_rep4 <- subset(
  x = mm10_3m_rep4,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
mm10_3m_rep4

mm10_3m_rep5 <- subset(
  x = mm10_3m_rep5,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
mm10_3m_rep5

mm10_3m_rep6 <- subset(
  x = mm10_3m_rep6,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
mm10_3m_rep6

mm10_18m_rep4 <- subset(
  x = mm10_18m_rep4,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
mm10_18m_rep4

mm10_18m_rep5 <- subset(
  x = mm10_18m_rep5,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
mm10_18m_rep5

mm10_18m_rep6 <- subset(
  x = mm10_18m_rep6,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
mm10_18m_rep6

##########################################

saveRDS(mm10_3m_rep4, "mm10_3M_rep4.rds")
saveRDS(mm10_3m_rep5, "mm10_3M_rep5.rds")
saveRDS(mm10_3m_rep6, "mm10_3M_rep6.rds")
saveRDS(mm10_18m_rep4, "mm10_18M_rep4.rds")
saveRDS(mm10_18m_rep5, "mm10_18M_rep5.rds")
saveRDS(mm10_18m_rep6, "mm10_18M_rep6.rds")

###########################################
mm10_3m_rep4 <- readRDS("../RDS_files/mm10_3M_rep4.rds")
mm10_3m_rep5 <- readRDS("../RDS_files/mm10_3M_rep5.rds")
mm10_3m_rep6 <- readRDS("../RDS_files/mm10_3M_rep6.rds")
mm10_18m_rep4 <- readRDS("../RDS_files/mm10_18M_rep4.rds")
mm10_18m_rep5 <- readRDS("../RDS_files/mm10_18M_rep5.rds")
mm10_18m_rep6 <- readRDS("../RDS_files/mm10_18M_rep6.rds")

##########################################

pdf("../Plots/ViolinPlot_mm10_3M_rep4_after_filter.pdf", width=20)
VlnPlot(
  object = mm10_3m_rep4,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("../Plots/ViolinPlot_mm10_3M_rep5_after_filter.pdf", width=20)
VlnPlot(
  object = mm10_3m_rep5,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("../Plots/ViolinPlot_mm10_3M_rep6_after_filter.pdf", width=20)
VlnPlot(
  object = mm10_3m_rep6,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("../Plots/ViolinPlot_mm10_18M_rep4_after_filter.pdf", width=20)
VlnPlot(
  object = mm10_18m_rep4,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("../Plots/ViolinPlot_mm10_18M_rep5_after_filter.pdf", width=20)
VlnPlot(
  object = mm10_18m_rep5,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("../Plots/ViolinPlot_mm10_18M_rep6_after_filter.pdf", width=20)
VlnPlot(
  object = mm10_18m_rep6,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

###########################################
mm10_3m_rep4_barcode_list <- names(mm10_3m_rep4@active.ident)
M3_rep4_singlecell_subset <- subset(metadata_3m_rep4, rownames(metadata_3m_rep4) %in% mm10_3m_rep4_barcode_list)
write.table(M3_rep4_singlecell_subset, "M3_rep4_singlecell.csv", sep=",", quote=FALSE)

mm10_3m_rep5_barcode_list <- names(mm10_3m_rep5@active.ident)
M3_rep5_singlecell_subset <- subset(metadata_3m_rep5, rownames(metadata_3m_rep5) %in% mm10_3m_rep5_barcode_list)
write.table(M3_rep5_singlecell_subset, "M3_rep5_singlecell.csv", sep=",", quote=FALSE)

mm10_3m_rep6_barcode_list <- names(mm10_3m_rep6@active.ident)
M3_rep6_singlecell_subset <- subset(metadata_3m_rep6, rownames(metadata_3m_rep6) %in% mm10_3m_rep6_barcode_list)
write.table(M3_rep6_singlecell_subset, "M3_rep6_singlecell.csv", sep=",", quote=FALSE)

mm10_18m_rep4_barcode_list <- names(mm10_18m_rep4@active.ident)
M18_rep4_singlecell_subset <- subset(metadata_18m_rep4, rownames(metadata_18m_rep4) %in% mm10_18m_rep4_barcode_list)
write.table(M18_rep4_singlecell_subset, "M18_rep4_singlecell.csv", sep=",", quote=FALSE)

mm10_18m_rep5_barcode_list <- names(mm10_18m_rep5@active.ident)
M18_rep5_singlecell_subset <- subset(metadata_18m_rep5, rownames(metadata_18m_rep5) %in% mm10_18m_rep5_barcode_list)
write.table(M18_rep5_singlecell_subset, "M18_rep5_singlecell.csv", sep=",", quote=FALSE)

mm10_18m_rep6_barcode_list <- names(mm10_18m_rep6@active.ident)
M18_rep6_singlecell_subset <- subset(metadata_18m_rep6, rownames(metadata_18m_rep6) %in% mm10_18m_rep6_barcode_list)
write.table(M18_rep6_singlecell_subset, "M18_rep6_singlecell.csv", sep=",", quote=FALSE)

# extract gene annotations from EnsDb
#Remove doublets

doublets_3m_rep4 = read.table("/Output_3M_rep4_new_doublet/DoubletCellIds_01.txt")
doublets_3m_rep5 = read.table("/Output_3M_rep5_new_doublet/DoubletCellIds_01.txt")
doublets_3m_rep6 = read.table("/Output_3M_rep6_new_doublet/DoubletCellIds_01.txt")
doublets_18m_rep4 = read.table("/Output_18M_rep4_new_doublet/DoubletCellIds_01.txt")
doublets_18m_rep5 = read.table("/Output_18M_rep5_new_doublet/DoubletCellIds_01.txt")
doublets_18m_rep6 = read.table("/Output_18M_rep6_new_doublet/DoubletCellIds_01.txt")

mm_3d_rep4 = as.vector(unlist(doublets_3m_rep4$V1))
mm_3d_rep5 = as.vector(unlist(doublets_3m_rep5$V1))
mm_3d_rep6 = as.vector(unlist(doublets_3m_rep6$V1))
mm_18d_rep4 = as.vector(unlist(doublets_18m_rep4$V1))
mm_18d_rep5 = as.vector(unlist(doublets_18m_rep5$V1))
mm_18d_rep6 = as.vector(unlist(doublets_18m_rep6$V1))

mm_3_rep4 = read.table("/OA20012/cellranger/singlecell.csv", sep=",", header=TRUE)
mm_3_rep5 = read.table("/OA20016/cellranger/singlecell.csv", sep=",", header=TRUE)
mm_3_rep6 = read.table("/OA20020/cellranger/singlecell.csv", sep=",", header=TRUE)
mm_18_rep4 = read.table("/OA20014/cellranger/singlecell.csv", sep=",", header=TRUE)
mm_18_rep5 = read.table("/OA20018/cellranger/singlecell.csv", sep=",", header=TRUE)
mm_18_rep6 = read.table("/OA20022/cellranger/singlecell.csv", sep=",", header=TRUE)

library(dplyr)
doub_3m_rep4 = dplyr::filter(mm_3_rep4, cell_id %in% mm_3d_rep4)
nondoublets_3m_rep4 = subset(mm10_3m_rep4, assay = NULL, cells = doub_3m_rep4$barcode, low.threshold = -Inf, high.threshold = Inf, invert=TRUE)

doub_3m_rep5 = dplyr::filter(mm_3_rep5, cell_id %in% mm_3d_rep5)
nondoublets_3m_rep5 = subset(mm10_3m_rep5, assay = NULL, cells = doub_3m_rep5$barcode, low.threshold = -Inf, high.threshold = Inf, invert=TRUE)

doub_3m_rep6 = dplyr::filter(mm_3_rep6, cell_id %in% mm_3d_rep6)
nondoublets_3m_rep6 = subset(mm10_3m_rep6, assay = NULL, cells = doub_3m_rep6$barcode, low.threshold = -Inf, high.threshold = Inf, invert=TRUE)

doub_18m_rep4 = dplyr::filter(mm_18_rep4, cell_id %in% mm_18d_rep4)
nondoublets_18m_rep4 = subset(mm10_18m_rep4, assay = NULL, cells = doub_18m_rep4$barcode, low.threshold = -Inf, high.threshold = Inf, invert=TRUE)

doub_18m_rep5 = dplyr::filter(mm_18_rep5, cell_id %in% mm_18d_rep5)
nondoublets_18m_rep5 = subset(mm10_18m_rep5, assay = NULL, cells = doub_18m_rep5$barcode, low.threshold = -Inf, high.threshold = Inf, invert=TRUE)

doub_18m_rep6 = dplyr::filter(mm_18_rep6, cell_id %in% mm_18d_rep6)
nondoublets_18m_rep6 = subset(mm10_18m_rep6, assay = NULL, cells = doub_18m_rep6$barcode, low.threshold = -Inf, high.threshold = Inf, invert=TRUE)

############################################################
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to mm10
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(nondoublets_3m_rep4) <- annotations
Annotation(nondoublets_3m_rep5) <- annotations
Annotation(nondoublets_3m_rep6) <- annotations
Annotation(nondoublets_18m_rep4) <- annotations
Annotation(nondoublets_18m_rep5) <- annotations
Annotation(nondoublets_18m_rep6) <- annotations

#Computing QC Metrics
# compute nucleosome signal score per cell
nondoublets_3m_rep4 <- NucleosomeSignal(object = nondoublets_3m_rep4)
nondoublets_3m_rep5 <- NucleosomeSignal(object = nondoublets_3m_rep5)
nondoublets_3m_rep6 <- NucleosomeSignal(object = nondoublets_3m_rep6)
nondoublets_18m_rep4 <- NucleosomeSignal(object = nondoublets_18m_rep4)
nondoublets_18m_rep5 <- NucleosomeSignal(object = nondoublets_18m_rep5)
nondoublets_18m_rep6 <- NucleosomeSignal(object = nondoublets_18m_rep6)

nondoublets_3m_rep4$nucleosome_group <- ifelse(nondoublets_3m_rep4$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
nondoublets_3m_rep5$nucleosome_group <- ifelse(nondoublets_3m_rep5$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
nondoublets_3m_rep6$nucleosome_group <- ifelse(nondoublets_3m_rep6$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
nondoublets_18m_rep4$nucleosome_group <- ifelse(nondoublets_18m_rep4$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
nondoublets_18m_rep5$nucleosome_group <- ifelse(nondoublets_18m_rep5$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
nondoublets_18m_rep6$nucleosome_group <- ifelse(nondoublets_18m_rep6$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# compute TSS enrichment score per cell
nondoublets_3m_rep4 <- TSSEnrichment(object = nondoublets_3m_rep4, fast = FALSE)
nondoublets_3m_rep5 <- TSSEnrichment(object = nondoublets_3m_rep5, fast = FALSE)
nondoublets_3m_rep6 <- TSSEnrichment(object = nondoublets_3m_rep6, fast = FALSE)
nondoublets_18m_rep4 <- TSSEnrichment(object = nondoublets_18m_rep4, fast = FALSE)
nondoublets_18m_rep5 <- TSSEnrichment(object = nondoublets_18m_rep5, fast = FALSE)
nondoublets_18m_rep6 <- TSSEnrichment(object = nondoublets_18m_rep6, fast = FALSE)

dir.create("./RDS_files")
saveRDS(nondoublets_3m_rep4, "./RDS_files/3m_rep4_TSS_enrich.rds")
saveRDS(nondoublets_3m_rep5, "./RDS_files/3m_rep5_TSS_enrich.rds")
saveRDS(nondoublets_3m_rep6, "./RDS_files/3m_rep6_TSS_enrich.rds")
saveRDS(nondoublets_18m_rep4, "./RDS_files/18m_rep4_TSS_enrich.rds")
saveRDS(nondoublets_18m_rep5, "./RDS_files/18m_rep5_TSS_enrich.rds")
saveRDS(nondoublets_18m_rep6, "./RDS_files/18m_rep6_TSS_enrich.rds")

####################################################################
# add blacklist ratio and fraction of reads in peaks
nondoublets_3m_rep4$pct_reads_in_peaks <- nondoublets_3m_rep4$peak_region_fragments / nondoublets_3m_rep4$passed_filters * 100
nondoublets_3m_rep4$blacklist_ratio <- nondoublets_3m_rep4$blacklist_region_fragments / nondoublets_3m_rep4$peak_region_fragments
nondoublets_3m_rep4$high.tss <- ifelse(nondoublets_3m_rep4$TSS.enrichment > 2, 'High', 'Low')

nondoublets_3m_rep5$pct_reads_in_peaks <- nondoublets_3m_rep5$peak_region_fragments / nondoublets_3m_rep5$passed_filters * 100
nondoublets_3m_rep5$blacklist_ratio <- nondoublets_3m_rep5$blacklist_region_fragments / nondoublets_3m_rep5$peak_region_fragments
nondoublets_3m_rep5$high.tss <- ifelse(nondoublets_3m_rep5$TSS.enrichment > 2, 'High', 'Low')

nondoublets_3m_rep6$pct_reads_in_peaks <- nondoublets_3m_rep6$peak_region_fragments / nondoublets_3m_rep6$passed_filters * 100
nondoublets_3m_rep6$blacklist_ratio <- nondoublets_3m_rep6$blacklist_region_fragments / nondoublets_3m_rep6$peak_region_fragments
nondoublets_3m_rep6$high.tss <- ifelse(nondoublets_3m_rep6$TSS.enrichment > 2, 'High', 'Low')

nondoublets_18m_rep4$pct_reads_in_peaks <- nondoublets_18m_rep4$peak_region_fragments / nondoublets_18m_rep4$passed_filters * 100
nondoublets_18m_rep4$blacklist_ratio <- nondoublets_18m_rep4$blacklist_region_fragments / nondoublets_18m_rep4$peak_region_fragments
nondoublets_18m_rep4$high.tss <- ifelse(nondoublets_18m_rep4$TSS.enrichment > 2, 'High', 'Low')

nondoublets_18m_rep5$pct_reads_in_peaks <- nondoublets_18m_rep5$peak_region_fragments / nondoublets_18m_rep5$passed_filters * 100
nondoublets_18m_rep5$blacklist_ratio <- nondoublets_18m_rep5$blacklist_region_fragments / nondoublets_18m_rep5$peak_region_fragments
nondoublets_18m_rep5$high.tss <- ifelse(nondoublets_18m_rep5$TSS.enrichment > 2, 'High', 'Low')

nondoublets_18m_rep6$pct_reads_in_peaks <- nondoublets_18m_rep6$peak_region_fragments / nondoublets_18m_rep6$passed_filters * 100
nondoublets_18m_rep6$blacklist_ratio <- nondoublets_18m_rep6$blacklist_region_fragments / nondoublets_18m_rep6$peak_region_fragments
nondoublets_18m_rep6$high.tss <- ifelse(nondoublets_18m_rep6$TSS.enrichment > 2, 'High', 'Low')

dir.create("./Plots")
pdf("./Plots/TSSPlot_3m_rep4.pdf")
TSSPlot(nondoublets_3m_rep4, group.by = 'high.tss') + NoLegend()
dev.off()
pdf("./Plots/TSSPlot_3m_rep5.pdf")
TSSPlot(nondoublets_3m_rep5, group.by = 'high.tss') + NoLegend()
dev.off()
pdf("./Plots/TSSPlot_3m_rep6.pdf")
TSSPlot(nondoublets_3m_rep6, group.by = 'high.tss') + NoLegend()
dev.off()
pdf("./Plots/TSSPlot_18m_rep4.pdf")
TSSPlot(nondoublets_18m_rep4, group.by = 'high.tss') + NoLegend()
dev.off()
pdf("./Plots/TSSPlot_18m_rep5.pdf")
TSSPlot(nondoublets_18m_rep5, group.by = 'high.tss') + NoLegend()
dev.off()
pdf("./Plots/TSSPlot_18m_rep6.pdf")
TSSPlot(nondoublets_18m_rep6, group.by = 'high.tss') + NoLegend()
dev.off()

nondoublets_3m_rep4$nucleosome_group <- ifelse(nondoublets_3m_rep4$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
nondoublets_3m_rep5$nucleosome_group <- ifelse(nondoublets_3m_rep5$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
nondoublets_3m_rep6$nucleosome_group <- ifelse(nondoublets_3m_rep6$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
nondoublets_18m_rep4$nucleosome_group <- ifelse(nondoublets_18m_rep4$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
nondoublets_18m_rep5$nucleosome_group <- ifelse(nondoublets_18m_rep5$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
nondoublets_18m_rep6$nucleosome_group <- ifelse(nondoublets_18m_rep6$nucleosome_signal > 10, 'NS > 10', 'NS < 10')

pdf("./Plots/Fragment_hist_3m_rep4.pdf")
FragmentHistogram(object = nondoublets_3m_rep4, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
dev.off()

pdf("./Plots/Fragment_hist_3m_rep5.pdf")
FragmentHistogram(object = nondoublets_3m_rep5, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
dev.off()

pdf("./Plots/Fragment_hist_3m_rep6.pdf")
FragmentHistogram(object = nondoublets_3m_rep6, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
dev.off()

pdf("./Plots/Fragment_hist_18m_rep4.pdf")
FragmentHistogram(object = nondoublets_18m_rep4, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
dev.off()

pdf("./Plots/Fragment_hist_18m_rep5.pdf")
FragmentHistogram(object = nondoublets_18m_rep5, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
dev.off()

pdf("./Plots/Fragment_hist_18m_rep6.pdf")
FragmentHistogram(object = nondoublets_18m_rep6, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
dev.off()

pdf("./Plots/ViolinPlot_mm10_3M_rep4.pdf", width=20)
VlnPlot(
  object = nondoublets_3m_rep4,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("./Plots/ViolinPlot_mm10_3M_rep5.pdf", width=20)
VlnPlot(
  object = nondoublets_3m_rep5,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("./Plots/ViolinPlot_mm10_3M_rep6.pdf", width=20)
VlnPlot(
  object = nondoublets_3m_rep6,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("./Plots/ViolinPlot_mm10_18M_rep4.pdf", width=20)
VlnPlot(
  object = nondoublets_18m_rep4,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("./Plots/ViolinPlot_mm10_18M_rep5.pdf", width=20)
VlnPlot(
  object = nondoublets_18m_rep5,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf("./Plots/ViolinPlot_mm10_18M_rep6.pdf", width=20)
VlnPlot(
  object = nondoublets_18m_rep6,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

########################################################################

saveRDS(nondoublets_3m_rep4, "./RDS_files/mm10_3M_rep4_nondoublets.rds")
saveRDS(nondoublets_3m_rep5, "./RDS_files/mm10_3M_rep5_nondoublets.rds")
saveRDS(nondoublets_3m_rep6, "./RDS_files/mm10_3M_rep6_nondoublets.rds")
saveRDS(nondoublets_18m_rep4, "./RDS_files/mm10_18M_rep4_nondoublets.rds")
saveRDS(nondoublets_18m_rep5, "./RDS_files/mm10_18M_rep5_nondoublets.rds")
saveRDS(nondoublets_18m_rep6, "./RDS_files/mm10_18M_rep6_nondoublets.rds")

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = nondoublets_3m_rep4,
  y = list(nondoublets_3m_rep5, nondoublets_3m_rep6, nondoublets_18m_rep4, nondoublets_18m_rep5, nondoublets_18m_rep6),
  add.cell.ids = c("mm_3m_rep4", "mm_3m_rep5", "mm_3m_rep6", "mm_18m_rep4", "mm_18m_rep5", "mm_18m_rep6")
)
combined[["ATAC"]]
saveRDS(combined, "./RDS_files/combined_all_samples_merged_nondoublets.rds")

## ChromatinAssay data with 89951 features for 21688 cells
#........................................................................
# define a convenient function to load all the data and create a Seurat object
combined <- readRDS("./RDS_files/combined_all_samples_merged_nondoublets.rds")
combined_tmp <- combined
combined_tmp <- RunTFIDF(combined_tmp)
combined_tmp <- FindTopFeatures(combined_tmp)
combined_tmp <- RunSVD(combined_tmp)
combined_tmp_expt <- combined_tmp
combined_tmp_expt <- RunUMAP(combined_tmp_expt, dims = 2:10, reduction = 'lsi')
pdf("./Plots/DimPlot_scATAC_dataset_nondoublets_PC10.pdf")
DimPlot(combined_tmp_expt, group.by = 'dataset', pt.size = 0.1)
dev.off()

pdf("./Plots/DimPlot_scATAC_age_nondoublets_PC10.pdf")
DimPlot(combined_tmp_expt, group.by = 'age', pt.size = 0.1)
dev.off()

#pdf("../Plots/DimPlot_scATAC_batch_nondoublets_PC10.pdf")
#DimPlot(combined_tmp_expt, group.by = 'batch', pt.size = 0.1)
#dev.off()

#pdf("CoveragePlot_test_top5.pdf")
#CoveragePlot(
#  object = combined_tmp_expt,
#  group.by = 'dataset',
#  region = "chr14-99700000-99760000"
#)
#dev.off()

combined_tmp_expt <- FindNeighbors(
  object = combined_tmp_expt,
  reduction = 'lsi',
  dims = 2:10
)
combined_tmp_expt <- FindClusters(
  object = combined_tmp_expt,
  algorithm = 3,
  verbose = FALSE, future.seed=TRUE
)

png("Dimplot_clusters_2_10_nondoublets_PC10.png")
DimPlot(object = combined_tmp_expt, label = TRUE) + NoLegend()
dev.off()

gene.activities <- GeneActivity(combined_tmp_expt)
combined_tmp_expt[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined_tmp_expt <- NormalizeData(
  object = combined_tmp_expt,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined_tmp_expt$nCount_RNA)
)

DefaultAssay(combined_tmp_expt) <- 'RNA'

#pdf("FeaturePlot_top5.pdf")
#FeaturePlot(
#  object = combined_tmp_expt,
#  features = c('Sst','Pvalb',"Gad2","Neurod6","Rorb","Syt6"),
#  pt.size = 0.1,
#  max.cutoff = 'q95',
#  ncol = 3
#)
#dev.off()

#pdf("CoveragePlot_luminal.pdf")
#CoveragePlot(
#  object = combined_tmp_expt,
#  region = c("Krt18", "Elf5"),
#  extend.upstream = 1000,
#  extend.downstream = 1000,
#  ncol = 1
#)
#dev.off()

saveRDS(combined_tmp_expt, "scRNA_ATAC_combined_nondoublets_before_transfer.rds")

#########################################################

combined_tmp_expt <- readRDS("scRNA_ATAC_combined_nondoublets_before_transfer.rds")

mouse_rna <- readRDS("seuratObject_rename_ident_May26_2021.rds")

transfer.anchors <- FindTransferAnchors(
  reference = mouse_rna,
  query = combined_tmp_expt,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = mouse_rna@active.ident,
  weight.reduction = combined_tmp_expt[['lsi']],
  dims = 2:10
)

combined_tmp_expt <- AddMetaData(object = combined_tmp_expt, metadata = predicted.labels)

saveRDS(combined_tmp_expt, "scRNA_ATAC_combined_nondoublets_June22_2021.rds")

plot1 <- DimPlot(
  #object = mouse_rna, cols = c("#FFCC33", "#746CB1", "#1BBDC1","#09713A", "#3DB54A", "#BD77B2", "#A0D082", "#8AB6E1", "#262262", "#F3766E"),
  object = mouse_rna, cols = c("#FFCC33", "#FF9933","#746CB1", "#1BBDC1", "#09713A", "#3DB54A", "#c00000", "#A0D082", "#8AB6E1", "#ff99cc", "#262262", "#BD77B2"),
  group.by = 'ident',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  #object = combined_tmp_expt, cols = c("#746CB1", "#BD77B2", "#1BBDC1", "#3DB54A", "#A0D082", "#09713A", "#262262", "#FFCC33", "#8AB6E1"),
  object = combined_tmp_expt, cols = c("#746CB1", "#c00000", "#1BBDC1", "#3DB54A", "#A0D082", "#09713A", "#262262", "#FF9933", "#FFCC33", "#8AB6E1"),
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

pdf("DimPlot_RNA_ATAC_nondoublets_PC10_with_doublets_June22_2021.pdf", width=15)
plot1 + plot2
dev.off()

###########--------------------------------------###################

combined_tmp_expt <- subset(combined_tmp_expt, subset = predicted.id == "Doublet", invert = TRUE)
mouse_rna <- subset(mouse_rna, ident = "Doublet", invert = TRUE)
#levels(Idents(mouse_rna)) <- c("Bcells", "Dendritic/Macrophages", "Doublet", "Fibroblasts", "Luminal-AV", "Luminal-HS", "Myoepithelial", "Pericytes", "T/NKcells", "Vascular")  
plot1 <- DimPlot(
  object = mouse_rna, cols = c("#FFCC33", "#FF9933","#746CB1", "#1BBDC1", "#09713A", "#3DB54A", "#c00000", "#A0D082", "#8AB6E1", "#262262", "#BD77B2"),
  group.by = 'ident',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = combined_tmp_expt, cols = c("#746CB1", "#c00000", "#1BBDC1", "#3DB54A", "#A0D082", "#09713A", "#262262", "#FF9933", "#FFCC33", "#8AB6E1"),
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

pdf("DimPlot_RNA_ATAC_nondoublets_PC10_doublets_rmvd_June22_2021.pdf", width=15)
plot1 + plot2
dev.off()

pdf("DimPlot_ATAC_nondoublets_PC10_doublets_rmvd_June22_2021.pdf")
plot2
dev.off()

#Idents(scRNA_seq_expt_new) <- scRNA_seq_expt@meta.data$orig.ident

Idents(combined_tmp_expt) <- combined_tmp_expt$age

combined_tmp_expt$age <- factor(x = combined_tmp_expt$age, levels = c("3M", "18M"))
p9 <- DimPlot(combined_tmp_expt, cols = c("blue", "blue"), reduction = "umap", split.by = "age", label = TRUE)

pdf("Dimplot_split_age_ATAC.pdf", height=10, width=20)
plot_grid(p9)
dev.off()

combined_tmp_expt <- RenameIdents(
  object = combined_tmp_expt,
  '0' = 'Bcells',
  '1' = 'Myoepithelial',
  '2' = 'Bcells',
  '3' = 'Bcells',
  '4' = 'Bcells',
  '5' = 'Myoepithelial',
  '6' = 'T/NKcells',
  '7' = 'T/NKcells',
  '8' = 'Luminal-AV',
  '9' = 'T/NKcells',
  '10' = 'T/NKcells',
  '11' = 'Fibroblasts',
  '12' = 'Luminal-AV',
  '13' = 'T/NKcells',
  '14' = 'Myoepithelial',
  '15' = 'Luminal-HS',
  '16' = 'T/NKcells',
  '17' = 'Luminal-HS',
  '18' = 'Luminal-AV',
  '19' = 'Dendritic/Macrophages',
  '20' = 'Dendritic/Macrophages',
  '21' = 'T/NKcells',
  '22' = 'Vascular',
  '23' = 'Fibroblasts',
  '24' = 'Pericytes',
  '25' = 'Dendritic/Macrophages',
  '26' = 'Luminal-AV'
)

#Find differentially accessible peaks for clusters.
# change back to working with peaks instead of gene activities
DefaultAssay(combined_tmp_expt) <- 'ATAC'
saveRDS(combined_tmp_expt, "combined_after_renaming.rds")

combined_tmp_expt <- readRDS("./RDS_files/combined_after_renaming.rds")
combined_tmp_expt$age <- factor(x = combined_tmp_expt$age, levels = c("3M", "18M"))
p9 <- DimPlot(combined_tmp_expt, reduction = "umap", split.by = "age", label = TRUE)
p10 <- DimPlot(combined_tmp_expt, cols = c("blue", "blue"), reduction = "umap", split.by = "age", label = TRUE, pt.size=2)

pdf("Dimplot_split_age_and_orig_cluster_after_harmony_new2.pdf", height=10, width=20)
plot_grid(p9)
dev.off()

pdf("Dimplot_split_age_and_orig_cluster_after_harmony_ptsize2.pdf", height=10, width=20)
plot_grid(p10)
dev.off()

#Find differentially accessible peaks for clusters.
# change back to working with peaks instead of gene activities

combined_tmp_expt$celltype.stim <- paste(Idents(combined_tmp_expt), combined_tmp_expt$age, sep = "_")
combined_tmp_expt$celltype <- Idents(combined_tmp_expt)
Idents(combined_tmp_expt) <- "celltype.stim"
saveRDS(combined_tmp_expt, "seuratObject_rename_ident_DE_ATAC.rds")

combined <- readRDS("seuratObject_rename_ident_DE_ATAC.rds")

da_peaks_T_or_NKcells <- FindMarkers(
  object = combined,
  ident.1 = "T/NKcells_18M",
  ident.2 = "T/NKcells_3M",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks_T_or_NKcells)
write.table(da_peaks_T_or_NKcells, "DA_peaks_T_or_NKcells.txt", sep="\t", quote=F)

da_peaks_T_or_NKcells <- da_peaks_T_or_NKcells[da_peaks_T_or_NKcells$p_val_adj < 0.05,]

open_l23_da_peaks_T_or_NKcells_FC_0.5 <- rownames(da_peaks_T_or_NKcells[da_peaks_T_or_NKcells$avg_log2FC > 0.5, ])
close_l456_da_peaks_T_or_NKcells_FC_0.5 <- rownames(da_peaks_T_or_NKcells[da_peaks_T_or_NKcells$avg_log2FC < -0.5, ])
closest_genes_open_l23_da_peaks_T_or_NKcells_FC0.5 <- ClosestFeature(combined, regions = open_l23_da_peaks_T_or_NKcells_FC_0.5)
closest_genes_close_l456_da_peaks_T_or_NKcells_FC0.5 <- ClosestFeature(combined, regions = close_l456_da_peaks_T_or_NKcells_FC_0.5)

open_l23_da_peaks_T_or_NKcells_pos <- rownames(da_peaks_T_or_NKcells[da_peaks_T_or_NKcells$avg_log2FC > 0, ])
close_l456_da_peaks_T_or_NKcells_neg <- rownames(da_peaks_T_or_NKcells[da_peaks_T_or_NKcells$avg_log2FC < 0, ])
closest_genes_open_l23_da_peaks_T_or_NKcells_pos <- ClosestFeature(combined, regions = open_l23_da_peaks_T_or_NKcells_pos)
closest_genes_close_l456_da_peaks_T_or_NKcells_neg <- ClosestFeature(combined, regions = close_l456_da_peaks_T_or_NKcells_neg)

write.table(closest_genes_open_l23_da_peaks_T_or_NKcells_FC0.5, "open_l23_T_or_NKcells_log2FC0.5.txt", sep="\t", quote=F)
write.table(closest_genes_close_l456_da_peaks_T_or_NKcells_FC0.5, "close_l456_T_or_NKcells_log2FC0.5.txt", sep="\t", quote=F)

write.table(closest_genes_open_l23_da_peaks_T_or_NKcells_pos, "open_l23_T_or_NKcells_pos.txt", sep="\t", quote=F)
write.table(closest_genes_close_l456_da_peaks_T_or_NKcells_neg, "close_l456_T_or_NKcells_close.txt", sep="\t", quote=F)

#########################################################

