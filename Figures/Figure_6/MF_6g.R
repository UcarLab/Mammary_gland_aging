library(tidyverse)
library(JASPAR2020)
library(TFBSTools)
library(universalmotif)
library(dendextend)
library(grid)
library(gridExtra)
library(ggside)
library(patchwork)

# read in the differential motifs dataframe from Signac Vignette
chromVAR.signature.motifs <- read.csv("~/path/to/chromvar/files", sep = '\t', header = T)
colnames(chromVAR.signature.motifs) <- c("ID", "p_val", "log_FC_Change", "pct.1", "pct.2", "adj_p_val", "CellType")
chromVAR.signature.motifs$CellType <- gsub(".txt","",chromVAR.signature.motifs$CellType)
chromVAR.signature.motifs$log_FC_Change <- as.numeric(chromVAR.signature.motifs$log_FC_Change)
chromVAR.signature.motifs <- chromVAR.signature.motifs[(chromVAR.signature.motifs$log_FC_Change > 0),]

# create a motif database in order to cluster motifs
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# convert motif ids to motif names using database
chromVAR.signature.motifs <- chromVAR.signature.motifs %>%
  rowwise() %>%
  mutate(motif = pfm[[ID]]@name) %>%
  as.data.frame()

#Choose list of top motifs to plot.
top_motifs <- chromVAR.signature.motifs %>%
  dplyr::group_by(CellType) %>%
  slice_min(order_by = adj_p_val, n = 20, with_ties = F) %>%
  as.data.frame() %>%
  dplyr::select(motif,ID)

top_motifs <- top_motifs %>%
  distinct()

# get the top motifs as Universal Motif objects so that we can filter and cluster them
motifs <- universalmotif::convert_motifs(pfm[top_motifs$ID])

# find euclidian distance of motifs
comp.matrix <- universalmotif::compare_motifs(motifs, method = "EUCL")

# convert to suitable format for dhc
comp.dist <- as.dist(comp.matrix)

# hierarchical clustering of the motifs for pheatmap#
dhc <- hclust(comp.dist)

# just get the top motifs from results and also convert adjusted p value to -log for plotting
# Adjust Celltypes for ordering

Celltypes = c('Mye-C1','Mye-C2','Mye-C3','Mye-C4','Mye-C5','Mye-C6','Mye-C7')
df <- chromVAR.signature.motifs %>%
  mutate(CellType = factor(CellType, levels = Celltypes)) %>%
  mutate(Log_Adj_P = -log10(adj_p_val)) %>%
  mutate(Log_Adj_P = ifelse(Log_Adj_P > 100, 100, Log_Adj_P)) %>%
  filter(motif %in% dhc$labels) %>%
  as.data.frame()

# convert data frame to matrix for pheatmap
a <- df[,c("CellType","motif","Log_Adj_P")] %>%
  distinct(CellType, motif, .keep_all = T) %>%
  spread(key = "CellType", value = "Log_Adj_P")

rownames(a) <- a$motif
b <- a[,-1]
c <- as.data.frame(t(b))

pdf("~/path/to/pdf", height = 15, width = 5)
pheatmap::pheatmap(b[dhc$labels,], 
                   cluster_rows = F,
                   cluster_cols = F,
                   scale = "none", 
                   na_col = "grey",
                   color = colorRampPalette(c("lightblue1","#56B1F7", "#132B43"))(50),
                   fontsize_number = 12,
                   show_colnames = T,
                   legend = T)
dev.off()


