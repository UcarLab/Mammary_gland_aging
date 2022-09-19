
#Script to perform differential expression analysis at single cell level.
########################################################

library(Seurat)
library(purrr)
library(dplyr)
dir.create("DE_analysis")

get_scDE_results <- function(num, A, B, x, y){
contrast = "18_vs_3"

celltype = paste0(mydata_all[num], "_", contrast)
celltype <- gsub("/", "_", celltype)

ident1 = paste0(mydata_all[num], "_", A)
ident2 = paste0(mydata_all[num], "_", B)

print(paste0("DE analysis for ", mydata_all[num], "...."))
var1 <- paste0(celltype, "_log1.markers")
mydata = y
var1.markers <- FindMarkers(mydata, ident.1 = ident1, ident.2 = ident2, test.use = 'LR', latent.vars = 'expt', verbose = FALSE, logfc.threshold = 0)

var1_file <- paste0("./DE_analysis/Diff_",var1, ".txt")
write.table(var1.markers, var1_file, sep="\t", quote=FALSE)

#celltype_log1.markers <- read.table("./DE_analysis/Diff_T_or_NKcells_18_vs_3mths_log1_LR.txt", sep="\t", header=T)
var1.markers$gene = rownames(var1.markers)
celltype_DE_log0_pos = var1.markers %>% filter(avg_log2FC > 0, p_val_adj <= 0.05)
celltype_DE_log0_neg = var1.markers %>% filter(avg_log2FC < 0, p_val_adj <= 0.05)
celltype_DE_log0_neg <- celltype_DE_log0_neg[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
celltype_DE_log0_pos <- celltype_DE_log0_pos[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]

file1_log0_neg <- paste0("./DE_analysis/", celltype, "_DE_neg_log0_LR.txt")
file1_log0_pos <- paste0("./DE_analysis/", celltype, "_DE_pos_log0_LR.txt")
write.table(celltype_DE_log0_neg, file1_log0_neg, sep="\t", quote=FALSE)
write.table(celltype_DE_log0_pos, file1_log0_pos, sep="\t", quote=FALSE)

DE_genes_heatmap <- c(celltype_DE_log0_neg$gene, celltype_DE_log0_pos$gene)

file_heatmap <- paste0("Heatmap_single_cell_", celltype, ".pdf")

pdf(file_heatmap)
DoHeatmap(mydata, features = DE_genes_heatmap) + NoLegend()
dev.off()

################---------------------#################################

celltype_log0.25.markers <- FindMarkers(mydata, ident.1 = ident1, ident.2 = ident2, test.use = 'LR', latent.vars = 'expt', verbose = FALSE)
var0.25 <- paste(celltype, "_log0.25.markers", sep="")
var0.25_file <- paste0("./DE_analysis/", "Diff_", celltype, "_log0.25_LR.txt")
write.table(celltype_log0.25.markers, var0.25_file, sep="\t", quote=FALSE)

celltype_log0.25.markers$gene = rownames(celltype_log0.25.markers)
celltype_DE_log0.25_pos = celltype_log0.25.markers %>% filter(avg_log2FC > 0, p_val_adj <= 0.05)
celltype_DE_log0.25_neg = celltype_log0.25.markers %>% filter(avg_log2FC < 0, p_val_adj <= 0.05)
celltype_DE_log0.25_neg <- celltype_DE_log0.25_neg[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
celltype_DE_log0.25_pos <- celltype_DE_log0.25_pos[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
file_DE_log0.25_neg <- paste0("./DE_analysis/", celltype, "_neg_log0.25_LR.txt")
file_DE_log0.25_pos <- paste0("./DE_analysis/", celltype, "_pos_log0.25_LR.txt")
write.table(celltype_DE_log0.25_neg, file_DE_log0.25_neg, sep="\t", quote=FALSE)
write.table(celltype_DE_log0.25_pos, file_DE_log0.25_pos, sep="\t", quote=FALSE)

celltype_DE_logFC_0.5_pos = celltype_log0.25.markers %>% filter(avg_log2FC > 0.5, p_val_adj <= 0.05)
celltype_DE_logFC_0.5_neg = celltype_log0.25.markers %>% filter(avg_log2FC < -0.5, p_val_adj <= 0.05)
file_DE_log0.5_neg <- paste0("./DE_analysis/", celltype, "_neg_log0.5_LR.txt")
file_DE_log0.5_pos <- paste0("./DE_analysis/", celltype, "_pos_log0.5_LR.txt")
write.table(celltype_DE_logFC_0.5_neg, file_DE_log0.5_neg, sep="\t", quote=FALSE)
write.table(celltype_DE_logFC_0.5_pos, file_DE_log0.5_pos, sep="\t", quote=FALSE)

celltype_DE_logFC_1.0_pos = celltype_log0.25.markers %>% filter(avg_log2FC > 1, p_val_adj <= 0.05)
celltype_DE_logFC_1.0_neg = celltype_log0.25.markers %>% filter(avg_log2FC < -1, p_val_adj <= 0.05)
file_DE_log1.0_neg <- paste0("./DE_analysis/", celltype, "_neg_log1.0_LR.txt")
file_DE_log1.0_pos <- paste0("./DE_analysis/", celltype, "_pos_log1.0_LR.txt")
write.table(celltype_DE_logFC_1.0_neg, file_DE_log1.0_neg, sep="\t", quote=FALSE)
write.table(celltype_DE_logFC_1.0_pos, file_DE_log1.0_pos, sep="\t", quote=FALSE)


}

#############--------------------------------########################
args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else
{
print("Running differential analysis...")
}

file_obj = readRDS(args[1])
mydata <- subset(file_obj, idents = c("Doublet"), invert = TRUE )
mydata$celltype.stim <- paste(Idents(mydata), mydata$orig.ident, sep = "_")
mydata$celltype <- Idents(mydata)
Idents(mydata) <- "celltype.stim"
saveRDS(mydata, "seuratObject_rename_ident_DE_tp.rds")

mydata_all <- names(table(mydata$celltype))

map(1:length(mydata_all), get_scDE_results, A= "mm10_18mths", B="mm10_3mths", x = mydata_all, y = mydata)


