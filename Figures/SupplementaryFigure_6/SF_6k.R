#Script to calculate score and generate Violin plots.

library(Seurat)
library(ggpubr)
library(rstatix)

scoring_plot <- function(mydata_new, module_type, module_filename){
	mydata_new$celltype.stim <- paste(Idents(mydata_new), mydata_new$orig.ident, sep = "_")
	mydata_new$celltype <- Idents(mydata_new)
	Idents(mydata_new) <- "celltype.stim"

	# read in inflammation genes
	inflamm <- read.delim(module_filename, header = T, check.names = F, stringsAsFactors = F)
	names(inflamm) <- c("Module.Name", "GeneName")
	# which inflammation genes are present in your seurat object
	length(which(inflamm$GeneName %in% rownames(mydata_new)))
	norm_data <- GetAssayData(mydata_new, slot = "data")
	
	id_inf <- which(inflamm$GeneName %in% rownames(mydata_new))
	# subset matrix by this data
	inflamm_sel <- inflamm[id_inf,]
	norm_sel <- norm_data[inflamm_sel$GeneName, ]
	
	# calculate sum of inflammation genes in each cell; divide by total number of transcripts
	norm_sel_df <- as.data.frame(norm_sel)
	norm_data_df <- as.data.frame(norm_data)
	cs_sel <- colSums(norm_sel_df)
	cs_all <- colSums(norm_data_df)
	cs_ratio <- cs_sel / cs_all
	# min-max scale between 0 and 1
	cs_min_max <- (cs_ratio-min(cs_ratio))/(max(cs_ratio)-min(cs_ratio))
	# add score to seurat object
	mydata_new <- AddMetaData(mydata_new, metadata = cs_min_max, col.name = module_type)

	mydata_new2 = mydata_new
	df <- as.data.frame(mydata_new2@meta.data)
	df_subset <- df[c("orig.ident", "celltype", module_type)]

	featureplot_filename <- paste0("./Scoring_plots/FeaturePlot_score_",module_type, ".pdf")
	pdf(featureplot_filename)
	print(FeaturePlot(mydata_new, features = module_type))
	dev.off()
	
	my_levels <- c("mm10_3mths", "mm10_18mths")
	mydata_new@meta.data$orig.ident <- factor(x = mydata_new@meta.data$orig.ident, levels = my_levels)
	
	my_levels_celltype <- c("Luminal-AV", "Luminal-HS", "Myoepithelial")
	
	mydata_new@meta.data$celltype <- factor(x = mydata_new@meta.data$celltype, levels = my_levels_celltype)

	#File with significance level.
	filename1 <- paste0("./Scoring_plots/VlnPlot_", module_type, ".pdf")
	pdf(filename1)
	print(VlnPlot(mydata_new, features = module_type, split.by = "orig.ident", group.by = "celltype", pt.size = 0, combine = TRUE) + stat_compare_means(method = "wilcox.test", label = "p.signif"))
	dev.off()
	
	#File with significance level.
	filename1_n <- paste0("./Scoring_plots/VlnPlot_", module_type, "_pformat.pdf")
	pdf(filename1_n)
	print(VlnPlot(mydata_new, features = module_type, split.by = "orig.ident", group.by = "celltype", pt.size = 0, combine = TRUE) + stat_compare_means(method = "wilcox.test", label = "p.format"))
	dev.off()
	
	#File with p-value	
	filename1_new <- paste0("./Scoring_plots/VlnPlot_", module_type, "_pval.pdf")
	pdf(filename1_new)
	print(VlnPlot(mydata_new, features = module_type, split.by = "orig.ident", group.by = "celltype", pt.size = 0, combine = TRUE) + stat_compare_means(method = "wilcox.test"))
	dev.off()
	
	#File split by age.
	filename2 <- paste0("./Scoring_plots/VlnPlot_", module_type, "_split", ".pdf")
	pdf(filename2)
	print(VlnPlot(mydata_new, features = module_type, split.by = "orig.ident", group.by = "celltype", pt.size = 0, combine = TRUE, split.plot=TRUE))
	dev.off()
	
	################################################################
	
	filename3 <- paste0("./Scoring_plots/VlnPlot_split_", module_type, "_log", ".pdf")
	pdf(filename3)
	print(VlnPlot(mydata_new, features = module_type, split.by = "orig.ident", group.by = "celltype", pt.size = 0, combine = FALSE, log=T, split.plot = TRUE))
	dev.off()
	
	filename4 <- paste0("./Scoring_plots/VlnPlot_", module_type, "log", ".pdf")
	pdf(filename4)
	print(VlnPlot(mydata_new, features = module_type, split.by = "orig.ident", group.by = "celltype", pt.size = 0, combine = FALSE, log=T))
	dev.off()
}

#Call function
dir.create("Scoring_plots")
args <- commandArgs(TRUE)
if (length(args) < 3) {
  stop("Seurat object with cluster annotations, module name and module filename must be supplied (input file).n", call.=FALSE)
} else
{
print("Generating Violin plots for module scoring...")
}

#mydata_new = readRDS("seuratObject_rename_ident_May26_2021.rds")
#module_type = "EMT_genelist"
#module_filename = "./Score_gene_plots/EMT_genelist_mouse.txt"

mydata = readRDS(args[1]) # load in your seurat dataset
module_type = args[2] # module name to be included in output filenames
module_filename = args[3] # path to module geneset file.

## Change Identities accordingly. For subsets, ensure right objects are loaded in and then identities are set to the right objects.
mydata_new <- subset(mydata, idents = c("Luminal-AV", "Luminal-HS", "Myoepithelial"))

#Call function
scoring_plot(mydata_new, module_type, module_filename) 


