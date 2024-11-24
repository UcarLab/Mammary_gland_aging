
#################-------------------------############################
library(Seurat)
library(ggplot2)
library(stringr)
library(reshape2)
library(dplyr)

###############################################################

scATAC_seq_expt <- readRDS("scRNA_ATAC_combined_nondoublets_June22_2021.rds")

	Idents(scATAC_seq_expt) <- scATAC_seq_expt$predicted.id
	scATAC_seq_expt_epithelial_only <- subset(scATAC_seq_expt, idents = c("Luminal-AV", "Luminal-HS", "Myoepithelial"))
	epithelial_cellfreq <- table(Idents(scATAC_seq_expt_epithelial_only), scATAC_seq_expt_epithelial_only$dataset)
	write.table(epithelial_cellfreq, "./Barplots_cell_composition/Epithelial_cellfreq_summary.txt", sep="\t", quote=F)

	test_epithelial_only <- prop.table(table(Idents(scATAC_seq_expt_epithelial_only), scATAC_seq_expt_epithelial_only$dataset), margin=2)
	write.table(test_epithelial_only, "./Barplots_cell_composition/Epithelial_prop_summary.txt", sep="\t", quote=F)

	test_epithelial_only_df <- as.data.frame(test_epithelial_only)
	names(test_epithelial_only_df) = c("Cell_type", "Sample", "prop")
	test_epithelial_only_df$Age <- gsub("_rep.*$", "",test_epithelial_only_df$Sample)

	LumHS <- filter(test_epithelial_only_df, Cell_type == "Luminal-HS")
	res_LumHS <- t.test(prop ~ Age, data = LumHS, paired=T, exact=FALSE)
        print(paste0("Statistical analysis for Luminal-HS"))
        print(res_LumHS)
	
	LumAV <- filter(test_epithelial_only_df, Cell_type == "Luminal-AV")
	res_LumAV <- t.test(prop ~ Age, data = LumAV, paired=T, exact=FALSE)
        print(paste0("Statistical analysis for Luminal-AV"))
        print(res_LumAV)

	Myoepithelial <- filter(test_epithelial_only_df, Cell_type == "Myoepithelial")
	res_Myoepithelial <- t.test(prop ~ Age, data = Myoepithelial, paired=T, exact=FALSE)
        print(paste0("Statistical analysis for Myoepithelial."))
        print(res_Myoepithelial)

	###############################################

	scATAC_seq_expt_immune_only <- subset(scATAC_seq_expt, idents = c("Dendritic/Macrophages", "Tcells_naive", "Tcells_mem", "Bcells"))
	immune_cellfreq <- table(scATAC_seq_expt_immune_only$predicted.id, scATAC_seq_expt_immune_only$dataset)
	write.table(immune_cellfreq, "./Barplots_cell_composition/Immune_cellfreq_summary.txt", sep="\t", quote=F)

	test_immune_only <- prop.table(table(Idents(scATAC_seq_expt_immune_only), scATAC_seq_expt_immune_only$dataset), margin=2)
	write.table(test_immune_only, "./Barplots_cell_composition/Immune_prop_summary.txt", sep="\t", quote=F)

	test_immune_only_df <- as.data.frame(test_immune_only)
	names(test_immune_only_df) = c("Cell_type", "Sample", "prop")
	test_immune_only_df$Age <- gsub("_rep.*$", "",test_immune_only_df$Sample)

	DC_Macro <- filter(test_immune_only_df, Cell_type == "Dendritic/Macrophages")
	res_DC_Macro <- t.test(prop ~ Age, data = DC_Macro, paired=T, exact=FALSE)
        print(paste0("Statistical analysis for Dendritic/Macro."))
        print(res_DC_Macro)

	Tcells_naive <- filter(test_immune_only_df, Cell_type == "Tcells_naive")
	res_Tcells_naive <- t.test(prop ~ Age, data = Tcells_naive, paired=T, exact=FALSE)
        print(paste0("Statistical analysis for Tcells_naive."))
        print(res_Tcells_naive)
	
	Tcells_mem <- filter(test_immune_only_df, Cell_type == "Tcells_mem")
	res_Tcells_mem <- t.test(prop ~ Age, data = Tcells_mem, paired=T, exact=FALSE)
        print(paste0("Statistical analysis for Tcells_mem."))
        print(res_Tcells_mem)

	Bcells <- filter(test_immune_only_df, Cell_type == "Bcells")
	res_Bcells <- t.test(prop ~ Age, data = Bcells, paired=T, exact=FALSE)
        print(paste0("Statistical analysis for Bcells."))
        print(res_Bcells)
	
	#############################################

	scATAC_seq_expt_stromal_only <- subset(scATAC_seq_expt, idents = c("Fibroblasts", "Vascular", "Pericytes"))
	stromal_cellfreq <- table(scATAC_seq_expt_stromal_only$predicted.id, scATAC_seq_expt_stromal_only$dataset)
	write.table(stromal_cellfreq, "./Barplots_cell_composition/Stromal_cellfreq_summary.txt", sep="\t", quote=F)
	
	test_stromal_only <- prop.table(table(Idents(scATAC_seq_expt_stromal_only), scATAC_seq_expt_stromal_only$dataset), margin=2)
	write.table(test_stromal_only, "./Barplots_cell_composition/Stromal_prop_summary.txt", sep="\t", quote=F)

	test_stromal_only_df <- as.data.frame(test_stromal_only)
	names(test_stromal_only_df) = c("Cell_type", "Sample", "prop")
	test_stromal_only_df$Age <- gsub("_rep.*$", "",test_stromal_only_df$Sample)
	
	Fibroblasts <- filter(test_stromal_only_df, Cell_type == "Fibroblasts")
	res_Fibro <- t.test(prop ~ Age, data = Fibroblasts, paired=T, exact=FALSE)
        print(paste0("Statistical analysis for Fibroblasts."))
        print(res_Fibro)
	
	Vascular <- filter(test_stromal_only_df, Cell_type == "Vascular")
	res_Vascular <- t.test(prop ~ Age, data = Vascular, paired=T, exact=FALSE)
        print(paste0("Statistical analysis for Vascular."))
        print(res_Vascular)
	
	Pericytes <- filter(test_stromal_only_df, Cell_type == "Pericytes")
	res_Pericytes <- t.test(prop ~ Age, data = Pericytes, paired=T, exact=FALSE)
        print(paste0("Statistical analysis for Pericytes."))
        print(res_Pericytes)


