#Script to generate cell proportions plots based on celltypes for scRNA-seq.

#################-------------------------############################
library(Seurat)
library(ggplot2)
library(stringr)
library(reshape2)

gen_barplots <- function(scRNA_seq_expt){
	cells_per_ident_per_cluster = table(Idents(scRNA_seq_expt), scRNA_seq_expt@meta.data$orig.ident)
	write.table(cells_per_ident_per_cluster, "./Barplots_cell_composition/Table_cells_per_ident_per_cluster.txt", sep="\t", quote=FALSE)
	
	######################################################################################
	cells_per_ident_per_cluster_reps =table(Idents(scRNA_seq_expt), scRNA_seq_expt@meta.data$replicate)
	write.table(cells_per_ident_per_cluster_reps, "./Barplots_cell_composition/Table_cells_per_replicate_per_cluster.txt", sep="\t", quote=FALSE)

	#Convert to long format.
	cells_per_ident_per_cluster = read.table("./Barplots_cell_composition/Table_cells_per_ident_per_cluster.txt", sep="\t", header=TRUE)
	df_cells_per_cluster = as.data.frame.matrix(cells_per_ident_per_cluster)

	cells_per_ident_per_cluster_reps = read.table("./Barplots_cell_composition/Table_cells_per_replicate_per_cluster.txt", sep="\t", header=TRUE)
	df_cells_per_cluster_reps = as.data.frame.matrix(cells_per_ident_per_cluster_reps)

	Cell_type = rownames(df_cells_per_cluster)
	df_cells_per_cluster <- cbind(Cell_type, df_cells_per_cluster)

	Cell_type = rownames(df_cells_per_cluster_reps)
	df_cells_per_cluster_reps <- cbind(Cell_type, df_cells_per_cluster_reps)

	#Convert to long format.
	df_cells_long_format_reps_all = melt(df_cells_per_cluster_reps, id.vars=c("Cell_type"))
	df_cells_long_format_reps_all$Time_point <- factor(df_cells_long_format_reps_all$variable, levels = c("M3_rep4", "M3_rep5", "M3_rep6", "M18_rep4", "M18_rep5", "M18_rep6", "M3_rep1", "M3_rep2", "M3_rep3", "M18_rep1", "M18_rep2", "M18_rep3"))

	df_cells_long_format_reps_all$Cell_type <- factor(df_cells_per_cluster$Cell_type, levels = c("Tcells_naive","Tcells_mem", "Bcells", "Plasma", "Dendritic/Macrophages", "Luminal-AV", "Luminal-HS", "Myoepithelial", "Fibroblasts","Vascular","Pericytes", "Doublet"))

	names(df_cells_long_format_reps_all) = c("Cell_type", "variable", "prop", "Time_point")
	df_cells_long_format_reps_all_new <- df_cells_long_format_reps_all[df_cells_long_format_reps_all$Cell_type != "Doublet", ]

	pdf("./Barplots_cell_composition/Barplot_cell_cluster_per_ident_reps_all.pdf", width=12)
	print(ggplot(df_cells_long_format_reps_all_new, aes(fill=Cell_type, y=prop, x=Time_point)) + geom_bar(position="fill", stat="identity", color="black", size=0.01) + scale_fill_manual(values=c("#FF9933", "#FFCC33", "#746CB1", "#BD77B2", "#c00000", "#3DB54A", "#A0D082", "#09713A", "#1BBDC1", "#262262", "#8AB6E1")) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)))
	dev.off()

	##########################################################
	cells_per_ident_per_cluster_reps = read.table("./Barplots_cell_composition/Table_cells_per_replicate_per_cluster.txt", sep="\t", header=TRUE)
	df_cells_per_cluster_reps = as.data.frame.matrix(cells_per_ident_per_cluster_reps)

	Cell_type = rownames(df_cells_per_cluster)
	df_cells_per_cluster <- cbind(Cell_type, df_cells_per_cluster)

	Cell_type = rownames(df_cells_per_cluster_reps)
	df_cells_per_cluster_reps <- cbind(Cell_type, df_cells_per_cluster_reps)

	#Convert to long format.
	df_cells_long_format_reps_all = melt(df_cells_per_cluster_reps, id.vars=c("Cell_type"))
	df_cells_long_format_reps_all$Time_point <- factor(df_cells_long_format_reps_all$variable, levels = c("M3_rep4", "M3_rep5", "M3_rep6", "M18_rep4", "M18_rep5", "M18_rep6", "M3_rep1", "M3_rep2", "M3_rep3", "M18_rep1", "M18_rep2", "M18_rep3"))

	df_cells_long_format_reps_all$Cell_type <- factor(df_cells_long_format_reps_all$Cell_type, levels = c("Fibroblasts", "Pericytes", "Vascular", "Dendritic/Macrophages", "T/NKcells","Bcells", "Luminal-AV", "Luminal-HS","Myoepithelial"))
	names(df_cells_long_format_reps_all) = c("Cell_type", "variable", "prop", "Time_point")
	df_cells_long_format_reps_all_new <- df_cells_long_format_reps_all[df_cells_long_format_reps_all$Cell_type != "Doublet", ]
	df_cells_long_format_reps_all_new <- df_cells_long_format_reps_all_new[complete.cases(df_cells_long_format_reps_all_new), ]

	pdf("./Barplots_cell_composition/Barplot_cell_cluster_per_ident_reps_all_new.pdf", width=12)
	print(ggplot(df_cells_long_format_reps_all_new, aes(fill=Cell_type, y=prop, x=Time_point)) + geom_bar(position="fill", stat="identity", color="black", size=0.01) + scale_fill_manual(values=c("#1BBDC1", "#262262", "#8AB6E1", "#BD77B2", "#FFCC33", "#746CB1", "#3DB54A", "#A0D082", "#09713A")) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)))
	dev.off()
	
	write.table(df_cells_long_format_reps_all_new, "./Barplots_cell_composition/Table_cells_per_replicate_per_cluster_long.txt", sep="\t", quote=F)

	###################################################
	#Cell composition barplot for all immune cells.
	df_cells_per_cluster_immune_reps <- cells_per_ident_per_cluster_reps[c("Tcells_naive", "Tcells_mem", "Bcells", "Plasma", "Dendritic/Macrophages"),]
	Cell_type = rownames(df_cells_per_cluster_immune_reps)
	df_cells_per_cluster_immune_reps <- cbind(Cell_type, df_cells_per_cluster_immune_reps)

	#Convert to long format.
	library(reshape2)
	df_cells_long_format_reps_all = melt(df_cells_per_cluster_immune_reps, id.vars=c("Cell_type"))
	df_cells_long_format_reps_all$Time_point <- factor(df_cells_long_format_reps_all$variable, levels = c("M3_rep4", "M3_rep5", "M3_rep6", "M18_rep4", "M18_rep5", "M18_rep6","M3_rep1", "M3_rep2", "M3_rep3", "M18_rep1", "M18_rep2", "M18_rep3"))

	df_cells_long_format_reps_all$Cell_type <- factor(df_cells_per_cluster_immune_reps$Cell_type, levels = c("Dendritic/Macrophages", "Tcells_naive", "Tcells_mem", "Bcells", "Plasma"))
	names(df_cells_long_format_reps_all) = c("Cell_type", "variable", "prop", "Time_point")

	#Cell composition barplot for all immune cells (order by experiment).
	pdf("./Barplots_cell_composition/Barplot_cell_cluster_per_ident_rep_immune_order_expt.pdf", width=12)
	print(ggplot(df_cells_long_format_reps_all, aes(fill=Cell_type, y=prop, x=Time_point)) + geom_bar(position="fill", stat="identity", color="black", size=0.01) + scale_fill_manual(values=c("#c00000", "#FF9933", "#FFCC33", "#746CB1", "#BD77B2")) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)))
	dev.off()

	#############################################
	df_cells_long_format_reps_all$Time_point <- factor(df_cells_long_format_reps_all$variable, levels = c("M3_rep4", "M3_rep5", "M3_rep6", "M3_rep1", "M3_rep2", "M3_rep3", "M18_rep4", "M18_rep5", "M18_rep6", "M18_rep1", "M18_rep2", "M18_rep3"))
	df_cells_long_format_reps_all$Cell_type <- factor(df_cells_per_cluster_immune_reps$Cell_type, levels = c("Dendritic/Macrophages", "Tcells_naive", "Tcells_mem", "Bcells", "Plasma"))
	names(df_cells_long_format_reps_all) = c("Cell_type", "variable", "prop", "Time_point")

	#Barplots for immune cells with labels (order by age).
	pdf("./Barplots_cell_composition/Barplot_cell_cluster_per_ident_rep_immune_order_age.pdf", width=12)
	print(ggplot(df_cells_long_format_reps_all, aes(fill=Cell_type, y=prop, x=Time_point)) + geom_bar(position="fill", stat="identity", colour="black") + scale_fill_manual(values=c("#c00000", "#FF9933", "#FFCC33", "#746CB1", "#BD77B2")) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)))
	dev.off()

	#Barplots for immune cells with labels (order by age) with labels.
	pdf("./Barplots_cell_composition/Barplot_cell_cluster_per_ident_rep_immune_order_age_with_label.pdf", width=12)
	print(ggplot(df_cells_long_format_reps_all, aes(fill=Cell_type, y=prop, x=Time_point)) + geom_bar(position="fill", stat="identity", colour="black") + scale_fill_manual(values=c("#c00000", "#FF9933", "#FFCC33", "#746CB1", "#BD77B2")) + geom_text(aes(label=prop), position=position_fill(vjust=0.5)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)))
	dev.off()

	#########################################
	#Cell composition barplot for epithelial cells.
	df_cells_per_cluster_epithelial_reps <- cells_per_ident_per_cluster_reps[c("Luminal-AV", "Luminal-HS", "Myoepithelial"),]
	Cell_type = rownames(df_cells_per_cluster_epithelial_reps)
	df_cells_per_cluster_epithelial_reps <- cbind(Cell_type, df_cells_per_cluster_epithelial_reps)

	#Convert to long format.
	library(reshape2)
	df_cells_long_format_reps_all = melt(df_cells_per_cluster_epithelial_reps, id.vars=c("Cell_type"))
	df_cells_long_format_reps_all$Time_point <- factor(df_cells_long_format_reps_all$variable, levels = c("M3_rep4", "M3_rep5", "M3_rep6", "M18_rep4", "M18_rep5", "M18_rep6","M3_rep1", "M3_rep2", "M3_rep3", "M18_rep1", "M18_rep2", "M18_rep3"))

	df_cells_long_format_reps_all$Cell_type <- factor(df_cells_per_cluster_epithelial_reps$Cell_type, levels = c("Luminal-AV", "Luminal-HS", "Myoepithelial"))

	#Cell composition barplot for epithelial cells ( by experiment).`
	names(df_cells_long_format_reps_all) = c("Cell_type", "variable", "prop", "Time_point")
	pdf("./Barplots_cell_composition/Barplot_cell_cluster_per_ident_rep_epithelial_order_expt.pdf", width=12)
	print(ggplot(df_cells_long_format_reps_all, aes(fill=Cell_type, y=prop, x=Time_point)) + geom_bar(position="fill", stat="identity", color="black", size=0.01) + scale_fill_manual(values=c("#3DB54A", "#A0D082", "#09713A")) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)))
	dev.off()

	#Cell composition barplot for epithelial cells ( by experiment) with labels.`
	pdf("./Barplots_cell_composition/Barplot_cell_cluster_per_ident_rep_epithelial_order_expt_with_label.pdf", width=12)
	print(ggplot(df_cells_long_format_reps_all, aes(fill=Cell_type, y=prop, x=Time_point)) + geom_bar(position="fill", stat="identity", color="black") + scale_fill_manual(values=c("#3DB54A", "#A0D082", "#09713A")) + geom_text(aes(label=prop), position=position_fill(vjust=0.5)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  + scale_y_continuous(expand = c(0,0))) 
	dev.off()

	############3---------------##################

	df_cells_long_format_reps_all$Time_point <- factor(df_cells_long_format_reps_all$variable, levels = c("M3_rep4", "M3_rep5", "M3_rep6", "M3_rep1", "M3_rep2", "M3_rep3", "M18_rep4", "M18_rep5", "M18_rep6", "M18_rep1", "M18_rep2", "M18_rep3"))
	df_cells_long_format_reps_all$Cell_type <- factor(df_cells_per_cluster_epithelial_reps$Cell_type, levels = c("Luminal-AV", "Luminal-HS", "Myoepithelial"))
	names(df_cells_long_format_reps_all) = c("Cell_type", "variable", "prop", "Time_point")

	#Cell composition barplot for all epithelial cells (order by age).
	pdf("./Barplots_cell_composition/Barplot_cell_cluster_per_ident_rep_epithelial_order_age.pdf", width=12)
	print(ggplot(df_cells_long_format_reps_all, aes(fill=Cell_type, y=prop, x=Time_point)) + geom_bar(position="fill", stat="identity", colour="black") + scale_fill_manual(values=c("#3DB54A", "#A0D082", "#09713A")) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)))
	dev.off()

	#Cell composition barplot for all epithelial cells (order by age) with number labels.
	pdf("./Barplots_cell_composition/Barplot_cell_cluster_per_ident_rep_epithelial_order_age_with_label.pdf", width=12)
	print(ggplot(df_cells_long_format_reps_all, aes(fill=Cell_type, y=prop, x=Time_point)) + geom_bar(position="fill", stat="identity", colour="black") + scale_fill_manual(values=c("#3DB54A", "#A0D082", "#09713A")) + geom_text(aes(label=prop), position=position_fill(vjust=0.5)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)))
	dev.off()

	################################################
	#Cell composition barplot for stromal cells.
	df_cells_per_cluster_stromal_reps <- cells_per_ident_per_cluster_reps[c("Fibroblasts", "Vascular", "Pericytes"),]
	Cell_type = rownames(df_cells_per_cluster_stromal_reps)
	df_cells_per_cluster_stromal_reps <- cbind(Cell_type, df_cells_per_cluster_stromal_reps)

	#Convert to long format.
	library(reshape2)
	df_cells_long_format_reps_all = melt(df_cells_per_cluster_stromal_reps, id.vars=c("Cell_type"))
	df_cells_long_format_reps_all$Time_point <- factor(df_cells_long_format_reps_all$variable, levels = c("M3_rep4", "M3_rep5", "M3_rep6", "M18_rep4", "M18_rep5", "M18_rep6","M3_rep1", "M3_rep2", "M3_rep3", "M18_rep1", "M18_rep2", "M18_rep3"))

	df_cells_long_format_reps_all$Cell_type <- factor(df_cells_per_cluster_stromal_reps$Cell_type, levels = c("Fibroblasts", "Pericytes", "Vascular"))
	names(df_cells_long_format_reps_all) = c("Cell_type", "variable", "prop", "Time_point")

	#Cell composition barplot for all stromal cells (order by experiment).
	pdf("./Barplots_cell_composition/Barplot_cell_cluster_per_ident_rep_stromal_order_expt.pdf", width=12)
	print(ggplot(df_cells_long_format_reps_all, aes(fill=Cell_type, y=prop, x=Time_point)) + geom_bar(position="fill", stat="identity", color="black", size=0.01) + scale_fill_manual(values=c("#1BBDC1", "#262262", "#8AB6E1")) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)))
	dev.off()

	#Cell composition barplot for all stromal cells with labels (order by expt).
	pdf("./Barplots_cell_composition/Barplot_cell_cluster_per_ident_rep_stromal_order_expt_with_label.pdf", width=12)
	print(ggplot(df_cells_long_format_reps_all, aes(fill=Cell_type, y=prop, x=Time_point)) + geom_bar(position="fill", stat="identity", color="black", size=0.01) + scale_fill_manual(values=c("#1BBDC1", "#262262", "#8AB6E1")) + geom_text(aes(label=prop), position=position_fill(vjust=0.5)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)))
	dev.off()

	df_cells_long_format_reps_all$Time_point <- factor(df_cells_long_format_reps_all$variable, levels = c("M3_rep4", "M3_rep5", "M3_rep6", "M3_rep1", "M3_rep2", "M3_rep3", "M18_rep4", "M18_rep5", "M18_rep6", "M18_rep1", "M18_rep2", "M18_rep3"))
	df_cells_long_format_reps_all$Cell_type <- factor(df_cells_per_cluster_stromal_reps$Cell_type, levels = c("Fibroblasts", "Pericytes", "Vascular"))
	names(df_cells_long_format_reps_all) = c("Cell_type", "variable", "prop", "Time_point")

	#Cell composition barplot for all other stromal cells (order by age).
	pdf("./Barplots_cell_composition/Barplot_cell_cluster_per_ident_rep_stromal_order_age.pdf", width=12)
	print(ggplot(df_cells_long_format_reps_all, aes(fill=Cell_type, y=prop, x=Time_point)) + geom_bar(position="fill", stat="identity", color="black") + scale_fill_manual(values=c("#1BBDC1", "#262262", "#8AB6E1")) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)))
	dev.off()

	#Cell composition barplot for all other stromal cells (order by age) with labels.
	pdf("./Barplots_cell_composition/Barplot_cell_cluster_per_ident_rep_stromal_order_age_with_label.pdf", width=12)
	print(ggplot(df_cells_long_format_reps_all, aes(fill=Cell_type, y=prop, x=Time_point)) + geom_bar(position="fill", stat="identity", color="black") + scale_fill_manual(values=c("#1BBDC1", "#262262", "#8AB6E1")) + geom_text(aes(label=prop), position=position_fill(vjust=0.5)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)))
	dev.off()

}

#Call function
dir.create("Barplots_cell_composition")
args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("Seurat object with cluster annotations must be supplied (input file).n", call.=FALSE)
} else
{
print("Generating cell composition barplots")
}

scRNA_seq_expt = readRDS(args[1]) #seuratObject_rename_ident.rds

#Call function to generate cell composition barplots.
gen_barplots(scRNA_seq_expt)

