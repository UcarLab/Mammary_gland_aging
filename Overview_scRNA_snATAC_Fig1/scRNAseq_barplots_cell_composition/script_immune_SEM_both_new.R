library(ggplot2)
library(plyr)
library(reshape2)
library(tidyr)
library(reshape2)
library(stringr)
library(Seurat)

#Call function
dir.create("Barplots_cell_composition")
args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("Seurat object with cluster annotations must be supplied (input file).n", call.=FALSE)
} else
{
print("Generating barplot with SEM for immune cells...")
}

scRNA_seq_expt = readRDS(args[1])
#scRNA_seq_expt <- readRDS("/projects/anczukow-lab/neerja_projects/Integrate_expt1_expt2_only/Analysis_PC_10_mt10/RDS_files/seuratObject_rename_ident.rds")
cells_per_ident_per_cluster_reps =table(Idents(scRNA_seq_expt), scRNA_seq_expt@meta.data$replicate)
write.table(cells_per_ident_per_cluster_reps, "./Barplots_cell_composition/Table_cells_per_replicate_per_cluster.txt", sep="\t", quote=FALSE)

#########################################
#Cell composition barplot for immune cells.
df_cells_per_cluster_immune_reps <- cells_per_ident_per_cluster_reps[c("Dendritic/Macrophages", "Tcells_naive", "Tcells_mem", "Bcells", "Plasma"),]
Cell_type = rownames(df_cells_per_cluster_immune_reps)
df_cells_per_cluster_immune_reps <- cbind(Cell_type, df_cells_per_cluster_immune_reps)
df_cells_per_cluster_immune_reps = as.data.frame.matrix(df_cells_per_cluster_immune_reps)

Cell_type = rownames(df_cells_per_cluster_immune_reps)
df_cells_per_cluster_immune_reps <- cbind(Cell_type, df_cells_per_cluster_immune_reps)
df_cells_per_cluster_immune_reps = as.data.frame.matrix(df_cells_per_cluster_immune_reps)

Cell_type = rownames(df_cells_per_cluster_immune_reps)
df_cells_per_cluster_immune_reps <- cbind(Cell_type, df_cells_per_cluster_immune_reps)

#Convert to long format.
df_cells_long_format_reps_all = melt(df_cells_per_cluster_immune_reps, id.vars=c("Cell_type"))
replacements <- str_replace_all(df_cells_long_format_reps_all[,2], "_rep.+","")
df_cells_long_format_reps_all$Time <- replacements

names(df_cells_long_format_reps_all) <- c("Cell_type", "Rep_age", "value", "Time")

data_wide <- spread(df_cells_long_format_reps_all, Cell_type, value)
names(data_wide) <- c("Rep_age", "Time", "Bcells", "DC_Macro", "Plasma", "Tcells_mem", "Tcells_naive")
data_wide$Bcells_percent <- as.numeric(data_wide$Bcells)/(as.numeric(data_wide$Bcells) + as.numeric(data_wide$DC_Macro) + as.numeric(data_wide$Tcells_naive) + as.numeric(data_wide$Tcells_mem) + as.numeric(data_wide$Plasma)) *100
data_wide$DC_Macro_percent <- as.numeric(data_wide$DC_Macro)/(as.numeric(data_wide$Bcells) + as.numeric(data_wide$DC_Macro) + as.numeric(data_wide$Tcells_naive) + as.numeric(data_wide$Tcells_mem)+ as.numeric(data_wide$Plasma)) *100
data_wide$Tcells_mem_percent <- as.numeric(data_wide$Tcells_mem)/(as.numeric(data_wide$Bcells) + as.numeric(data_wide$DC_Macro) + as.numeric(data_wide$Tcells_naive) + as.numeric(data_wide$Tcells_mem) + as.numeric(data_wide$Plasma)) *100
data_wide$Tcells_naive_percent <- as.numeric(data_wide$Tcells_naive)/(as.numeric(data_wide$Bcells) + as.numeric(data_wide$DC_Macro) + as.numeric(data_wide$Tcells_naive) + as.numeric(data_wide$Tcells_mem) + as.numeric(data_wide$Plasma)) *100
data_wide$Plasma_percent <- as.numeric(data_wide$Plasma)/(as.numeric(data_wide$Bcells) + as.numeric(data_wide$DC_Macro) + as.numeric(data_wide$Tcells_naive) + as.numeric(data_wide$Tcells_mem) + as.numeric(data_wide$Plasma)) *100
data_wide_new <- subset(data_wide, select = c(Rep_age, Time, Bcells_percent, DC_Macro_percent, Plasma_percent, Tcells_mem_percent, Tcells_naive_percent))

#Calculate means for both groups
melted <- melt(data_wide_new, id.vars=c("Rep_age", "Time"))
means <- ddply(melted, c("variable", "Time"), summarise,
           mean=mean(value))

means$variable <- factor(means$variable, levels =c("DC_Macro_percent", "Tcells_naive_percent", "Tcells_mem_percent", "Bcells_percent", "Plasma_percent"))
means$Time <- factor(means$Time, levels = c("M3", "M18"))

plot <- ggplot(data=means, aes(x=Time, y=mean, fill=variable)) + scale_fill_manual(values=c("#c00000", "#FF9933", "#FFCC33", "#746CB1", "#BD77B2")) + geom_bar(stat="identity", width = 0.4, colour = "black") + scale_y_continuous(breaks = seq(0, 100, by = 25), expand = c(0,0)) + 
  xlab(" ") + ylab("Percentage (%)") + 
  theme_classic(base_size = 16, base_family = "Helvetica") + 
  theme(axis.text.y=element_text(size=16, face="bold")) + 
  theme(axis.title.y=element_text(size=16, face="bold", vjust=1)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1, size=16, face="bold")) +
  theme(legend.position="right")

# Calc SEM  
means.sem <- ddply(melted, c("variable", "Time"), summarise,
                   mean=mean(value), sem=sd(value)/sqrt(length(value)))
means.sem <- transform(means.sem, lower=mean, upper=mean+sem)
means.sem[means.sem$variable=='DC_Macro_percent',5:6] <- means.sem[means.sem$variable=='Plasma_percent',3] + means.sem[means.sem$variable=='DC_Macro_percent',5:6] + means.sem[means.sem$variable=='Tcells_naive_percent',5:6] + means.sem[means.sem$variable=='Tcells_mem_percent',5:6] + means.sem[means.sem$variable=='Bcells_percent',5:6]

means.sem[means.sem$variable=='Tcells_naive_percent',5:6] <- means.sem[means.sem$variable=='Plasma_percent',3] + means.sem[means.sem$variable=='Tcells_naive_percent',5:6] + means.sem[means.sem$variable=='Tcells_mem_percent', 5:6] + means.sem[means.sem$variable=='Bcells_percent', 5:6]

means.sem[means.sem$variable=='Tcells_mem_percent',5:6] <- means.sem[means.sem$variable=='Plasma_percent',3] + means.sem[means.sem$variable=='Tcells_mem_percent',5:6] + means.sem[means.sem$variable=='Bcells_percent', 5:6]

# Add SEM & change appearance of barplot
plotSEM <- plot + geom_errorbar(data=means.sem, aes(ymax=upper,  ymin=lower), 
                                width=0.15)
pdf("./Barplots_cell_composition/Plot_immune_SEM.pdf")
plotSEM
dev.off()

###################--------------------------#########################
plot <- ggplot(data=means, aes(x=Time, y=mean, fill=variable)) + scale_fill_manual(values=c("#c00000", "#FF9933", "#FFCC33", "#746CB1", "#BD77B2")) + geom_bar(stat="identity", width = 0.4, colour = "black") + scale_y_continuous(breaks = seq(0, 100, by = 25), expand = c(0,0)) +
  xlab(" ") + ylab("Percentage (%)") + 
  theme_classic(base_size = 16, base_family = "Helvetica") + 
  theme(axis.text.y=element_text(size=16, face="bold")) + 
  theme(axis.title.y=element_text(size=16, face="bold", vjust=1)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1, size=16, face="bold")) +
  theme(legend.position="right")

# Calc SEM  
means.sem <- ddply(melted, c("variable", "Time"), summarise,
                   mean=mean(value), sem=sd(value)/sqrt(length(value)))
means.sem <- transform(means.sem, lower=mean-sem, upper=mean+sem)
means.sem[means.sem$variable=='DC_Macro_percent',5:6] <- means.sem[means.sem$variable=='Plasma_percent',3] + means.sem[means.sem$variable=='DC_Macro_percent',5:6] + means.sem[means.sem$variable=='Tcells_naive_percent',5:6] + means.sem[means.sem$variable=='Tcells_mem_percent', 5:6] + means.sem[means.sem$variable=='Bcells_percent', 5:6]

means.sem[means.sem$variable=='Tcells_naive_percent',5:6] <- means.sem[means.sem$variable=='Plasma_percent',3] + means.sem[means.sem$variable=='Tcells_naive_percent',5:6] + means.sem[means.sem$variable=='Tcells_mem_percent', 5:6] + means.sem[means.sem$variable=='Bcells_percent', 5:6]

means.sem[means.sem$variable=='Tcells_mem_percent',5:6] <- means.sem[means.sem$variable=='Plasma_percent',3] + means.sem[means.sem$variable=='Tcells_mem_percent', 5:6] + means.sem[means.sem$variable=='Bcells_percent', 5:6]

# Add SEM & change appearance of barplot
plotSEM <- plot + geom_errorbar(data=means.sem, aes(ymax=upper,  ymin=lower), 
                                width=0.15)
pdf("./Barplots_cell_composition/Plot_immune_SEM_both.pdf")
plotSEM
dev.off()


