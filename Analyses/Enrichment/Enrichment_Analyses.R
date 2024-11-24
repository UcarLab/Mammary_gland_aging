setwd("~/path/to/output")
# Hypergeometric P-value Enrichment Analyses

# Required libraries
library(dplyr)
library(purrr)
library(cinaR)
library(cinaRgenesets)
library(ggplot2)
library(plyr)
#Load HPEA
source("~/path/to/HPEA.R")

## The following steps are for running the Wikipathways/ImmuneModules dataset. Please comment out everything until Fold change and p-value cutoffs and change geneset name accordingly 
## Change name to vp2008 to run for for Immune Modules.
## Change name to wp to run for for Wikipathways.
load("~/path/to/enrichment_analysis.Rdata")
geneset.name <- "wp"

# Selected modules 
# There are multiple genesets consisting of multiple modules.
# We select one geneset per test and correct the p-values for it.
geneset <- selected_genesets[[geneset.name]]

# Module annotations
geneset.label <- selected_genesets_labels[[geneset.name]]

# Merge annotations via modules
geneset.merged <- merge(geneset, geneset.label, by = "Module.ID")

# Fold change and p-value cutoffs
cutoff.fc <- 0
cutoff.p  <- 0.1
cutoff.tss<- 50e3

# Contrast file names
filelist = list.files(pattern = "*txt") 
filelist 

## Check contrast list to ensure the vector of names for plot headers is acceptable.
contrast.list <- (filelist) %>% strsplit("_", fixed = T) %>% 
  sapply(function(x){paste(x[1:2], collapse = "_")}) %>% 
  strsplit(".", fixed = T) %>% sapply(function(x){x[1]})
contrast.list

## This function reads in all the text files present in the given directory. Please save the text files accordingly in one place and name them by the conditions the differential analysis was run in.
## Update the contrast.list function above if the filenames are saved differently.

HPEA.list <- lapply(filelist, function(file){
  
  contrast <- file %>% 
    strsplit("_", fixed = T) %>% 
    sapply(function(x){paste(x[1:2], collapse = "_")})
  
  # read the DE results
  peaks <- read.csv(file, sep = "\t")
  peaks$geneSymbol <- peaks$gene_name
  # filter out the genes bigger than tss cutoff
  peaks.filtered <- peaks[abs(peaks$distanceToTSS) < cutoff.tss,]
  
  # order according to tss
  peaks.filtered <- peaks.filtered %>% arrange(abs(distanceToTSS))
  
  # remove duplicated gene names
  peaks.filtered <- peaks.filtered[!duplicated(peaks.filtered$geneSymbol),]
  
  #peaks.filtered$logFC <- peaks.filtered$avg_logFC
  # For RNA-seq analyses you need up and down regulated genes (and opening/closing peaks for ATACseq).
  # The sign of this regulation is determined via fold-change (FC) or log(FC) of the genes, 
  # which are calculated with differential analyses. Then, these genes are filtered further
  # with respect to their p-values (or adjusted p-values). 
  
  # DE genes up-regulated
  DE.genes.ur <- peaks.filtered[peaks.filtered$logFC > cutoff.fc & peaks.filtered$FDR < cutoff.p,]
  # DE genes down-regulated
  DE.genes.dr <- peaks.filtered[peaks.filtered$logFC < cutoff.fc & peaks.filtered$FDR < cutoff.p,]
  
  # If both up and down regulated genes do not exist skip the contrast
  if(nrow(DE.genes.ur) > 0 | nrow(DE.genes.dr) > 0){
    
    gene.names.ur <- DE.genes.ur$geneSymbol
    gene.names.dr <- DE.genes.dr$geneSymbol
    
    result.ur <- HPEA(geneset.merged, gene.names.ur, geneset.name = geneset.name, contrast = contrast)
    result.dr <- HPEA(geneset.merged, gene.names.dr, geneset.name = geneset.name, contrast = contrast)
    
    return(list(UR = result.ur, DR = result.dr))
  }
  return(NULL)
  
})

names(HPEA.list) <- contrast.list

# Merge UR/DR lists
HPEA.list <- lapply(HPEA.list, function(x){
  if (!is.null(x)){
    rbind(
      cbind(x[["UR"]], Status = "Up"), 
      cbind(x[["DR"]], Status = "Down")
    )
  }
})

# Make it a table
HPEA.table <- map_df(HPEA.list, ~as.data.frame(.x), .id="contrast")

# Remove empty rows
HPEA.table.filtered <- HPEA.table[HPEA.table$p < 1,]
write.csv(HPEA.table.filtered, file = paste0("./", geneset.name, "_hpea-list.csv"))
df.plot <- HPEA.table.filtered


filter.pathways <- TRUE
fdr.cutoff <- 0.1

if(filter.pathways){
  if (sum(df.plot$adj.p < fdr.cutoff) == 0){
    stop("You can't filter because there are no pathways to be displayed!")
  }
  df.plot <- subset(df.plot, adj.p < fdr.cutoff)
}

df.plot$Status <- factor(df.plot$Status, levels=c("Up", "Down"))

## Set order for contrast level for Dendritic Cells/Macrophage groups and/or TCells.
#df.plot$contrast <- factor(df.plot$contrast, levels=c("Naive CD4-1", "Naive CD4-2","Memory CD4","Naive CD8-1","Naive CD8-2","CD8 GZMK","CD8 GZMM","CD8 ISG15","γδ T","NK"))
#df.plot$contrast <- factor(df.plot$contrast, levels=c("DC_C1", "DC_C2","DC_C3","DC_C4","DC_C5","DC_C6","DC_C7"))

# create ggplot
color_values <- color_values[c(4,7)]
plot.dot <- ggplot2::ggplot(df.plot,
                            ggplot2::aes(x = contrast,
                                         y = module.name,
                                         size = ifelse(adj.p < fdr.cutoff, -log(adj.p), NA),
                                         color = Status))

plot.dot <- plot.dot + ggplot2::geom_point()

plot.dot <- plot.dot + facet_grid(scales='free', space = "free")
plot.dot <- plot.dot + ggplot2::labs(x = "Contrast",
                                     y = "Pathways",
                                     color = "Sign",
                                     size = "-log10(adj.p)",
                                     caption = paste0("FDR < ", fdr.cutoff))
plot.dot <- plot.dot + ggplot2::scale_color_manual(values = color_values)
plot.dot <- plot.dot + theme(axis.text =element_text(size=13),
                             axis.title=element_text(size=15,face="bold"))
plot.dot <- plot.dot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot.dot <- plot.dot + theme(legend.text = element_text(size=16))
plot.dot <- plot.dot + facet_grid(Module ~ ., space = "free", scales = "free")


pdf("~/path/to/Dotplot.pdf", height = 10, width = 14)
plot.dot
dev.off()
