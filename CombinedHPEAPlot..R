setwd("~/path/to/output")
# Required libraries
library(dplyr)
library(purrr)
library(cinaR)
library(cinaRgenesets)
library(tidyr)
library(readxl)
library(stringr)
library(ggplot2)

#Load HPEA
source("~/Desktop/Mice_BC/DiffPeaks/HPEA.R")

## Load both immune modules and Wikipathways output files.
ImmMod <- read.csv("vp2008_hpea-list.csv")
WP <- read.csv("./wp_hpea-list.csv")

## Add Module column
ImmMod$Module <- "ImmuneModules"
WP$Module <- "Wikipathways"


ImmMod <- ImmMod[grep("Inflammation", ImmMod$module.name), ]
WP <- WP[grep("WP530|WP2703|WP205|WP1984|WP1836|WP585|WP1839|WP382|WP2775|WP75|WP585|WP615|WP619",WP$module.name),]
ImmMod <- ImmMod[(ImmMod$adj.p <= 0.1),]
WP <- WP[(WP$adj.p <= 0.1),]
# A <- as.vector(ImmMod$overlapping.genes)
# 
# ImmMod <- ImmMod[grep("Inflammation|U_P53 signaling|Cytotoxic cells", ImmMod$module.name), ]
# WP <- WP[grep("WP437|WP2704|WP382|WP69|WP585|WP49|WP205|WP3404|WP366|WP2112|WP195|WP75",WP$module.name),]
HPEA.table.filtered <- rbind(ImmMod, WP)
df.plot <- HPEA.table.filtered


filter.pathways <- TRUE
fdr.cutoff <- 0.1

if(filter.pathways){
  if (sum(df.plot$adj.p < fdr.cutoff) == 0){
    stop("You can't filter because there are no pathways to be displayed!")
  }
  df.plot <- subset(df.plot, adj.p < fdr.cutoff)
}

df.plot$contrast <- gsub("IL4I1_TRAF1_FLT3_CCR7_CCL22","mregDC", df.plot$contrast)

df.plot$Status <- factor(df.plot$Status, levels=c("Up", "Down"))
#df.plot$contrast <- factor(df.plot$contrast, levels=c("CD4_CCR7_JUN", "NaiveCD4_LEF1","MemoryCD4_CD44","CD8_CCR7","CD8_CCR9","CD8_GZMK","CD8_GZMM","CD8_CD69_ISG15","IL7R_ICOS","NK"))
df.plot$contrast <- factor(df.plot$contrast, levels=c("ActivatedMonocytes", "M1_Macrophages","M2_Macrophages","mregDC","FCGR1_NAPSA","DC","cDC1"))

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


pdf("./Combined.pdf", height = 10, width = 14)
plot.dot
dev.off()
