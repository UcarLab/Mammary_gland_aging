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

## Fig7a
## Load both immune modules and Wikipathways output files. ## Provided in analyses files. Can also be loaded in using cinarGenesets.
ImmMod <- read.csv("~/path/to/vp2008_hpea-list.csv")
WP <- read.csv("~/path/to/wp_hpea-list.csv")

## Add Module column
ImmMod$Module <- "ImmuneModules"

## Selecting Modules
ImmMod <- ImmMod[grep("Inflammation|Cytotoxic Cells", ImmMod$module.name), ]
## These modules were picked manually. Were renamed to actual pathways while editing the panels later. 
ImmMod <- ImmMod[(ImmMod$adj.p <= 0.1),]
df.plot <- ImmMod

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
df.plot$contrast <- factor(df.plot$contrast, levels=c("Memory CD4", "CD8 GZMK+","CD8 GZMM+","CD8 ISG15+","GD and MAIT"))

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


pdf("~/path/to/pdf", height = 10, width = 14)
plot.dot
dev.off()
