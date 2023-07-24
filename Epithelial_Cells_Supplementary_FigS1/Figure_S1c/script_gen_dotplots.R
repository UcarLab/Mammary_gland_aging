#Script to generate DotPlots for marker genes.

library(Seurat)
library(ggplot2)

dir.create("./DotPlots")
gen_dotplots <- function(mydata_new){
	mydata_new_rmv_doublet <- subset(mydata_new, idents = c("Doublet"), invert = TRUE)
	marker_list_orig <- c("Pecam1", "Eng", "Cdh5", "Rgs5", "Notch3", "Des", "Col1a1", "Fn1", "Ighg1", "Ighg3", "Cd27", "Jchain", "Mzb1", "Gzmb", "Ms4a2", "Ighd", "Blnk","Cd79a", "Cd8a", "Cd8b1", "Cd3d", "Cd4", "Sell", "Cd44", "Itgam", "Csf1r", "Cd86", "Cd163", "Ly75", "Fcgr2b", "Itgax", "Aif1","Myl9", "Acta2", "Krt17", "Prlr", "Cited1", "Esr1",  "Mfge8", "Trf", "Csn3", "Pi16", "Clec4a", "Cd38", "IgG", "Tnfrsf13b", "Tnfrsf17", "Il7r", "Ccr7", "Cd19", "Ms4a1", "Cd34", "B3gat1", "Klrg1", "Eomes", "Prdm1", "Ptprc", "Id2", "Stat4", "T-bet", "Tcf1", "Bcl6", "Id3", "Stat3", "Cx3cr1")
	
	marker_list_subset <- c("Pecam1", "Eng", "Cdh5", "Rgs5", "Notch3", "Des", "Col1a1", "Fn1", "Jchain", "Mzb1", "Cd79a", "Blnk", "Cd44", "Il7r", "Sell", "Cd8a", "Cd8b1", "Cd4", "Cd3d", "Itgam", "Csf1r", "Cd86", "Cd163", "Fcgr2b", "Itgax", "Aif1","Myl9", "Acta2", "Krt17", "Prlr", "Cited1", "Esr1",  "Mfge8", "Trf", "Csn3")

	mydata_new_rmv_doublet@active.ident <- factor(mydata_new_rmv_doublet@active.ident,
        	                    levels=c("Luminal-AV", "Luminal-HS", "Myoepithelial", "Dendritic/Macrophages", "Tcells_naive", "Tcells_mem", "Bcells", "Plasma", "Fibroblasts", "Pericytes", "Vascular"))

	pdf("../DotPlots/DotPlot_markers_all_celltypes_flip_new.pdf", width=12, height=18)
	print(DotPlot(object = mydata_new_rmv_doublet, features = marker_list_subset, cols = c("lightgrey", "red"), dot.scale=10) + coord_flip() + theme(legend.position = "bottom", legend.box = "vertical", legend.key.width = unit(1, "cm")) + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) )
	dev.off()

}

####################################
dir.create("DotPlots")
args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else
{
print("Generating dotplots using marker genes...")
}

mydata_new = readRDS(args[1])

gen_dotplots(mydata_new)

