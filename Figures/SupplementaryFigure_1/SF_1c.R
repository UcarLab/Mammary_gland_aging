library(Seurat)
library(ggplot2)
mydata_new <- readRDS("/path/to/object")

mydata_new_rmv_doublet <- subset(mydata_new, idents = c("Doublet"), invert = TRUE)
marker_list_subset <- c("Pecam1", "Eng", "Cdh5", "Rgs5", "Notch3", "Des", "Col1a1", "Fn1", "Jchain", "Mzb1", "Cd79a", "Blnk", "Cd44", "Il7r", "Sell", "Cd8a", "Cd8b1", "Cd4", "Cd3d", "Itgam", "Csf1r", "Cd86", "Cd163", "Fcgr2b", "Itgax", "Aif1","Myl9", "Acta2", "Krt17", "Prlr", "Cited1", "Esr1",  "Mfge8", "Trf", "Csn3")
mydata_new_rmv_doublet@active.ident <- factor(mydata_new_rmv_doublet@active.ident,levels=c("Luminal-AV", "Luminal-HS", "Myoepithelial", "Dendritic/Macrophages", "Tcells_naive", "Tcells_mem", "Bcells", "Plasma", "Fibroblasts", "Pericytes", "Vascular"))

pdf("path/to/pdf", width=12, height=18)
print(DotPlot(object = mydata_new_rmv_doublet, features = marker_list_subset, cols = c("lightgrey", "red"), dot.scale=10) + coord_flip() + theme(legend.position = "bottom", legend.box = "vertical", legend.key.width = unit(1, "cm")) + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) )
dev.off()