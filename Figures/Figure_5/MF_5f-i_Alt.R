setwd("~/Desktop/Mice_BC/FinalFigurePanels")
#Script to generate expression plots using seurat object.
library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(gridBase)
library(gridExtra)

mydata = readRDS("~/path/to/Seurat/Object")
gene_name = #Gene of Choice
Idents(mydata) <- #Annotation Label
table(Idents(mydata))
mydata$celltype.stim <- paste(Idents(mydata), mydata$orig.ident, sep = "_") #orig.ident should be age labels
mydata$celltype <- Idents(mydata)
Idents(mydata) <- "celltype.stim"
table(Idents(mydata))

## Example shown here for Luminal AV. ## Redo for each cell type and gene. Change colors and gene and cell type above.
mydata_Lum_3M <- subset(mydata, idents = c("Luminal_AV_mm10_3mths"))
mydata_Lum_18M <- subset(mydata, idents = c("Luminal_AV_mm10_18mths"))
mydata_new <- subset(mydata, subset = Lag3 > 0)
mydata_Lum_3M_no_zero <- subset(mydata_new, idents = c("Luminal_AV_mm10_3mths"))
mydata_Lum_18M_no_zero <- subset(mydata_new, idents = c("Luminal_AV_mm10_18mths"))

subset_LumAV <- subset(mydata_new, idents = c("Luminal_AV_mm10_3mths", "Luminal_AV_mm10_18mths"))
#my_levels <- c("Luminal-AV_mm10_3mths", "Luminal-AV_mm10_18mths")
#subset_LumAV@meta.data$orig.ident <- factor(x = subset_LumAV@meta.data$orig.ident, levels = my_levels)

# exprn_val <- subset(subset_LumAV[["RNA"]]@data, rownames(subset_LumAV[["RNA"]]@data) %in% gene_name)
exprn_val <- FetchData(subset_LumAV, gene_name)
exprn_val_df <- as.data.frame(exprn_val)
rownames(exprn_val_df) <- colnames(subset_LumAV[["RNA"]]@data)
celltype_age <- subset_LumAV@meta.data$celltype.stim
celltype <- subset_LumAV@meta.data$celltype
#exprn_df_new <- as.data.frame(t(exprn_val_df))
df_final_AV <- cbind(exprn_val_df, celltype_age, celltype)
names(df_final_AV) <- c("exprn", "celltype_age", "celltype")

Celltype = "Luminal_AV"
var_18m <- paste0(Celltype, "_mm10_18mths")
var_3m <- paste0(Celltype, "_mm10_3mths")

df_final_AV_18m <- df_final_AV[df_final_AV$celltype_age == var_18m,]
df_final_AV_3m <- df_final_AV[df_final_AV$celltype_age == var_3m,]
test_sig <- t.test(df_final_AV_18m$exprn, df_final_AV_3m$exprn)
print(Celltype)
print(test_sig)
print(test_sig$p.value)

#######---------------------############
#df_table <- table(Idents(subset_LumAV))
m1 = ncol(mydata_Lum_3M_no_zero)
m2 = ncol(mydata_Lum_18M_no_zero)
m1_remain = ncol(mydata_Lum_3M) - m1
m2_remain = ncol(mydata_Lum_18M) - m2

df_3M_AV <- data.frame("Cells" = c("Expressed", "Not expressed"), "Value" = c(m1,m1_remain))
df_18M_AV <- data.frame("Cells" = c("Expressed", "Not expressed"), "Value" = c(m2,m2_remain))


## Need this if combining multiple cell types together

# df_combined <- rbind(df_final_AV, df_final_HS, df_final_Myo, df_final_Myo1)
# level_order <- c("Luminal-HS_mm10_3mths", "Luminal-HS_mm10_18mths", "Luminal-AV_mm10_3mths", "Luminal-AV_mm10_18mths", "Myoepithelial_mm10_3mths", "Myoepithelial_mm10_18mths", "Fibroblasts_mm10_3mths", "Fibroblasts_mm10_18mths")
# #df_combined <- rbind(df_final_HS)
# #level_order <- c("Luminal-AV_mm10_3mths", "Luminal-AV_mm10_18mths")
# #scale_fill_manual(values=c("firebrick3", "firebrick4", "purple", "plum3", "dodgerblue1", "dodgerblue3", "goldenrod2", "goldenrod4")) + geom_jitter(shape=16, position=position_jitter(0.2),cex=5) +

## Modify colors based on celltype used. Switch out df_final_AV with df_combined from above if plotting multiple cell types together.
colors_celltype <- c("firebrick3", "firebrick4")
library(cowplot)
pdf("~/path/to/fig", width=88, height=40)
vlnplot_n <- df_final_AV %>%
	ggplot(aes(y=exprn, x = factor(celltype_age, level = level_order), fill = factor(celltype_age, level = level_order))) +
        geom_violin(position = position_dodge(0.8), alpha=0.5) +
  scale_fill_manual(values=colors_celltype) + geom_jitter(shape=16, position=position_jitter(0.2),cex=5) +
        xlab("Sample")+ylab("Expression")+
        theme(plot.margin = unit(c(4,4,20,10), "cm"))+ggtitle(gene_name)+theme(plot.title = element_text(hjust = 0.5, size=80))+ theme(axis.text.x=element_text(size=60, angle=90), axis.text.y = element_text(size = 60), axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60), legend.title=element_text(size=60), legend.text=element_text(size=60))+theme(legend.key.size = unit(2, "cm")) + theme(panel.background = element_blank(), axis.line = element_line(color = "black"))
+ theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank(), axis.line = element_line(color = "black"))
bp1<- ggplot(df_3M_AV, aes(x="", y=Value, fill=Cells)) + geom_bar(width = 0.2, stat = "identity") + theme_void()+theme(legend.position = "none")
pie1 <- bp1 + coord_polar("y", start=0) + scale_fill_grey()
bp2<- ggplot(df_18M_AV, aes(x="", y=Value, fill=Cells)) + geom_bar(width = 0.2, stat = "identity") + theme_void()+ theme(legend.position = "none")
pie2 <- bp2 + coord_polar("y", start=0) + scale_fill_grey()
ggdraw() +draw_plot(vlnplot_n, x=0, y=0, width=1, height=1)+draw_plot(pie1, x=0.04, y=-0.001, width=0.2, height=0.2)+draw_plot(pie2, x=0.15, y=-0.001, width=0.2, height=0.2)
dev.off()