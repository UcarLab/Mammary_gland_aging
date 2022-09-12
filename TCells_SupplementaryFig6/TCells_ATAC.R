setwd("~/Desktop/Mice_BC/FinalFigurePanels")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape)
library(readxl)
library(purrr)
library(cinaR)
library(cinaRgenesets)
library(tidyr)
library(stringr)
library(scater)
library(Seurat)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

TCells <- readRDS("./TCells_ATAC.rds")
Idents(TCells) <- "predicted.id"
#TCells <- RunUMAP(TCells, dims = 1:30, reduction = "har30", reduction.name="humap33", reduction.key = "humap3key")
TCells <- subset(TCells, idents = "Doublets", invert = T)
my_levels <- c("Naive CD4 LEF1+","CD4 CCR7+ JUN+","Memory CD4 CD44+", "CD8 CCR7+ ", "CD8 CCR9+", "CD8 GZMK+", "CD8 GZMM +", "IL7R+ ICOS+","NK Cells")
#my_levels <- rev(my_levels)
TCells@meta.data$predicted.id <- factor(x = TCells@meta.data$predicted.id, levels = my_levels)

##FIGURE 1A
#colr <- c("orangered1","deepskyblue","lightblue2","darkblue","purple","dodgerblue3", "goldenrod2", "firebrick4", "lightpink", "springgreen4")
colr <- c("springgreen4","goldenrod2","dodgerblue3","purple","lightblue2", "deepskyblue", "firebrick4","orangered1", "lightpink")
p1 <- DimPlot(TCells, reduction = "umap3", group.by = "predicted.id", cols = rev(colr),label = F, label.size = 8)
pdf("./TCells_ATAC_Fig1A.pdf", height = 6, width = 8)
p1
dev.off()

##FIGURE 1B
Idents(TCells) <- "age"
OldTCells <- subset(TCells, idents = "18M")
YoungTCells <- subset(TCells, idents = "3M")
p1 <- DimPlot(OldTCells, reduction = "umap3", group.by = "predicted.id", cols = rev(colr),label = F, label.size = 8)
p2 <- DimPlot(YoungTCells, reduction = "umap3", group.by = "predicted.id", cols = rev(colr),label = F, label.size = 8)
pdf("./TCells_ATAC_Fig1B.pdf", height = 4, width = 6)
p1
p2
dev.off()

##FIGURE 1C
# cd_genes <- list('TCells'=c('Cd3e','Cd4','Cd8a'),
#                  'Naive T'=c('Sell','Il7r','Lef1', 'Ccr7', "Ccr9"),
#                  'Memory T'=c('Cd44', "Itgae"), 
#                  'Cytotoxic T'=c( 'Gzma', 'Gzmk', 'Gzmm'),
#                  'Activated T'=c('Cd69',  'Icos', 'Jun'),
#                  'Interferon'=c('Isg15',"Irf1","Irf7"), 
#                  'Gamma Delta'=c("Il17a",'Trdc',"Trgv2"),
#                  'NK'=c('Ncr1','Nkg7', "Klrc1", "Xcl1", "Ctsw", "Car2"))
# p1 <- DotPlot(TCells,features = cd_genes, cols = c("lightgrey","red"), group.by = "NewAnnotations", dot.scale = 4, cluster.idents = F) + RotatedAxis() + theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
# pdf("./TCells_RNA_Fig1C.pdf", height = 8, width = 16)
# p1
# dev.off()

##FIGURE 1C
Meta <- TCells@meta.data
Meta <- Meta[,c(22,23,36)]
Meta <- Meta %>% group_by_all() %>% summarise(COUNT = n())
#Meta$replicate <- gsub("M3","M18",Meta$replicate)
NewMeta = cast(Meta, rep + predicted.id ~ age, sum) 
NewMeta$Ratio <- log2(NewMeta$`18M`/NewMeta$`3M`)
NewMeta <- NewMeta[,c(1,2,5)]
NewMeta$predicted.id <- as.character(NewMeta$predicted.id)
NewMeta[nrow(NewMeta) + 1,] <- c("rep4","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("rep5","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("rep6","Base",0)
NewMeta$predicted.id <- as.factor(NewMeta$predicted.id)
NewMeta$Ratio <- as.numeric(NewMeta$Ratio)
colr <- c("black","springgreen4","goldenrod2","dodgerblue3","purple","lightblue2", "deepskyblue", "firebrick4","orangered1", "lightpink")

p1 <- ggplot(NewMeta,aes(x = reorder(predicted.id, -Ratio), y = Ratio, fill = predicted.id)) + 
  geom_boxplot()+
  scale_fill_manual(values=colr) +
  xlab("CellTypes") +
  ylab("log2 FC (Old/Young)") +
  ggtitle("Cell Type Proportions") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  stat_compare_means(method = "anova", label.y = 10)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = "Base", hide.ns = T, label.y = 8)   

stat1 <- compare_means(Ratio ~ predicted.id,  data = NewMeta, ref.group = "Base", method = "t.test")
write.table(stat1, "./ATAC_FoldChange_Stats_TCells.txt", row.names = F, sep = '\t')


p1                                                                            
pdf("./TCells_ATAC_Fig1C.pdf", height = 6, width = 8)
p1
dev.off()

# ##FIGURE 1E
DefaultAssay(TCells_RNA) <- "RNA"
p1 <- FeaturePlot(TCells, features = c("Gzmk","Pdcd1","Ctla4", "Lag3", "Ccl5","Itgae"), reduction = "umap3", cols = c("lightgrey", "red3"),label = F, ncol = 3, keep.scale="all", slot = "data")
pdf("./TCells_ATAC_Fig1F.pdf", height = 6, width = 12)
p1
dev.off()

##FIGURE 1F
A <- read_xlsx("~/Desktop/Mice_BC/Subclustering_TCells/Counts1.xlsx")
A <- A[-1,]
A <- A[,c(1,4,5)]
colnames(A) <- c("CellType", "Opening", "Closing")
A <- A[-c(4,9),]

Closing <- A[,c(1,2)]
Opening <- A[,c(1,3)]
colnames(Opening) <- c("CellType","Counts")
colnames(Closing) <- c("CellType","Counts")
Opening$Direction <- "Opening"
Closing$Direction <- "Closing"
Closing$Counts <- as.numeric(Closing$Counts) * -1
FinalCounts <- rbind(Opening, Closing)
long_DF <- FinalCounts %>% gather("Type", "Count", Counts)
my_levels <- c("NaiveCD4_LEF1","CD4_CCR7_JUN","MemoryCD4_CD44", "CD8_CCR7", "CD8_CCR9", "CD8_GZMK", "CD8_GZMM", "IL7R_ICOS")
long_DF$CellType <- factor(long_DF$CellType, levels=my_levels)
long_DF$Count <- as.numeric(long_DF$Count)

pdf("./TCells_ATAC_Fig1D.pdf", height = 9, width = 14)
ggplot(long_DF[order(long_DF$Type, decreasing = F),], aes(y=Count, x=CellType, fill = Direction)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("lightgrey","black")) +
  xlab("Gene Counts") +
  theme(axis.text =element_text(size=16),axis.title=element_text(size=16,face="bold")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Cell Types")+
  guides(fill=guide_legend(title="Direction"))
dev.off()

##FIGURE 1E
ImmMod <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/ATAC/DiffPeaks_Type/Chipseeker_Annots/vp2008_hpea-list.csv")
WP <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/ATAC/DiffPeaks_Type/Chipseeker_Annots//wp_hpea-list.csv")
ImmMod$Module <- "ImmuneModules"
WP$Module <- "Wikipathways"

ImmMod <- ImmMod[grep("Inflammation|Cytotoxic cells", ImmMod$module.name), ]
WP <- WP[grep("WP437|WP2704|WP382|WP69|WP585|WP49|WP205)|WP3404|WP366|WP2112|WP75",WP$module.name),]
ImmMod <- ImmMod[(ImmMod$adj.p <= 0.1),]
WP <- WP[(WP$adj.p <= 0.1),]
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

df.plot$Status <- factor(df.plot$Status, levels=c("Up", "Down"))
df.plot$contrast <- factor(df.plot$contrast, levels=c("NaiveCD4_LEF1","CD4_CCR7_JUN","MemoryCD4_CD44","CD8_CCR7","CD8_CCR9","CD8_GZMK","CD8_GZMM","IL7R_ICOS","NK"))

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
pdf("./TCells_ATAC_Fig1E.pdf", height = 10, width = 12)
plot.dot
dev.off()


##Figure 1G - Motifs
Motifs <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/ATAC/DiffPeaks_Type/Chromvar/FinalMotifs1.txt", sep = '\t', header = T)
colnames(Motifs) <- c("ID", "p_val", "log_FC_Change", "pct.1", "pct.2", "adj_p_val", "CellType")
Motifs$CellType <- gsub(".txt","",Motifs$CellType)
Motifs$log_FC_Change <- as.numeric(Motifs$log_FC_Change)
MotifNames <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/ATAC/DiffPeaks/Motifs/Final_Total_Motifs1.txt", sep = '\t', header = F)
MotifNames <- MotifNames[-1,]
MotifNames <- MotifNames[,c(2,3)]
colnames(MotifNames) <- c("ID", "MotifName")
MotifNames <- unique(MotifNames)
Motifs <- merge(MotifNames, Motifs)
Motifs2 <- Motifs
Motifs <- Motifs[(Motifs$log_FC_Change > 0),]

Motifs <- Motifs[(Motifs$log_FC_Change >= 1),]
Motifs <- Motifs[(Motifs$adj_p_val <= 0.05),]

Motifs1 <- Motifs %>% group_by(CellType) %>% top_n(20, abs(log_FC_Change))
MotifsList <- as.data.frame(unique(Motifs1$MotifName))
colnames(MotifsList) <- "MotifName"

FinalMotifs <- merge(Motifs2, MotifsList)
FinalMotifs <- FinalMotifs[,c(1,4,7,8)]
FinalMotifs$Significance <- "*"

for (i in 1:nrow(FinalMotifs)) {
  if (FinalMotifs$adj_p_val[i] > 0.05) {
    FinalMotifs$Significance[i] = ""
  }
}
Clustering <- FinalMotifs[,c(1,2,4)]
Clustering1 <- cast(Clustering, MotifName ~ CellType, mean, value = 'log_FC_Change')
rownames(Clustering1) <- Clustering1$MotifName
Clustering1 <- Clustering1[,-1]
data <- scale(t(Clustering1))

ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
ord
FinalMotifs$CellType <- factor(FinalMotifs$CellType, levels = c("MemoryCD4_CD44","IL7R_ICOS","NK","CD4_CCR7_JUN","NaiveCD4_LEF1","CD8_GZMK","CD8_GZMM", "CD8_CCR7", "CD8_CCR9"))
RS <- FinalMotifs[,c(1,5)]
RS1 <- RS %>% group_by_all() %>% summarise(COUNT = n())
RS1 <- RS1[!(RS1$Significance=="" & RS1$COUNT==9),]
RS1 <- as.data.frame(unique(RS1$MotifName))
colnames(RS1) <- "MotifName"

FinalMotifs <- merge(RS1, FinalMotifs)

pdf("./TCells_ATAC_Fig1G_Poster.pdf", height = 12, width = 12)
ggplot(FinalMotifs, aes(CellType, MotifName, fill= log_FC_Change)) +
  geom_tile() +
  scale_fill_distiller(palette = "PuOr", ) +
  geom_text(aes(label=Significance), color="black", size=4) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

## FIGURE 1I - Motifs by Age
Motifs <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/ATAC/DiffPeaks/Motifs/Final_Total_Motifs1.txt", sep = '\t', header = F)
Motifs <- Motifs[-1,]
Motifs <- Motifs[,-1]
colnames(Motifs) <- c("ID","MotifName", "p_val", "log_FC_Change", "pct.1", "pct.2", "adj_p_val", "CellType")
Motifs$CellType <- gsub(".txt","",Motifs$CellType)
Motifs$log_FC_Change <- as.numeric(Motifs$log_FC_Change)
Motifs2 <- Motifs
Motifs <- Motifs[(Motifs$log_FC_Change >= 1 | Motifs$log_FC_Change <= -1),]
Motifs <- Motifs[(Motifs$adj_p_val <= 0.05),]

Motifs1 <- Motifs %>% group_by(CellType) %>% top_n(25, abs(log_FC_Change))
MotifsList <- as.data.frame(unique(Motifs1$MotifName))
colnames(MotifsList) <- "MotifName"

FinalMotifs <- merge(Motifs2, MotifsList)
FinalMotifs <- FinalMotifs[,c(1,4,7,8)]
FinalMotifs$Significance <- "*"

for (i in 1:nrow(FinalMotifs)) {
  if (FinalMotifs$adj_p_val[i] > 0.05) {
    FinalMotifs$Significance[i] = ""
  }
}
Clustering <- FinalMotifs[,c(1,2,4)]
Clustering1 <- cast(Clustering, MotifName ~ CellType, mean, value = 'log_FC_Change')
rownames(Clustering1) <- Clustering1$MotifName
Clustering1 <- Clustering1[,-1]
data <- scale(t(Clustering1))

ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
ord

FinalMotifs$CellType <- factor(FinalMotifs$CellType, levels = c("MemoryCD4_CD44","NaiveCD4_LEF1","CD4_CCR7_JUN","IL7R_ICOS","NK","CD8_CCR7","CD8_CCR9", "CD8_GZMK", "CD8_GZMM"))
RS <- FinalMotifs[,c(1,5)]
RS1 <- RS %>% group_by_all() %>% summarise(COUNT = n())
RS1 <- RS1[!(RS1$Significance=="" & RS1$COUNT==3),]
RS1 <- as.data.frame(unique(RS1$MotifName))
colnames(RS1) <- "MotifName"

FinalMotifs <- merge(RS1, FinalMotifs)

pdf("./TCells_ATAC_Fig1H.pdf", height = 12, width = 12)
ggplot(FinalMotifs, aes(CellType, MotifName, fill= log_FC_Change)) +
  geom_tile() +
  scale_fill_distiller(palette = "PuOr", ) +
  geom_text(aes(label=Significance), color="black", size=4) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()





##FIGURE 1H - CD8 GZMM
normalized_counts <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/Pseudobulk/Cluster6/DifferentialGeneCounts_All_Cluster6.txt", sep = '\t')
sig_genes <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/Pseudobulk/Combined/CD8_GZMM.txt", sep = '\t')

sig_genes_down <- c("Jun","Fos","Jund","Cxcr5","Ccr9","Irf4","Tcf7","Cd27","Lef1")
#sig_genes_down <- c("Ccr7","Ifi214","Klf2","Sell","Lef1","Foxp1","Bach2")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

sig_genes_up <- c("Gzmk","Gzmm","Gzmb","Prf1","Nkg7","S100a4","S100a6","S100a10","S100a13","Ccr2","Cxcr6","Ccr5","Ccl5","Il7r","Ctla2a","Cdk4","Fasl", "Eomes") 
#sig_genes_up <- c("Cxcl2","Cxcl3","Ccr2","Ccr1","Ccl5","Cxcr6","S100a4","S100a6","S100a11","Gzmk","Fosb","Itga1","Lmna","Ctla2a","Klrg1","Il7r","Casp1","Casp4","Lag3","Dusp1","Junb","Pdcd1","Anxa1","Cdk6")
sig_genes_up <- as.data.frame(sig_genes_up)
colnames(sig_genes_up) <- "genes"
sig_genes_up$Direction <- "Opening"

sig_genes <- rbind(sig_genes_down, sig_genes_up)
rownames(sig_genes) <- sig_genes$genes

sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_genes$genes)

sig_norm <- sig_norm[,c(1,8:13,2:7)]

# Run pheatmap using the metadata data frame for the annotation
heat_colors <- rev(brewer.pal(9, "RdBu"))
df <- colnames(sig_norm[,2:length(colnames(sig_norm))])
df <- as.data.frame(df)
colnames(df) <- "group_id"
df$group_id[1:6] <- "Young"
df$group_id[7:12] <- "Old"
rownames(df) <- colnames(sig_norm[,2:length(colnames(sig_norm))])

# df_row <- as.data.frame(sig_genes$Role)
# rownames(df_row) <- rownames(sig_genes)
# colnames(df_row) <- "Role"

rownames(sig_norm) <- sig_norm$gene
sig_norm <- sig_norm[,-1]
ann_colors = list(group_id = c(Young = "lightgray", Old = "black"))

pdf("./TCells_RNA_Fig1H.pdf", height = 6, width = 6)
pheatmap(sig_norm, 
         color = heat_colors, 
         annotation_colors = ann_colors,
         cluster_rows = T, 
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         treeheight_row = 0, treeheight_col = 0,
         annotation_col = df, 
         fontsize_row = 13,
         scale = "row")      
dev.off()

##FIGURE 1I - Memory CD4 CD44
normalized_counts <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/Pseudobulk/Cluster9/DifferentialGeneCounts_All_Cluster9.txt", sep = '\t')
sig_genes <- read.csv("~/Desktop/Mice_BC/Subclustering_TCells/Pseudobulk/Combined/MemoryCD4_CD44.txt", sep = '\t')

#sig_genes_down <- c("Jun","Fos","Jund","Cxcr5","Ccr9","Irf4","Tcf7","Cd27","Lef1")
sig_genes_down <- c("Ccr7","Ifi214","Klf2","Sell","Lef1","Foxp1","Bach2")
sig_genes_down <- as.data.frame(sig_genes_down)
colnames(sig_genes_down) <- "genes"
sig_genes_down$Direction <- "Closing"

#sig_genes_up <- c("Gzmk","Gzmm","Gzmb","Prf1","Nkg7","S100a4","S100a6","S100a10","S100a13","Ccr2","Cxcr6","Ccr5","Ccl5","Il7r","Ctla2a","Cdk4","Fasl", "Eomes") 
sig_genes_up <- c("Cxcl2","Cxcl3","Ccr2","Ccr1","Ccl5","Cxcr6","S100a4","S100a6","S100a11","Gzmk","Fosb","Itga1","Lmna","Ctla2a","Klrg1","Il7r","Casp1","Casp4","Lag3","Dusp1","Junb","Pdcd1","Anxa1","Cdk6")
sig_genes_up <- as.data.frame(sig_genes_up)
colnames(sig_genes_up) <- "genes"
sig_genes_up$Direction <- "Opening"

sig_genes <- rbind(sig_genes_down, sig_genes_up)
rownames(sig_genes) <- sig_genes$genes

sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_genes$genes)

sig_norm <- sig_norm[,c(1,8:13,2:7)]

# Run pheatmap using the metadata data frame for the annotation
heat_colors <- rev(brewer.pal(9, "RdBu"))
df <- colnames(sig_norm[,2:length(colnames(sig_norm))])
df <- as.data.frame(df)
colnames(df) <- "group_id"
df$group_id[1:6] <- "Young"
df$group_id[7:12] <- "Old"
rownames(df) <- colnames(sig_norm[,2:length(colnames(sig_norm))])

# df_row <- as.data.frame(sig_genes$Role)
# rownames(df_row) <- rownames(sig_genes)
# colnames(df_row) <- "Role"

rownames(sig_norm) <- sig_norm$gene
sig_norm <- sig_norm[,-1]
ann_colors = list(group_id = c(Young = "lightgray", Old = "black"))

pdf("./TCells_RNA_Fig1I.pdf", height = 7, width = 7)
pheatmap(sig_norm, 
         color = heat_colors, 
         annotation_colors = ann_colors[1],
         cluster_rows = T, 
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         treeheight_row = 0, treeheight_col = 0,
         annotation_col = df, 
         fontsize_row = 13,
         scale = "row")      
dev.off()

