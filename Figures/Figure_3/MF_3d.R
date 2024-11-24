#####Luminal Markers
library(Seurat)
library(dplyr)
library(ggplot2)

scRNA_seq_expt_comb_harmony_umap_2 <- readRDS("path/to/object")
table(Idents(scRNA_seq_expt_comb_harmony_umap_2))
Meta <- scRNA_seq_expt_comb_harmony_umap_2@meta.data
Meta %>% dplyr::filter(Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C1"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C2"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C3"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C4"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C5"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C6"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C7") -> Luminal_Cells
dim(Luminal_Cells)
cell_to_keep <- row.names(Luminal_Cells)
length(cell_to_keep)

subset_Luminal = subset(scRNA_seq_expt_comb_harmony_umap_2, cells = cell_to_keep)
table(Idents(subset_Luminal))

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/Dotplot_Luminal.pdf", height=20, width=20)
DotPlot(subset_Luminal, features = c("Epcam", "Prlr", "Cited1", "Pgr", "Prom1", "Esr1", "Ptn", "Gpx3", "Wfdc2", "Ly6a",
                                     "Rcan1", "Fxyd2", "Tph1", "Fam3c", "Capsl",
                                     "Kit", "Aldh1a3", "Cd14", "Elf5", 
                                     "Mfge8", "Trf", "Csn3", "Wfdc18", "Ltf", "Cst3", "Igfbp5",
                                     "Hey1", "Txnip", "Crip1", "Alox15", "Alox12e", "Palmd", "Rspo1",
                                     "Mcm2", "Mcm3", "Mcm5", "Cdk1", "Mki67", "Top2a", "E2f1", "Hmgb2", "Tubb5", "Stmn1"),
        dot.scale = 5) + coord_flip() + theme(text = element_text(size=12))
dev.off()

####Luminal HS Proportions
Meta %>% dplyr::filter(Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C1"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C2"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C3") -> Luminal_HS_Cells
dim(Luminal_HS_Cells)
cell_to_keep_HS <- row.names(Luminal_HS_Cells)
length(cell_to_keep_HS)

subset_Luminal_HS = subset(scRNA_seq_expt_comb_harmony_umap_2, cells = cell_to_keep_HS)
table(Idents(subset_Luminal_HS))

###Create a proportions table 
y=prop.table(x=table(subset_Luminal_HS$orig.ident, as.character(subset_Luminal_HS$new_cluster_ids)), margin=2)

###Adding proportions for all of the cells
nbr_row <- length(unique(subset_Luminal_HS$orig.ident))
ALL <- vector("numeric", nbr_row) # prepare a container

for (i in 1:nbr_row) {
  res <- table(subset_Luminal_HS$orig.ident)[[i]]/dim(subset_Luminal_HS)[2]
  ALL[i] <- res         # change to assignment
}
ALL

y1 <- cbind(y, All=ALL)
y1

#Change to long version of the table and rename the columns
df1 <- melt(y1)
head(df1)
length(unique(subset_Luminal_HS$new_cluster_ids))
colnames(df1) <- c("orig.ident", "Clusters", "Percentage")
head(df1)

#To order the clusters
ordered_clusters <- c('All', 'C1', 'C2', 'C3')

#To color the ages
col_orig.ident = c('mm10_3mths'= "#9e9ac8",
                   'mm10_18mths'= "#66c2a4")

##Plot with ggplot
p_Age <- ggplot(data=df1, aes(x=as.character(Clusters), y=Percentage, 
                              fill=factor(orig.ident, levels=c('mm10_3mths',"mm10_18mths"))))+ 
  geom_bar(stat="identity", color="black", position =position_fill(reverse = TRUE)) +
  scale_fill_manual(values=col_orig.ident) + 
  scale_x_discrete(limits=ordered_clusters) +
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y =  element_text(size = 20))+
  theme(axis.title = element_text(size = 20))+
  labs(fill='orig.ident') + 
  ylab("Proportion of cells") +
  xlab("Clusters") +
  ggtitle("orig.ident")+
  geom_hline(yintercept = 0.6525856,linetype="dotted", size=1) 

print(p_Age)

ggsave("C:/Users/angarb/Box/Aging mammary gland paper/Figures/Figure S4 - Epi Subclusters/LuminalHS_Prop.pdf", 
       p_Age , width=3, height=1.7,  units="in", scale=3)

####Luminal AV Proportions
Meta %>% dplyr::filter(Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C4"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C5"| 
                         Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C6"| Idents(scRNA_seq_expt_comb_harmony_umap_2) == "C7") -> Luminal_AV_Cells
dim(Luminal_AV_Cells)
cell_to_keep_AV <- row.names(Luminal_AV_Cells)
length(cell_to_keep_AV)

subset_Luminal_AV = subset(scRNA_seq_expt_comb_harmony_umap_2, cells = cell_to_keep_AV)
table(Idents(subset_Luminal_AV))

###Create a proportions table
y=prop.table(x=table(subset_Luminal_AV$orig.ident, as.character(subset_Luminal_AV$new_cluster_ids)), margin=2)

###Adding proportions for all of the cells
nbr_row <- length(unique(subset_Luminal_AV$orig.ident))
ALL <- vector("numeric", nbr_row) # prepare a container

for (i in 1:nbr_row) {
  res <- table(subset_Luminal_AV$orig.ident)[[i]]/dim(subset_Luminal_AV)[2]
  ALL[i] <- res         # change to assignment
}
ALL

y1 <- cbind(y, All=ALL)
y1

#Change to long version of the table and rename the columns
df1 <- melt(y1)
head(df1)
length(unique(subset_Luminal_AV$new_cluster_ids))
colnames(df1) <- c("orig.ident", "Clusters", "Percentage")
head(df1)

#To order the clusters
ordered_clusters <- c('All', 'C4', 'C5', 'C6', 'C7')

#To color the ages
col_orig.ident = c('mm10_3mths'= "#9e9ac8",
                   'mm10_18mths'= "#66c2a4")

##Plot with ggplot
p_Age <- ggplot(data=df1, aes(x=as.character(Clusters), y=Percentage, 
                              fill=factor(orig.ident, levels=c('mm10_3mths',"mm10_18mths"))))+ 
  geom_bar(stat="identity", color="black", position =position_fill(reverse = TRUE)) +
  scale_fill_manual(values=col_orig.ident) + 
  scale_x_discrete(limits=ordered_clusters) +
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y =  element_text(size = 20))+
  theme(axis.title = element_text(size = 20))+
  labs(fill='orig.ident') + 
  ylab("Proportion of cells") +
  xlab("Clusters") +
  ggtitle("orig.ident")+
  geom_hline(yintercept = 0.5493166,linetype="dotted", size=1) 

print(p_Age)

ggsave("path/to/pdf", 
       p_Age , width=3, height=1.7,  units="in", scale=3)
