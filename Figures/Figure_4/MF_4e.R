library(Seurat)
library(dplyr)
library(ggplot2)
scRNA_seq_expt_comb_harmony_umap_3 <- readRDS("path/to/seuratobject")
###Dotplot
DotPlot(subset_Fibro, features = c("Penk", "Fbln1","Tnc", "Fmo1", "Postn", "Mfap4", "Spon1","Crabp1","Col15a1", 
                                   "Col6a3", "Col5a3", "Col4a1", "Col4a2", "Lpl", "Fabp4", "Pparg", "Fap",
                                   "Gdf10", "F3", "Hmcn2", "Gpc3", "Dpt", "Ly6a",
                                   "Pi16", "Dpp4", "Ly6c1", "Tek", "Aldh1a3","Fn1", "Col14a1","Sema3c", "Esr1",
                                   "Ccl19", "Cxcl12", "Cd74"),
        dot.scale = 5) + coord_flip() + theme(text = element_text(size=12))

#################################################################################################
###Create a proportions table
y=prop.table(x=table(subset_Fibro$orig.ident, as.character(subset_Fibro$seurat_clusters)), margin=2)

###Adding proportions for all of the cells
nbr_row <- length(unique(subset_Fibro$orig.ident))
ALL <- vector("numeric", nbr_row) # prepare a container

for (i in 1:nbr_row) {
  res <- table(subset_Fibro$orig.ident)[[i]]/dim(subset_Fibro)[2]
  ALL[i] <- res         # change to assignment
}
ALL

y1 <- cbind(y, All=ALL)
y1

#Change to long version of the table and rename the columns
df1 <- melt(y1)
head(df1)
length(unique(subset_Fibro$seurat_clusters))
colnames(df1) <- c("orig.ident", "Clusters", "Percentage")
head(df1)

#To order the clusters
ordered_clusters <- c('All', '0', '6', '1', '4', '3','10','7')

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
  geom_hline(yintercept = 0.2683663,linetype="dotted", size=1) 

print(p_Age)

ggsave("path/to/pdf", 
       p_Age , width=3, height=1.7,  units="in", scale=3)

