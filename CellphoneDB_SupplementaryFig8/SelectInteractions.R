library(ggplot2)
library(tidyverse)
library(dplyr)
setwd("~/path/to/directory")

## Import the means and pvalues files for young and old mice obtained from CellphoneDB results.
means_young <- read.csv("./means_young.txt", sep = '\t')
pvals_young <- read.csv("./pvals_young.txt", sep = '\t')
means_old <- read.csv("./means_old.txt", sep = '\t')
pvals_old <- read.csv("./pvals_old.txt", sep = '\t')

## Select Interactions and Celltypes. CT selects the first 11 column names representing the metadata. This stays the same in the means and pvalues file. 
## Below we have 3 sets of interactions and celltype groups based on figures. You can change to plot for different celltypes.

I1 <- c("VCAM1_a4b1 complex","VCAM1_a4b7 complex","VCAM1_a9b1 complex","VCAM1_aDb2 complex","PLAUR_aVb3 complex", "PLAUR_a4b1 complex")
I2 <- c("CCL5_ACKR4", "CCL3_IDE",  "CCL4_SLC7A1", "PDCD1_FAM3C")
I3 <- c("NOTCH1_JAG2", "NOTCH2_JAG2", "NOTCH3_JAG2", "NOTCH4_JAG2", "NOTCH1_JAG1", "JAG1_NOTCH2", "JAG1_NOTCH3", "JAG1_NOTCH4")
C1 <- c("Luminal_AV.CD8_GZMK","Luminal_HS.CD8_GZMK","Myoepithelial.CD8_GZMK","Fibroblasts.CD8_GZMK")
C2 <- c("CD8_GZMK.Luminal_AV", "CD8_GZMK.Luminal_HS", "CD8_GZMK.Myoepithelial", "CD8_GZMK.Fibroblasts","MemoryCD4_CD44.Luminal_AV", "MemoryCD4_CD44.Luminal_HS", "MemoryCD4_CD44.Myoepithelial", "MemoryCD4_CD44.Fibroblasts")
C3 <- c("IL7R_ICOS.Luminal_AV", "IL7R_ICOS.Luminal_HS", "IL7R_ICOS.Myoepithelial", "IL7R_ICOS.Fibroblasts")
CT <- colnames(means_old[,c(1:11)])
C1 <- c(CT,C1)
C2 <- c(CT,C2)
C3 <- c(CT,C3)

MO1 <- means_old[means_old$interacting_pair %in% I1,colnames(means_old) %in% C1]
MY1 <- means_young[means_young$interacting_pair %in% I1,colnames(means_young) %in% C1]
PO1 <- pvals_old[pvals_old$interacting_pair %in% I1,colnames(pvals_old) %in% C1]
PY1 <- pvals_young[pvals_young$interacting_pair %in% I1,colnames(pvals_young) %in% C1]

MO2 <- means_old[means_old$interacting_pair %in% I2,colnames(means_old) %in% C2]
MY2 <- means_young[means_young$interacting_pair %in% I2,colnames(means_young) %in% C2]
PO2 <- pvals_old[pvals_old$interacting_pair %in% I2,colnames(pvals_old) %in% C2]
PY2 <- pvals_young[pvals_young$interacting_pair %in% I2,colnames(pvals_young) %in% C2]

MO3 <- means_old[means_old$interacting_pair %in% I3,colnames(means_old) %in% C3]
MY3 <- means_young[means_young$interacting_pair %in% I3,colnames(means_young) %in% C3]
PO3 <- pvals_old[pvals_old$interacting_pair %in% I3,colnames(pvals_old) %in% C3]
PY3 <- pvals_young[pvals_young$interacting_pair %in% I3,colnames(pvals_young) %in% C3]

## Add old and young to the column names to differentiate between the values for interactions in young and old mice. 
colnames(MO1)[12:15] <- paste("Old", colnames(MO1)[12:15], sep = '_')
colnames(MO2)[12:19] <- paste("Old", colnames(MO2)[12:19], sep = '_')
colnames(MO3)[12:15] <- paste("Old", colnames(MO3)[12:15], sep = '_')

colnames(PO1)[12:15] <- paste("Old", colnames(PO1)[12:15], sep = '_')
colnames(PO2)[12:19] <- paste("Old", colnames(PO2)[12:19], sep = '_')
colnames(PO3)[12:15] <- paste("Old", colnames(PO3)[12:15], sep = '_')

colnames(MY1)[12:15] <- paste("Young", colnames(MY1)[12:15], sep = '_')
colnames(MY2)[12:19] <- paste("Young", colnames(MY2)[12:19], sep = '_')
colnames(MY3)[12:15] <- paste("Young", colnames(MY3)[12:15], sep = '_')

colnames(PY1)[12:15] <- paste("Young", colnames(PY1)[12:15], sep = '_')
colnames(PY2)[12:19] <- paste("Young", colnames(PY2)[12:19], sep = '_')
colnames(PY3)[12:15] <- paste("Young", colnames(PY3)[12:15], sep = '_')


## Merges the two old and young files together so they can be plotted in one dotplot. 
M1 <- merge(MO1, MY1)
P1 <- merge(PO1, PY1)
M2 <- merge(MO2, MY2)
P2 <- merge(PO2, PY2)
M3 <- merge(MO3, MY3)
P3 <- merge(PO3, PY3)

## Save text files for cellphonedb plot dot_plot function.
write.table(M1, "M1.txt", sep = '\t',row.names = F)
write.table(P1, "P1.txt", sep = '\t',row.names = F)

write.table(M2, "M2.txt", sep = '\t',row.names = F)
write.table(P2, "P2.txt", sep = '\t',row.names = F)

write.table(M3, "M3.txt", sep = '\t',row.names = F)
write.table(P3, "P3.txt", sep = '\t',row.names = F)
