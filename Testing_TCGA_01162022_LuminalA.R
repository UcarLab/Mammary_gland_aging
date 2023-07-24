# 01162023
# Attempting to generate a DGE list from  Luminal A tumors in TCGA vs all normal samples
# Pulling counts from TCGA biolinks

## http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#### Sites on TCGA Biolinks
## https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/analysis.html
## https://gist.github.com/mbk0asis/
## https://www.biostars.org/p/9518549/
## https://support.bioconductor.org/p/9143118/#9143119
## https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/analysis.html
## https://gist.github.com/mbk0asis/5db572a9182075c012c1a04567bbfe42

#### Sites on correlations
## https://bookdown.org/daniel_dauber_io/r4np_book/correlations.html
## https://stackoverflow.com/questions/70656226/multiple-lines-in-ggscatter-plot

options(stringsAsFactors = FALSE)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(biomaRt)
library(vsn)
library(ggplot2)
library(RColorBrewer)
library(readxl)
library(pheatmap)
library(rtracklayer)
library(data.table)
library(enrichR)
library(ggrepel)
library(apeglm)
library(plotly)
library(tibble)
library(dplyr)
library(ggpubr)
library(correlation)

## Uploading gene name conversion list (with biotype)
geneName_ID_and_biotype <- read.csv("C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/geneName_ID_and_biotype.csv")

#### Downloading samples from TCGA Biolinks
## Using the 1186 samples that passed LU's quality filters
## The metadata from samples found here: C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/AllTumors_07262022/Inputs/BRCA_metadata.csv
listSamples <- scan("C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/AllTumors_07262022/Inputs/BRCA_TCGA.txt", what = "", quiet=TRUE)

project_name <- "TCGA-BRCA"
query <- GDCquery(
  project = project_name, 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  barcode = listSamples
)

metadata <- query[[1]][[1]]
# write.csv(metadata, "C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/AllTumors_07262022/Intermediate_files/GDC_metadata.csv") 

GDCdownload(query, method = "api", files.per.chunk = 5)

BRCA.Rnaseq <- GDCprepare(query)
BRCAMatrix_Initial <- as.data.frame(assay(BRCA.Rnaseq,"unstranded"))
# write.csv(BRCAMatrix_Initial, "C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/AllTumors_07262022/Intermediate_files/GDC_counts_matrix.csv") 

### Reformat count matrix to get the desired format
# Remove rows with "PAR_Y", these are genes on both X and Y 
# (https://www.gencodegenes.org/pages/faq.html)
counts_initial<- subset(BRCAMatrix_Initial, grepl("PAR_Y", rownames(BRCAMatrix_Initial))==FALSE)

# Reformat gencode rownames
reformat_ensembl_ID <- function(val){
  new_val <- strsplit(val, "[.]")[[1]][1]
  return(new_val)
}
counts_initial$new_id <- sapply(rownames(counts_initial), function(x) reformat_ensembl_ID(x))
rownames(counts_initial) <- counts_initial$new_id
BRCAMatrix <- counts_initial[,c(-1187)]

#### Load in LU data from BRCA TCGA
dds_colData <- read.csv("C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/BRCA_metadata.csv")

# Remove metastatic samples (ShortLetterCode = 'TM') and males (gender = male and Blank))
brca_dds_colData_final <- subset(dds_colData, dds_colData$ShortLetterCode != 'TM' & dds_colData$gender == "female")

# Make new column with better "Normal" annotations 
# Labeling NT as "Normal Tissue" and Normal TP as "Normal-like"
# Overwriting the "LumA" NT as "Normal Tissue"
brca_dds_colData_final<-brca_dds_colData_final%>%mutate(New_Column = case_when(
  Subtype_mRNA=="Basal" & ShortLetterCode == "TP" ~ "Basal",
  Subtype_mRNA=="Her2" & ShortLetterCode == "TP" ~ "Her2",
  Subtype_mRNA=="LumA" & ShortLetterCode == "TP" ~ "LumA",
  Subtype_mRNA=="LumB" & ShortLetterCode == "TP" ~ "LumB",
  Subtype_mRNA=="Normal" & ShortLetterCode == "TP" ~ "Normal-like",
  ShortLetterCode == "NT" ~ "Normal Tissue"
))
table(brca_dds_colData_final$New_Column)
table(brca_dds_colData_final$Subtype_mRNA)
table(brca_dds_colData_final$ShortLetterCode)

# Subset Normal and Luminal A
brca_dds_colData_final <- subset(brca_dds_colData_final, brca_dds_colData_final$New_Column == 'LumA' | brca_dds_colData_final$ShortLetterCode == "NT")

# Subsetting the data
tumor_quant <- subset(brca_dds_colData_final, brca_dds_colData_final$ShortLetterCode == 'TP')
table(tumor_quant$New_Column)

normal_quant <- subset(brca_dds_colData_final, brca_dds_colData_final$ShortLetterCode == 'NT')
table(normal_quant$New_Column)

# Looking at age histogram
age_1 <- ggplot(tumor_quant, aes(x=age_at_index,fill = Subtype_mRNA)) + geom_histogram(colour= 'black', bins = 5)
print(age_1)
age_2 <- ggplot(tumor_quant, aes(x=age_at_index)) + geom_histogram(colour= 'black', bins = 5)
print(age_2)
age_3 <- density(tumor_quant$age_at_index) # returns the density data
plot(age_3)

hist(tumor_quant$age_at_index)

#Subset counts data based on metadata and validating in correct order
brca_all_counts <- BRCAMatrix[,colnames(BRCAMatrix) %in% brca_dds_colData_final$associated_entities.entity_submitter_id]
dim(brca_all_counts)

dds_colData_all <- brca_dds_colData_final[order(match(brca_dds_colData_final$associated_entities.entity_submitter_id,colnames(brca_all_counts))),]
row.names(dds_colData_all) = dds_colData_all$associated_entities.entity_submitter_id
all(rownames(dds_colData_all) %in% colnames(brca_all_counts))
all(rownames(dds_colData_all) == colnames(brca_all_counts))
    
# create dds object
model.matrix(~ShortLetterCode, dds_colData_all)
dds <- DESeqDataSetFromMatrix(countData=brca_all_counts, 
                              colData=dds_colData_all, 
                              design = ~ ShortLetterCode)
dds$ShortLetterCode <- factor(dds$ShortLetterCode, levels = c("NT","TP"))

# pre-filtering
# note that this is really only to help speed, filtering is part of DESeq2 downstream
row_sums <- as.data.frame(rowSums(as.data.frame(counts(dds))))
keep <- row_sums>=10
dds <- dds[keep,]


## PCA plots ---- 
brca_vsd_blind <- vst(dds, blind=TRUE)
brca_vsd_NOTblind <- vst(dds, blind=FALSE)
meanSdPlot(assay(brca_vsd_blind))
meanSdPlot(assay(brca_vsd_NOTblind))

pca.res <- prcomp(t(assay(brca_vsd_NOTblind)))
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)

#classification plot with specific colors
plot_LU_pca <- function(pca_df, var_name){
  ggplot(pca_df) +
    aes(x=PC1, y=PC2, 
        label= dds_colData_all$sample_name, 
        color = dds_colData_all[,var_name]) +
    geom_point(size=4) +
    xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
    ylab(paste0("PC2 (",pc.per[2],"%",")")) +
    labs(title="PCA plot",
         caption=paste0("produced on ", Sys.time())) +
    coord_fixed() +
    theme_bw()+
    theme(text=element_text(size=20))
}
plot_LU_pca(pca.res.df, var_name="ShortLetterCode")
plot_LU_pca(pca.res.df, var_name="New_Column")
plot_LU_pca(pca.res.df, var_name="Subtype_mRNA")
plot_LU_pca(pca.res.df, var_name="age_at_index")
plot_LU_pca(pca.res.df, var_name="tumor_stage")

# run DESeq
BRCA_dds_postDDS <- DESeq(dds)
resultsNames(BRCA_dds_postDDS)

#### getting results 
get_results <- function(df, comparison){
  new_df <- results(df, name = comparison, alpha = 0.05)
  print(summary(new_df))
  return(new_df)
}
res_all <- get_results(BRCA_dds_postDDS, "ShortLetterCode_TP_vs_NT")
mcols(res_all)
sum(res_all$padj < 0.01, na.rm=TRUE)

# convert res to df and adding ensembl ids to new column
res_df <- as.data.frame(res_all)
res_df$ensembl_gene_id <- rownames(res_df)

# get res dataframe and merge with gene names 
res_df <- merge(res_df, geneName_ID_and_biotype, by.x="ensembl_gene_id", by.y="ensembl_gene_id", all.x=TRUE )

res_df$significant <- ifelse(res_df$padj<=0.05 & abs(res_df$log2FoldChange) > 0.5, "Yes", "No")
res_df$direction <- ifelse(res_df$log2FoldChange > 0, "up", "down")
res_df$ensembl_gene_id <- as.character(res_df$ensembl_gene_id)
table(res_df[,c("significant", "direction")])
#    down    up
#No  18507 14910
#Yes 10704  9495
res_significant <- subset(res_df, significant =="Yes")
sum(res_significant$log2FoldChange<0)
#9467
sum(res_significant$log2FoldChange>0)
#9495

write.csv(res_df, "C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/AllTumors_07262022/Outputs/TCGA_LuminalAvsNormal_DE.csv") 

# plot counts for specific gene
e <- plotCounts(BRCA_dds_postDDS, gene= "ENSG00000073111", intgroup= "New_Column", returnData=TRUE)
e$log2 <- log2(e$count+1)
e <- e %>% mutate(New_Column=factor(New_Column,levels=c("Normal Tissue", "LumA")))

p <- ggplot(e, aes(x= New_Column, y=log2)) + 
  geom_violin(trim=FALSE)+
  labs(title= "MCM2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0), text = element_text(size=20))

p+ geom_jitter(shape=16, position=position_jitter(0.2))

f <- plotCounts(BRCA_dds_postDDS, gene= "ENSG00000073111", intgroup= "ShortLetterCode", returnData=TRUE)
f$log2 <- log2(f$count+1)

p <- ggplot(f, aes(x= ShortLetterCode, y=log2)) + 
  geom_violin(trim=FALSE)+
  labs(title= "MCM2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0), text = element_text(size=20))

p+ geom_jitter(shape=16, position=position_jitter(0.2)) + stat_compare_means(method = "t.test")

# extract normalized counts 
BRCA_norm_counts <- counts(BRCA_dds_postDDS, normalized=TRUE)
write.csv(BRCA_norm_counts, "C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/AllTumors_07262022/Outputs/Normcounts_TCGA_LumAvsNormal.csv") 

BRCA_norm_counts.log <- log(BRCA_norm_counts+1)
write.csv(BRCA_norm_counts.log2, "C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/AllTumors_07262022/Outputs/Log_Normcounts_TCGA_TumorvsNormal.csv") 
