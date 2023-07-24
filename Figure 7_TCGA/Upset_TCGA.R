#01242023
#Making UpSet plots
#https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html

library(UpSetR)
library(ComplexHeatmap)
library(dplyr)
library(purrr)
library(RVenn)
library(ggplot2)
library(readxl)

#Loading Tumor Data
All_Tumors <- read.csv("C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/AllTumors_07262022/Outputs/TCGA_TumorvsNormal_DE.csv")
LuminalA <- read.csv("C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/AllTumors_07262022/Outputs/TCGA_LuminalAvsNormal_DE.csv")
LuminalB <- read.csv("C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/AllTumors_07262022/Outputs/TCGA_LuminalBvsNormal_DE.csv")
Basal <- read.csv("C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/AllTumors_07262022/Outputs/TCGA_BasalvsNormal_DE.csv")
Her2 <- read.csv("C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/AllTumors_07262022/Outputs/TCGA_Her2vsNormal_DE.csv")
Normallike <- read.csv("C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/TCGA_MatchDGE/AllTumors_07262022/Outputs/TCGA_Normal-likevsNormal_DE.csv")

mousetohuman <- read.csv("C:/Users/angarb/OneDrive - The Jackson Laboratory/Documents/Anczukow Lab/Data/human_and_mouse_conversion.csv")

#Filter for significant FDR and separate into positive logFC
upregulated_genes <- function(df){
  new_df <- df
  sig_df <- subset(new_df, new_df$significant == 'Yes' & new_df$direction == 'up')
  print(sprintf("There were %i significant upregulated genes", nrow(sig_df)))
  mouse_sig_df <- merge(sig_df, mousetohuman[c("mouse_gene_name","human_ensembl")],by.x=c("ensembl_gene_id"),by.y=c("human_ensembl"))
  print(sprintf("There were %i significant upregulated mouse genes", nrow(mouse_sig_df)))
  gene_list <- mouse_sig_df$mouse_gene_name
}

All_Tumors_SigUp_Genes <- upregulated_genes(All_Tumors)
LuminalA_SigUp_Genes <- upregulated_genes(LuminalA)
LuminalB_SigUp_Genes <- upregulated_genes(LuminalB)
Basal_SigUp_Genes <- upregulated_genes(Basal)
Her2_SigUp_Genes <- upregulated_genes(Her2)
Normallike_SigUp_Genes <- upregulated_genes(Normallike)

downregulated_genes <- function(df){
  new_df <- df
  sig_df <- subset(new_df, new_df$significant == 'Yes' & new_df$direction == 'down')
  print(sprintf("There were %i significant downregulated genes", nrow(sig_df)))
  mouse_sig_df <- merge(sig_df, mousetohuman[c("mouse_gene_name","human_ensembl")],by.x=c("ensembl_gene_id"),by.y=c("human_ensembl"))
  print(sprintf("There were %i significant downregulated mouse genes", nrow(mouse_sig_df)))
  gene_list <- mouse_sig_df$mouse_gene_name
}

All_Tumors_SigDown_Genes <- downregulated_genes(All_Tumors)
LuminalA_SigDown_Genes <- downregulated_genes(LuminalA)
LuminalB_SigDown_Genes <- downregulated_genes(LuminalB)
Basal_SigDown_Genes <- downregulated_genes(Basal)
Her2_SigDown_Genes <- downregulated_genes(Her2)
Normallike_SigDown_Genes <- downregulated_genes(Normallike)

DE <- read_excel("C:/Users/angarb/Box/Aging mammary gland paper/Manuscript/Tables/Supplementary Table 3 - DE and DA epithelial.xlsx", sheet = 2, skip = 3)
Sig_DE <- subset(DE, DE$`Significant in`== 'Both'|DE$`Significant in`== 'Pseudobulk'|DE$`Significant in`== 'Single_cell')

Upset_List_Up <- function(cluster,direction){
  new_df <- subset(Sig_DE, Sig_DE$Cluster == cluster & Sig_DE$Direction == direction)
  new_list <- new_df$`Gene Name`
  list_of_lists = list(Luminal_A = LuminalA_SigUp_Genes,
                  Luminal_B = LuminalB_SigUp_Genes,
                  Basal = Basal_SigUp_Genes,
                  Her2 = Her2_SigUp_Genes,
                  Normal_Like = Normallike_SigUp_Genes,
                  Myoepithelial = new_list)}

Upset_List_Down <- function(cluster,direction){
  new_df <- subset(Sig_DE, Sig_DE$Cluster == cluster & Sig_DE$Direction == direction)
  new_list <- new_df$`Gene Name`
  list_of_lists = list(Luminal_A = LuminalA_SigDown_Genes,
                       Luminal_B = LuminalB_SigDown_Genes,
                       Basal = Basal_SigDown_Genes,
                       Her2 = Her2_SigDown_Genes,
                       Normal_Like = Normallike_SigDown_Genes,
                       Myoepithelial = new_list)}

lt_up_hs <- Upset_List_Up('Luminal-HS','upregulated')
lt_down_hs <- Upset_List_Down('Luminal-HS','downregulated')

lt_up_av <- Upset_List_Up('Luminal-AV','upregulated')
lt_down_av <- Upset_List_Down('Luminal-AV','downregulated')

lt_up_myo <- Upset_List_Up('Myoepithelial','upregulated')
lt_down_myo <- Upset_List_Down('Myoepithelial','downregulated')

#upset(fromList(lt_up_myo), sets = c("Myoepithelial", "Luminal_A", "Luminal_B", "Basal", "Her2", "Normal_Like"), 
#    nintersects = NA, group.by = "sets",  keep.order = TRUE, order.by = "freq", empty.intersections = "on", text.scale=0.3)

Genes_upset <- function(df,param){upset(fromList(df), sets = c(param, "Luminal_A", "Luminal_B", "Basal", "Her2", "Normal_Like"), 
                                                    intersections = list(list(param,"Luminal_A"), 
                                                                         list(param,"Luminal_B"),
                                                                         list(param,"Basal"),
                                                                         list(param,"Her2"),
                                                                         list(param,"Normal_Like"),
                                                                         list(param, "Luminal_A", "Luminal_B", "Basal", "Her2", "Normal_Like"),
                                                                         list(param, "Luminal_A", "Luminal_B", "Basal", "Her2"),
                                                                         list(param, "Luminal_A", "Luminal_B", "Basal", "Normal_Like"),
                                                                         list(param, "Luminal_A", "Luminal_B", "Her2", "Normal_Like"),
                                                                         list(param, "Luminal_A", "Basal", "Her2", "Normal_Like"),
                                                                         list(param, "Luminal_B", "Basal", "Her2", "Normal_Like"),
                                                                         list(param, "Luminal_A", "Luminal_B", "Basal"),
                                                                         list(param, "Luminal_A", "Luminal_B", "Her2"),
                                                                         list(param, "Luminal_A", "Luminal_B", "Normal_Like"),
                                                                         list(param, "Luminal_A", "Basal", "Normal_Like"),
                                                                         list(param, "Luminal_A", "Basal", "Her2"),
                                                                         list(param, "Luminal_A", "Normal_Like", "Her2"),
                                                                         list(param, "Luminal_B", "Normal_Like", "Her2"),
                                                                         list(param, "Luminal_B", "Basal", "Her2"),
                                                                         list(param, "Luminal_B", "Basal", "Normal_Like"),
                                                                         list(param, "Basal", "Her2", "Normal_Like"),
                                                                         list(param, "Luminal_A", "Luminal_B"),
                                                                         list(param, "Luminal_A", "Basal"),
                                                                         list(param, "Luminal_A", "Her2"),
                                                                         list(param, "Luminal_A", "Normal_Like"),
                                                                         list(param, "Luminal_B", "Basal"),
                                                                         list(param, "Luminal_B", "Her2"),
                                                                         list(param, "Luminal_B", "Normal_Like"),
                                                                         list(param, "Basal", "Normal_Like"),
                                                                         list(param, "Basal", "Her2"),
                                                                         list(param, "Her2", "Normal_Like"),
                                                                         list(param)), 
                                                    order.by = "freq")
  pdf(file = outfile2, width=30, height=10)
  print(ggdraw()+draw_plot(pie_only))
  dev.off()}

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Manuscript/Figures/Figure 6 Age related genes in tumors/Upset Plots/Luminal_HS_Up.pdf", height=5, width=6)
Genes_upset(lt_up_hs, "Luminal_HS")
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Manuscript/Figures/Figure 6 Age related genes in tumors/Upset Plots/Luminal_HS_Down.pdf", height=5, width=6)
Genes_upset(lt_down_hs, "Luminal_HS")
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Manuscript/Figures/Figure 6 Age related genes in tumors/Upset Plots/Luminal_AV_Up.pdf", height=5, width=6)
Genes_upset(lt_up_av, "Luminal_AV")
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Manuscript/Figures/Figure 6 Age related genes in tumors/Upset Plots/Luminal_AV_Down.pdf", height=5, width=6)
Genes_upset(lt_down_av, "Luminal_AV")
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Manuscript/Figures/Figure 6 Age related genes in tumors/Upset Plots/Myo_Up.pdf", height=5, width=6)
Genes_upset(lt_up_myo, "Myoepithelial")
dev.off()

pdf("C:/Users/angarb/Box/Aging mammary gland paper/Manuscript/Figures/Figure 6 Age related genes in tumors/Upset Plots/Myo_Down.pdf", height=5, width=6)
Genes_upset(lt_down_myo, "Myoepithelial")
dev.off()
