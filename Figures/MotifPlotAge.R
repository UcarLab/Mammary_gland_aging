## FIGURE 2D - Motifs by Age
Motifs <- read.csv("path/to/motifs", sep = '\t', header = F)
Motifs <- Motifs[-1,]
Motifs <- Motifs[,-1]
colnames(Motifs) <- c("ID","MotifName", "p_val", "log_FC_Change", "pct.1", "pct.2", "adj_p_val", "CellType")
Motifs$CellType <- gsub(".txt","",Motifs$CellType)
Motifs$log_FC_Change <- as.numeric(Motifs$log_FC_Change)
Motifs <- Motifs[- grep("::", Motifs$MotifName),]
Motifs2 <- Motifs
Motifs <- Motifs[(Motifs$log_FC_Change >= 0.5 | Motifs$log_FC_Change <= -0.5),]
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
write.table(FinalMotifs,"path/to/output", sep = '\t', row.names = F)
Clustering <- FinalMotifs[,c(1,2,4)]
Clustering1 <- cast(Clustering, MotifName ~ CellType, mean, value = 'log_FC_Change')
rownames(Clustering1) <- Clustering1$MotifName
Clustering1 <- Clustering1[,-1]
data <- scale(t(Clustering1))

ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
ord

#FinalMotifs$CellType <- factor(FinalMotifs$CellType, levels = c("NaiveCD4_LEF1","CD4_CCR7_JUN","IL7R_ICOS","NK","CD8_CCR7","CD8_CCR9","MemoryCD4_CD44", "CD8_GZMK", "CD8_GZMM"))
FinalMotifs$CellType <- factor(FinalMotifs$CellType, levels = c("Luminal-AV","Luminal-HS","Myoepithelial"))

RS <- FinalMotifs[,c(1,5)]
RS1 <- RS %>% group_by_all() %>% summarise(COUNT = n())
RS1 <- RS1[!(RS1$Significance=="" & RS1$COUNT==3),]
RS1 <- as.data.frame(unique(RS1$MotifName))
colnames(RS1) <- "MotifName"

FinalMotifs <- merge(RS1, FinalMotifs)

pdf("path/to/pdf", height = 12, width = 12)
ggplot(FinalMotifs, aes(MotifName, CellType, fill= log_FC_Change)) +
  geom_tile() +
  scale_fill_distiller(palette = "PuOr", ) +
  geom_text(aes(label=Significance), color="black", size=4) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
