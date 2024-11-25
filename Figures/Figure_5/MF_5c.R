library(Seurat)
TCells <- readRDS("/path/to/object")
Idents(TCells) <- "predicted.id"
colr <- c("springgreen4","goldenrod2","dodgerblue3","purple","lightblue2", "deepskyblue", "firebrick4","orangered1", "lightpink")
p1 <- DimPlot(TCells, reduction = "umap3", group.by = "predicted.id", cols = rev(colr),label = F, label.size = 8)
pdf("./TCells_ATAC_Fig1A.pdf", height = 6, width = 8)
p1
dev.off()


##For the barplot
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

pdf("path/to/pdf", height = 6, width = 8)
p1
dev.off()

## For the significance values
stat1 <- compare_means(Ratio ~ predicted.id,  data = NewMeta, ref.group = "Base", method = "t.test")
write.table(stat1, "path/to/table", row.names = F, sep = '\t')
