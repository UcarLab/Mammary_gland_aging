library(Seurat)
TCells <- readRDS("/path/to/object")
Idents(TCells) <- "NewAnnotations"
colr <- c("springgreen4","goldenrod2","darkblue","dodgerblue3","purple","lightblue2", "deepskyblue", "firebrick4","orangered1", "lightpink")
p1 <- DimPlot(TCells, reduction = "humap33", group.by = "NewAnnotations", cols = colr,label = T, label.size = 5)
pdf("path/to/pdf", height = 6, width = 8)
p1
dev.off()

## For the barplot
Meta <- TCells@meta.data
Meta <- Meta[,c(1,4,29)]
Meta <- Meta %>% group_by_all() %>% summarise(COUNT = n())
Meta$replicate <- gsub("M3","M18",Meta$replicate)
NewMeta = cast(Meta, replicate + NewAnnotations ~ orig.ident, sum) 
NewMeta$Ratio <- log2(NewMeta$mm10_18mths/NewMeta$mm10_3mths)
NewMeta <- NewMeta[,c(1,2,5)]
NewMeta$NewAnnotations <- as.character(NewMeta$NewAnnotations)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep1","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep2","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep3","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep4","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep5","Base",0)
NewMeta[nrow(NewMeta) + 1,] <- c("M18_rep6","Base",0)
NewMeta$NewAnnotations <- as.factor(NewMeta$NewAnnotations)
NewMeta$Ratio <- as.numeric(NewMeta$Ratio)

#colr1 <- c("orangered1","deepskyblue","lightblue2","darkblue","purple","dodgerblue3", "goldenrod2", "firebrick4", "lightpink", "springgreen4")
colr <- c("black","springgreen4","goldenrod2","darkblue","dodgerblue3","purple","lightblue2", "deepskyblue", "firebrick4","orangered1", "lightpink")

p1 <- ggplot(NewMeta,aes(x = reorder(NewAnnotations, -Ratio), y = Ratio, fill = NewAnnotations)) + 
  geom_boxplot()+
  scale_fill_manual(values=(colr)) +
  xlab("CellTypes") +
  ylab("log2 FC (Old/Young)") +
  ggtitle("Cell Type Proportions") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  stat_compare_means(method = "anova", label.y = 10)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = "Base", hide.ns = T, label.y = 8)   
p1                                                                            
pdf("path/to/pdf", height = 7, width = 8)
p1
dev.off()

## For the significance values
stat1 <- compare_means(Ratio ~ NewAnnotations,  data = NewMeta, ref.group = "Base", method = "t.test")
write.table(stat1, "path/to/table", row.names = F, sep = '\t')
