# plot counts for specific gene
e <- plotCounts(BRCA_dds_postDDS, gene= "ENSG00000164530", intgroup= "New_Column", returnData=TRUE)
e$log2 <- log2(e$count+1)
e <- subset(e, New_Column=="Normal Tissue" | New_Column=="LumA" | New_Column=="LumB")
e <- e %>% mutate(New_Column=factor(New_Column,levels=c("Normal Tissue", "LumA", "LumB")))

p <- ggplot(e, aes(x= New_Column, y=log2)) + 
  geom_violin(trim=FALSE)+
  labs(title= "PI16_ENSG00000164530") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0), text = element_text(size=20))

pdf("path/to/pdf", height=10, width=10)
p+ geom_jitter(shape=16, position=position_jitter(0.2)) + stat_compare_means(method = "t.test", 
                                                                             comparisons = list(c("Normal Tissue", "LumA"),
                                                                                                c("Normal Tissue", "LumB")))
dev.off()
