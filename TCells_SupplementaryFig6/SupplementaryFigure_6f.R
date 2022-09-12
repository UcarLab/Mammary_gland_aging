library(VennDiagram)
setwd("~path/to/directory")
GZMK <- read.csv("~/path/to/CD8_GZMK_6f.txt", sep = '\t')
GZMM <- read.csv("~/path/to/CD8_GZMM_6f.txt", sep = '\t')
CD44 <- read.csv("~/path/to/Memory_CD4_6f.txt", sep = '\t')

GZMK <- as.data.frame(GZMK$gene)
GZMM <- as.data.frame(GZMM$gene)
CD44 <- as.data.frame(CD44$gene)

colnames(GZMK) <- "gene"
colnames(GZMM) <- "gene"
colnames(CD44) <- "gene"

myCol <- c("purple", "dodgerblue3","firebrick4")
p1 <- venn.diagram(
  x = list(GZMK$gene, GZMM$gene, CD44$gene),
  category.names = c("GZMK" , "GZMM" , "CD44"),
  filename = NULL,
  output=TRUE,
  
  # Output features
  imagetype="pdf" ,
  height = 300 , 
  width = 300 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


pdf("~/path/to/SFig_6f.pdf", height = 4, width = 4)
grid.draw(p1)
dev.off()
