library(VennDiagram)
setwd("~path/to/directory")
GZMK <- read.csv("~/path/to/CD8_GZMK_6i.txt", sep = '\t')
GZMM <- read.csv("~/path/to/CD8_GZMM_6i.txt", sep = '\t')
CD44 <- read.csv("~/path/to/Memory_CD4_6i.txt", sep = '\t')

GZMK <- as.data.frame(paste(GZMK$seqnames,GZMK$start,GZMK$end,sep = '_'))
GZMM <- as.data.frame(paste(GZMM$seqnames,GZMM$start,GZMM$end,sep = '_'))
CD44 <- as.data.frame(paste(CD44$seqnames,CD44$start,CD44$end,sep = '_'))

colnames(GZMK) <- "peak"
colnames(GZMM) <- "peak"
colnames(CD44) <- "peak"

myCol <- c("purple", "dodgerblue3","firebrick4")
p1 <- venn.diagram(
  x = list(GZMK$peak, GZMM$peak, CD44$peak),
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


pdf("~/path/to/SFig_6i.pdf", height = 4, width = 4)
grid.draw(p1)
dev.off()
