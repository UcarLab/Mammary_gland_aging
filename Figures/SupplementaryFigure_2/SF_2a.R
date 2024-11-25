library(VennDiagram)
## Needs list of differential genes
LumAV <- read.csv("/path/to/Luminal-AV diff genes", sep = '\t')
LumHS <- read.csv("/path/to/Luminal-HS diff genes", sep = '\t')
Myo <- read.csv("/path/to/Myoepithelial diff genes", sep = '\t')

LumAV <- as.data.frame(paste(LumAV$seqnames,LumAV$start,LumAV$end,sep = '_'))
LumHS <- as.data.frame(paste(LumHS$seqnames,LumHS$start,LumHS$end,sep = '_'))
Myo <- as.data.frame(paste(Myo$seqnames,Myo$start,Myo$end,sep = '_'))

colnames(LumAV) <- "gene"
colnames(LumHS) <- "gene"
colnames(Myo) <- "gene"

myCol <- c("purple", "dodgerblue3","firebrick4")
p1 <- venn.diagram(
  x = list(LumAV$gene, LumHS$gene, Myo$gene),
  category.names = c("LumAV" , "LumHS" , "Myo"),
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


#grid.draw(p1)
# 
pdf(file="path/to/pdf", height = 4, width = 4)
grid.draw(p1)
dev.off()
