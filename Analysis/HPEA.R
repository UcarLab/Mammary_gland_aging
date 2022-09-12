HPEA <- function(modules, genes, 
                 union.size = 20e3, geneset.name = "", contrast = "") {
  
  # Create empty p value, module name vectors
  module_names        <- unique(modules[,"Module.Name"])
  module_count        <- module_names %>% length
  vector_p            <- numeric(0) 
  overlapping_genes   <- character(0)
  
  for (module_name in module_names){
    
    # Get the module genes
    module      <- modules[modules[,"Module.Name"] == module_name,"GeneName"]
    
    # Module gene count
    k <- unique(module) %>% length
    
    # gene count
    n <- length(genes)
    
    # overlapped up-regulated genes with the module
    q <- genes %in% module %>% sum 
    
    
    # here we calculate the probability of having a bigger intersection
    # than the count of overlapping genes given the module size and the total gene count.
    # we substract 1 for removing the equality when the lower.tail = F, which changes P(X<x) to 1-P(X>=x).
    vector_p[module_name] <- phyper(q-1, k, union.size - k, n, lower.tail = F, log.p = F)
    
    # take the overlapping genes for the modules
    og <- genes[genes %in% module]
    overlapping_genes [module_name] <- ifelse(length(og) <= 0, "", paste(og, collapse = ","))
  }
  
  df_modules   <- data.frame(geneset.name = geneset.name,
                             contrast = contrast,
                             module.name = module_names,
                             overlapping.genes = overlapping_genes,
                             p = vector_p,
                             stringsAsFactors = F)
  
  # correct p-values for multiple hypothesis
  df_modules$adj.p   <- p.adjust(p = df_modules$p  , method = "fdr")
  
  # sort according to adjusted p-values and then to p-values
  df_modules   <- df_modules  [order(df_modules$adj.p, df_modules$p),]
  
  return(df_modules)
}
