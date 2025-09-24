## ---- 0. Install & Load Packages -----------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "enrichplot",   
  "DOSE",         
  "org.Hs.eg.db", 
  "GO.db",       
  "AnnotationDbi"
), update = TRUE, ask = FALSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GO.db)
library(AnnotationDbi)

files <- c(
  "diff_expression_t0_t24.csv",
  "diff_expression_t0_t48.csv",
  "diff_expression_t0_t72.csv",
  "diff_expression_t24_t48.csv",
  "diff_expression_t24_t72.csv",
  "diff_expression_t48_t72.csv"
)

gsea_results <- lapply(files, function(f) {
  

  res <- read.csv(f)
  
  gene.df <- bitr(res$symbol,
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = org.Hs.eg.db)
  
  res2 <- merge(res, gene.df, by.x = "symbol", by.y = "SYMBOL")
  
  res2 <- res2 %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::summarize(log2FoldChange = log2FoldChange[which.max(abs(log2FoldChange))]) %>%
    dplyr::ungroup()
  
  geneList <- res2$log2FoldChange
  names(geneList) <- res2$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  
  if (length(geneList) == 0) {
    message(f, ": no valid genes after mapping")
    return(NULL)
  }
  

  gsea <- gseGO(
    geneList     = geneList,
    OrgDb        = org.Hs.eg.db,
    ont          = "BP",
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    verbose      = FALSE,
    BPPARAM      = BiocParallel::SerialParam()
  )
  
  comp <- sub("diff_expression_(.*)\\.csv", "\\1", basename(f))
  out_file <- paste0("CP_enrichment_results_", comp, ".csv")
  write.csv(as.data.frame(gsea), file = out_file, row.names = FALSE)
  
  return(gsea)
})


names(gsea_results) <- basename(files)

ridgeplot(gsea_results[[6]])
