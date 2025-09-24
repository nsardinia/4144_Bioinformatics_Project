## ---- 0. Install & Load Packages -----------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "enrichplot",   # for ridgeplot, dotplot, gseaplot2
  "DOSE",         # optional, for Disease Ontology
  "org.Hs.eg.db", # human gene annotations
  "GO.db",        # GO term definitions
  "AnnotationDbi" # for mapping IDs
), update = TRUE, ask = FALSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GO.db)
library(AnnotationDbi)

## ---- 1. Files to Process -------------------------------------------------
files <- c(
  "diff_expression_t0_t24.csv",
  "diff_expression_t0_t48.csv",
  "diff_expression_t0_t72.csv",
  "diff_expression_t24_t48.csv",
  "diff_expression_t24_t72.csv",
  "diff_expression_t48_t72.csv"
)

## ---- 2. GSEA Loop --------------------------------------------------------
gsea_results <- lapply(files, function(f) {
  
  # Read differential expression results
  res <- read.csv(f)
  
  # Convert gene symbols to Entrez IDs
  gene.df <- bitr(res$symbol,
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = org.Hs.eg.db)
  
  # Merge to keep only rows with valid Entrez IDs
  res2 <- merge(res, gene.df, by.x = "symbol", by.y = "SYMBOL")
  
  # Collapse duplicates: keep gene with max absolute log2FC
  res2 <- res2 %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::summarize(log2FoldChange = log2FoldChange[which.max(abs(log2FoldChange))]) %>%
    dplyr::ungroup()
  
  # Create ranked gene list
  geneList <- res2$log2FoldChange
  names(geneList) <- res2$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Skip if geneList is empty
  if (length(geneList) == 0) {
    message(f, ": no valid genes after mapping")
    return(NULL)
  }
  
  # Run GSEA for GO Biological Process
  gsea <- gseGO(
    geneList     = geneList,
    OrgDb        = org.Hs.eg.db,
    ont          = "BP",
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    verbose      = FALSE,
    BPPARAM      = BiocParallel::SerialParam() # disable parallelization
  )
  
  # Save results to CSV
  comp <- sub("diff_expression_(.*)\\.csv", "\\1", basename(f))
  out_file <- paste0("CP_enrichment_results_", comp, ".csv")
  write.csv(as.data.frame(gsea), file = out_file, row.names = FALSE)
  
  return(gsea)
})

# Name the list elements with the file names
names(gsea_results) <- basename(files)


## ---- 3. Example Plot -----------------------------------------------------
# View a ridgeplot for the first comparison
ridgeplot(gsea_results[[1]])
