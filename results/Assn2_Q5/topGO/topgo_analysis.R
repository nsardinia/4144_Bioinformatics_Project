if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("topGO", quietly = TRUE))
  BiocManager::install("topGO", force = TRUE)

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

if (!requireNamespace("Rgraphviz", quietly = TRUE))
  BiocManager::install("Rgraphviz")

library(topGO)
library(org.Hs.eg.db) 
library(dplyr)
library(Rgraphviz)

files <- c("diff_expression_t0_t24.csv",
           "diff_expression_t0_t48.csv",
           "diff_expression_t0_t72.csv",
           "diff_expression_t24_t48.csv",
           "diff_expression_t24_t72.csv",
           "diff_expression_t48_t72.csv")

alpha <- 0.05

for (f in files) {
  deg <- read.csv(f)
  
  deg <- deg %>% mutate(isSig = ifelse(padj < alpha, 1, 0))
  
  geneList <- deg$isSig
  names(geneList) <- deg$symbol  # use gene symbols
  
  GOdata <- new("topGOdata",
                description = "GO enrichment",
                ontology = "BP",   # BP = Biological Process, MF, CC also possible
                allGenes = geneList,
                geneSel = function(x) x == 1,
                nodeSize = 10,
                annot = annFUN.org,
                mapping = "org.Hs.eg.db",
                ID = "symbol")
  
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  
  allRes <- GenTable(GOdata,
                     classicFisher = resultFisher,
                     elimFisher = resultFisher.elim,
                     orderBy = "elimFisher",
                     topNodes = 30)
  
  print(allRes)
  
  comp <- sub("diff_expression_(.*)\\.csv", "\\1", basename(f))
  
  out_file <- paste0("GO_enrichment_results_", comp, ".csv")
  write.csv(allRes, file = out_file, row.names = FALSE)
  
  # Optional: visualize
  pdf(paste0("GO_plot_", comp, ".pdf"))
  showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = "all")
  dev.off()
}

