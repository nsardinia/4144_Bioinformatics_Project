if (!requireNamespace("gprofiler2", quietly = TRUE))
  install.packages("gprofiler2")

library("dplyr")
library("gprofiler2")
library("tools") 

files <- c("diff_expression_t0_t24.csv",
           "diff_expression_t0_t48.csv",
           "diff_expression_t0_t72.csv",
           "diff_expression_t24_t48.csv",
           "diff_expression_t24_t72.csv",
           "diff_expression_t48_t72.csv")

alpha <- 0.05

for (f in files) {
  
  deg <- read.csv(f)
  
  sigGenes <- deg %>%
    filter(padj < alpha) %>%
    pull(symbol) %>%
    unique()
  
  gp_run <- gost(
    query = sigGenes,
    organism = "hsapiens", 
    sources = c("GO:BP"),
    correction_method = "fdr"
  )
  
   results_flat <- gp_run$result %>%
      mutate(across(where(is.list), ~sapply(., function(x) paste(x, collapse = ",")))) %>%
      arrange(p_value) %>%       
      slice_head(n=30)             
    

    top_n = 30
    baseName <- file_path_sans_ext(basename(f))
    comp <- sub("diff_expression_(.*)\\.csv", "\\1", basename(f))
    outFile <- paste0("GProf_enrichment_results_", comp, ".csv")
    write.csv(results_flat, outFile, row.names = FALSE)
}
