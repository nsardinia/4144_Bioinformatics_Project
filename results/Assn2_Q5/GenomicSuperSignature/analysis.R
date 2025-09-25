library(GenomicSuperSignature)
library(data.table)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

# Load the normalized expression data
expr <- fread("gene_expression.tsv")

# Get a list of gene ids
gene_ids <- expr$Gene

# Construct the expression matrix
expr_matrix <- as.matrix(expr[, -1, with = FALSE])

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gene_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Use symbols as row names
rownames(expr_matrix) <- gene_symbols

# Remove unresolved gene ids
expr_matrix <- expr_matrix[!is.na(rownames(expr_matrix)), ]

# Download publicly available RAV model
RAVmodel <- getModel("C2", version = "latest", load = TRUE)

# Run GenomicSuperSignature validation (correlation)
val_all <- validate(expr_matrix, RAVmodel)

# Write top 50 correlated RAVs to csv
val_all_ordered <- arrange(val_all, desc(val_all$score))
top_50 <- head(val_all_ordered, 50)
write.csv(top_50, "GenomicSuperSignature_results.csv")

# Plot validation results
heatmapTable(val_all, RAVmodel = RAVmodel, num.out = 10, swCutoff = 0)
plotValidate(val_all, interactive = FALSE)

validated_ind <- validatedSignatures(val_all, RAVmodel, num.out = 3, 
                                     swCutoff = 0, indexOnly = TRUE)

set.seed(1)
drawWordcloud(RAVmodel, validated_ind[1])

