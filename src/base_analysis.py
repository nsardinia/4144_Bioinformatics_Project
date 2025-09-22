"""
Functions that calculate the address questions (median, # of genes, variation, generate density plot) 
asked in Part 1 c of Assignment 2.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("data/with_gene_names.tsv", sep="\t")

log_scaled = np.log2(df.select_dtypes(include=[np.number]) + 1)

median_expression_range = log_scaled.median(axis=1, numeric_only=True)

median_expression_range.plot.density()
plt.xlim(0, 5)
plt.show()

"""
The data sets includes 43405 genes from 658 samples with little variance. The data had a lower variance of
about 0.85 with the total range of log scaled data spanning from about 0.45 to 7.58. A vast majority of the
data is centralized between 0 and 1 (about 80%) with some outliers skewing the data towards the right.
"""