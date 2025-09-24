"""
Function to run a differential analysis on our data
"""
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def dif_analysis(datapath, metadatapath):
    df = pd.read_csv(datapath, sep="\t")
    meta_df = pd.read_csv(metadatapath, sep="\t", index_col=0)

    if df.columns[0].lower().startswith("ensg") or not pd.api.types.is_numeric_dtype(df.iloc[:, 0]):
        df = df.set_index(df.columns[0])

    transposed = df.T  # Transpose the data

    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)  # Convert strings to numeric data
    numeric_data = numeric_data.loc[:, numeric_data.sum(axis=0) > 0]

    dds = DeseqDataSet(counts=numeric_data, metadata=meta_df, design_factors="refinebio_time")
    dds.deseq2()

    timepoints = [tp for tp in meta_df["refinebio_time"].unique() if tp != "t0"]
    results = {}
    for tp in timepoints:

        stats_res = DeseqStats(dds, n_cpus=1, contrast=("refinebio_time", tp, "t0"))
        stats_res.summary()

        df_res = stats_res.results_df.copy()
        df_res = df_res[pd.to_numeric(df_res['log2FoldChange'], errors='coerce').notnull()]
        df_res = df_res[pd.to_numeric(df_res['pvalue'], errors='coerce').notnull()]
        df_res['gene'] = df_res.index
        results[tp] = df_res

    return results
# Differential analysis comparing t0 to all other times.


# da = dif_analysis("data/SRP120552.tsv", "data/metadata_SRP120552.tsv")