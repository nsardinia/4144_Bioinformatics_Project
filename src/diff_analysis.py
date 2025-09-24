from sanbomics.plots import volcano

from pathlib import Path
import pandas as pd

import seaborn as sns
import numpy as np

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def differential_analysis(data_path: Path, metadata_path: Path) -> Path:

    # Load expression matrix
    df = pd.read_csv(data_path, sep="\t")

    # Extract key columns
    gene_info = df[["symbol", "name", "Gene"]]
    expr_data = df.drop(columns=["symbol", "name"])

    counts = expr_data.T

    # Temporary until we get the counts dataset
    counts = counts.map(lambda x: x if isinstance(x,str) else int(x)).drop("Gene")


    # Load metadata
    meta = pd.read_csv(metadata_path, sep="\t")
    metadata = pd.DataFrame(
        zip(
            meta["refinebio_accession_code"], 
            meta["refinebio_time"],
        ),
        columns=["sample", "condition"]
    )
    metadata = metadata.set_index("sample")


    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design_factors="condition"
    )

    dds.deseq2()

    # Define comparisons
    comparisons = [("t72", "t0"), ("t72", "t24"), ("t72", "t48")]

    for g1, g2 in comparisons:
        stat_res = DeseqStats(dds, n_cpus=1, contrast= ("condition", g1, g2))
        stat_res.summary()
        res = stat_res.results_df

        res = res[res.baseMean >= 10]

        sigs = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > 0.5)]

        # Heatmap
        dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])

        dds_sigs = dds[:, sigs.index]

        grapher = pd.DataFrame(
            dds_sigs.layers['log1p'].T,
            index=dds_sigs.var_names,
            columns=dds_sigs.obs_names
        )

        sns.clustermap(grapher, z_score=0)