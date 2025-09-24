from pathlib import Path

import pandas as pd
import mygene


def create_table_with_genes(data_path: Path) -> Path:
    """Add gene names to table and save as new tsv file."""
    save_path = data_path.parent / f"{data_path.name}_with_genes.tsv"
    
    if save_path.exists():
        return save_path
    
    df = pd.read_csv(data_path, sep="\t")

    ensembl_ids = df.iloc[:, 0].dropna().tolist()

    mg = mygene.MyGeneInfo()
    results = mg.querymany(ensembl_ids,
                        scopes="ensembl.gene",
                        fields="symbol,name",
                        species="human")

    mapped_df = pd.DataFrame(results)[["query", "symbol", "name"]]

    final_df = df.merge(mapped_df,left_on=df.columns[0], right_on="query", how="left")

    cols = ["symbol", "name"] + [c for c in final_df.columns if c not in ["symbol", "name", "query"]]
    final_df = final_df[cols]


    final_df.to_csv(save_path, sep="\t", index=False)

    return save_path
