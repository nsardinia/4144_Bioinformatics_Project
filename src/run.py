import pandas as pd
import mygene

df = pd.read_csv("data/SRP063948.tsv", sep="\t")

ensembl_ids = df.iloc[:, 0].dropna().tolist()

mg = mygene.MyGeneInfo()
results = mg.querymany(ensembl_ids,
                       scopes="ensembl.gene",
                       fields="symbol,name",
                       species="human")

mapped_df = pd.DataFrame(results)[["query", "symbol", "name"]]

final_df = df.merge(mapped_df,left_on=df.columns[0], right_on="query", how="left")

final_df.to_csv("data/with_gene_names.tsv", sep="\t", index=False)

print("Done! Saved with_gene_names.tsv")
