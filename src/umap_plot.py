"""
Function to generate a UMAP plot for our data.
"""

import pandas as pd
from sklearn import preprocessing
import seaborn as sns
import matplotlib.pyplot as plt
import umap


def umapplot(datapath, metadatapath):

    df = pd.read_csv(datapath, sep="\t")
    meta_df = pd.read_csv(metadatapath, sep="\t", index_col=0)
    if df.columns[0].lower().startswith("ensg") or not pd.api.types.is_numeric_dtype(df.iloc[:, 0]):
        df = df.set_index(df.columns[0])

    transposed = df.T  # Transpose the data

    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0)  # Convert strings to numeric data
    scaled_data = preprocessing.scale(numeric_data)  # Scale data

    # Create umap
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=32)
    u_features = reducer.fit_transform(scaled_data)

    umap_df = pd.DataFrame(u_features, index=numeric_data.index, columns=['x', 'y'])
    umap_df = umap_df.join(meta_df["refinebio_time"])

    sns.scatterplot(x='PC1', y='PC2', hue='refinebio_time', data=umap_df, palette='tab10')
    plt.title("UMAP plot by time point")
    plt.show()


# umapplot("data/with_gene_names.tsv", "data/metadata_SRP120552.tsv")