"""
Function to generate a t-SNE plot for our data.
"""
import pandas as pd
from sklearn.manifold import TSNE
import seaborn as sns
import matplotlib.pyplot as plt


def tsneplot(datapath, metadatapath):

    df = pd.read_csv(datapath, sep="\t")
    meta_df = pd.read_csv(metadatapath, sep="\t", index_col=0)

    if df.columns[0].lower().startswith("ensg") or not pd.api.types.is_numeric_dtype(df.iloc[:, 0]):
        df = df.set_index(df.columns[0])

    transposed = df.T  # Transpose the data

    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0)  # Convert strings to numeric data

    # Make and color code TSNE Plot
    m = TSNE(learning_rate=50, random_state=32)
    tsne_features = m.fit_transform(numeric_data)
    tsne_df = pd.DataFrame(tsne_features, columns=['x', 'y'], index=numeric_data.index)
    tsne_df = tsne_df.join(meta_df["refinebio_time"])

    sns.scatterplot(x='PC1', y='PC2', hue='refinebio_time', data=tsne_df)
    plt.title("t-SNE plot by time point")
    plt.show()


# tsneplot("data/with_gene_names.tsv", "data/metadata_SRP120552.tsv")