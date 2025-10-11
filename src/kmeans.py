import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA


def kmeans(datapath, metadatapath):

    df = pd.read_csv(datapath, sep="\t")
    meta_df = pd.read_csv(metadatapath, sep="\t", index_col=0)

    if df.columns[0].lower().startswith("ensg") or not pd.api.types.is_numeric_dtype(df.iloc[:, 0]):
        df = df.set_index(df.columns[0])

    transposed = df.T  # Transpose the data

    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0)  # Convert strings to numeric data

    # Make and color code TSNE Plot
    k = 3
    km = KMeans(n_clusters=k)
    clusters = km.fit_predict(numeric_data)
    clustered_df = transposed.copy()
    clustered_df['Cluster'] = clusters

    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(numeric_data)

    # Plot
    sns.scatterplot(x=pca_result[:, 0], y=pca_result[:, 1], hue=clusters, palette='Set2')
    plt.title(f"K-means clustering (k={k})")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.show()


kmeans("data/with_gene_names.tsv", "data/metadata_SRP120552.tsv")
