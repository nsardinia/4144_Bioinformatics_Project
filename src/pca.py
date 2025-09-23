"""
Function to generate a pca plot for our data.
"""
import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns


def pcaplot(datapath, metadatapath):

    df = pd.read_csv(datapath, sep="\t")

    if df.columns[0].lower().startswith("ensg") or not pd.api.types.is_numeric_dtype(df.iloc[:, 0]):
        df = df.set_index(df.columns[0])

    transposed = df.T  # Transpose the data

    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0)  # Convert strings to numeric data

    scaled_data = preprocessing.scale(numeric_data)  # Scale data

    pca = PCA()
    pca_data = pca.fit_transform(scaled_data)  # Fit & transform data

    # Calculate the importance of each PC (1 was the overwhelming majority of difference so comparision is between
    # PC1 & PC2
    per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
    labels = ['PC' + str(x) for x in range(1, len(per_var) + 1)]

    meta_df = pd.read_csv(metadatapath, sep="\t", index_col=0)
    pca_df = pd.DataFrame(pca_data, index=numeric_data.index, columns=labels)
    pca_df = pca_df.join(meta_df["refinebio_time"])

    plt.xlabel("PC1 = {0}%".format(per_var[0]))
    plt.ylabel("PC2 = {0}%".format(per_var[1]))

    sns.scatterplot(x='PC1', y='PC2', hue='refinebio_time', data=pca_df)
    plt.title("PCA graph by time point")
    plt.show()


# pcaplot("data/with_gene_names.tsv", "data/metadata_SRP120552.tsv")