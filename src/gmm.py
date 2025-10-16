import pandas as pd
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import chi2_contingency
import seaborn as sns
from typing import Optional


def plotGMM(labels: list[int], save_path: Optional[str]="data/img/gmm_clusters.png") -> None:
    plt.figure(figsize=(8, 6))
    
    # Histogram of cluster label frequencies
    sns.countplot(x=labels)
    
    plt.title('Frequency of Cluster Labels')
    plt.xlabel('Cluster Label')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)



#Use this for generic / further computation with the same standardized dataframe.
def gmm(top_5000_df: pd.DataFrame, k:int, num_of_genes:int) -> list[int]:
    #run preliminary PCA to reduce dimensionality before GMM (otherwise, too large)
    if num_of_genes < 500:
        pca = PCA(n_components=num_of_genes, random_state=42)
    else:
        pca = PCA(n_components=500, random_state=42)
    X_reduced = pca.fit_transform(top_5000_df)


    gmodel = GaussianMixture(n_components=k, random_state=42)
    gmodel.fit(X_reduced)
    labels = gmodel.predict(X_reduced)

    return labels

def processDF(num_of_genes: int) -> None:
    df = pd.read_csv("data/with_gene_names.tsv", sep="\t")
    if df.columns[0].lower().startswith("ensg") or not pd.api.types.is_numeric_dtype(df.iloc[:, 0]):
        df = df.set_index(df.columns[0])

    transposed = df.T  

    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0)  
    gene_var = numeric_data.var(axis=1)
    genes = gene_var.sort_values(ascending=False).head(num_of_genes).index
    data = numeric_data.loc[genes]

    return data

# #g_10 = plotGMM(gmm(pd.read_csv("data/with_gene_names.tsv", sep="\t"), 2, 10), "data/img/gmm_frq_10_2.png")
# g_100 = plotGMM(gmm(processDF(100), 2, 100), "data/img/gmm_frq_100_2.png")
# g_1000 = plotGMM(gmm(processDF(1000), 2, 1000), "data/img/gmm_frq_1000_2.png")
# g_5000 = plotGMM(gmm(processDF(5000), 2, 5000), "data/img/gmm_frq_5000_2.png")
# g_10000 = plotGMM(gmm(processDF(10000), 2, 10000), "data/img/gmm_frq_10000_2.png")

# #g_10 = plotGMM(gmm(pd.read_csv("data/with_gene_names.tsv", sep="\t"), 4, 10), "data/img/gmm_frq_10_4.png")
# g_100 = plotGMM(gmm(processDF(100), 4, 100), "data/img/gmm_frq_100_4.png")
# g_1000 = plotGMM(gmm(processDF(1000), 4, 1000), "data/img/gmm_frq_1000_4.png")
# g_5000 = plotGMM(gmm(processDF(5000), 4, 5000), "data/img/gmm_frq_5000_4.png")
# g_10000 = plotGMM(gmm(processDF(10000), 4, 10000), "data/img/gmm_frq_10000_4.png")

# #g_10 = plotGMM(gmm(pd.read_csv("data/with_gene_names.tsv", sep="\t"), 8, 10), "data/img/gmm_frq_10_8.png")
# g_100 = plotGMM(gmm(processDF(100), 8, 100), "data/img/gmm_frq_100_8.png")
# g_1000 = plotGMM(gmm(processDF(1000), 8, 1000), "data/img/gmm_frq_1000_8.png")
# g_5000 = plotGMM(gmm(processDF(5000), 8, 5000), "data/img/gmm_frq_5000_8.png")
# g_10000 = plotGMM(gmm(processDF(10000), 8, 10000), "data/img/gmm_frq_10000_8.png")

# #g_10 = plotGMM(gmm(pd.read_csv("data/with_gene_names.tsv", sep="\t"), 16, 10), "data/img/gmm_frq_10_16.png")
# g_100 = plotGMM(gmm(processDF(100), 16, 100), "data/img/gmm_frq_100_16.png")
# g_1000 = plotGMM(gmm(processDF(1000), 16, 1000), "data/img/gmm_frq_1000_16.png")
# g_5000 = plotGMM(gmm(processDF(5000), 16, 5000), "data/img/gmm_frq_5000_16.png")
# g_10000 = plotGMM(gmm(processDF(10000), 16, 10000), "data/img/gmm_frq_10000_16.png")

# #g_10 = plotGMM(gmm(pd.read_csv("data/with_gene_names.tsv", sep="\t"), 32, 10), "data/img/gmm_frq_10_32.png")
# g_100 = plotGMM(gmm(processDF(100), 32, 100), "data/img/gmm_frq_100_32.png")
# g_1000 = plotGMM(gmm(processDF(1000), 32, 1000), "data/img/gmm_frq_1000_32.png")
# g_5000 = plotGMM(gmm(processDF(5000), 32, 5000), "data/img/gmm_frq_5000_32.png")
# g_10000 = plotGMM(gmm(processDF(10000), 32, 10000), "data/img/gmm_frq_10000_32.png")

