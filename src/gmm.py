import pandas as pd
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import chi2_contingency
import seaborn as sns



def gmm(datapath, k, num_of_genes):
    df = pd.read_csv(datapath, sep="\t")

    if df.columns[0].lower().startswith("ensg") or not pd.api.types.is_numeric_dtype(df.iloc[:, 0]):
        df = df.set_index(df.columns[0])

    transposed = df.T  

    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0)  
    gene_var = numeric_data.var(axis=1)
    genes = gene_var.sort_values(ascending=False).head(num_of_genes).index
    data = numeric_data.loc[genes]

    pca = PCA(n_components=50, random_state=42)
    X_reduced = pca.fit_transform(data)

    gmodel = GaussianMixture(n_components=k, random_state=42)
    gmodel.fit(X_reduced)
    labels = gmodel.predict(X_reduced)

    pca_vis = PCA(n_components=2, random_state=0)
    X_2D = pca_vis.fit_transform(X_reduced)

    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=X_2D[:, 0], y=X_2D[:, 1], hue=labels, palette='Set1', s=60, edgecolor='k')
    plt.title('GMM Clustering (Visualized in 2D via PCA)')
    plt.xlabel('PCA Component 1')
    plt.ylabel('PCA Component 2')
    plt.legend(title='Cluster')
    plt.tight_layout()
    plt.savefig("data/img/gmm_clusters10000.png", dpi=300)


#g_10 = gmm("data/with_gene_names.tsv", 8, 10)
#g_100 = gmm("data/with_gene_names.tsv", 8, 100)
#g_1000 = gmm("data/with_gene_names.tsv", 1000)
#g_5000 = gmm("data/with_gene_names.tsv", 8, 5000)
g_10000 = gmm("data/with_gene_names.tsv", 8, 10000)


#gmm() -> 