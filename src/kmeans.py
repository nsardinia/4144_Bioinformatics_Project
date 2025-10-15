import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.stats import chi2_contingency
from sklearn.preprocessing import StandardScaler

def optimize(datapath, max_k):
    means = []
    inertias = []
    df = pd.read_csv(datapath, sep="\t")

    if df.columns[0].lower().startswith("ensg") or not pd.api.types.is_numeric_dtype(df.iloc[:, 0]):
        df = df.set_index(df.columns[0])

    transposed = df.T  # Transpose the data

    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0)  # Convert strings to numeric data
    gene_var = numeric_data.var(axis=0)
    genes = gene_var.sort_values(ascending=False).head(5000).index
    data = numeric_data[genes]
    for k in range(1, max_k):
        kmeans = KMeans(n_clusters=k)
        kmeans.fit(data)
        means.append(k)
        inertias.append(kmeans.inertia_)
    fig = plt.subplots(figsize=(10,5))
    plt.plot(means, inertias, 'o-')
    plt.xlabel('Numbers of Clusters')
    plt.ylabel('Inertia')
    plt.grid(True)
    plt.show()


def kmeans_for_heatmap(datapath, k, num_of_genes) -> list[int]:

    df = pd.read_csv(datapath, sep="\t")

    if df.columns[0].lower().startswith("ensg") or not pd.api.types.is_numeric_dtype(df.iloc[:, 0]):
        df = df.set_index(df.columns[0])

    transposed = df.T  # Transpose the data

    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0)  # Convert strings to numeric data
    gene_var = numeric_data.var(axis=1)
    genes = gene_var.sort_values(ascending=False).head(num_of_genes).index
    data = numeric_data.loc[genes]
    scaled_data = StandardScaler().fit_transform(data)
    km = KMeans(n_clusters=k, random_state=42)
    clusters = km.fit_predict(scaled_data)
    clustered_df = data.copy()
    clustered_df['Cluster'] = clusters

    pca = PCA(n_components=2, random_state=42)
    pca_result = pca.fit_transform(scaled_data)

    # Plot
    sns.scatterplot(x=pca_result[:, 0], y=pca_result[:, 1], hue=clusters, palette='Set2')
    plt.title(f"K-means clustering (k={k})")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.show()
    return clusters
# Optimize("data/with_gene_names.tsv", 100)
# Optimal seems to be around k = 8, as the steepest drops in the inertia curve occur around k = 5 through k = 8.
# Afterward the inertia values get smaller in about a linear rate.


def kmeans(datapath, k, num_of_genes):

    df = pd.read_csv(datapath, sep="\t")

    if df.columns[0].lower().startswith("ensg") or not pd.api.types.is_numeric_dtype(df.iloc[:, 0]):
        df = df.set_index(df.columns[0])

    transposed = df.T  # Transpose the data

    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0)  # Convert strings to numeric data
    gene_var = numeric_data.var(axis=1)
    genes = gene_var.sort_values(ascending=False).head(num_of_genes).index
    data = numeric_data.loc[genes]
    scaled_data = StandardScaler().fit_transform(data)
    km = KMeans(n_clusters=k, random_state=42)
    clusters = km.fit_predict(scaled_data)
    clustered_df = data.copy()
    clustered_df['Cluster'] = clusters

    pca = PCA(n_components=2, random_state=42)
    pca_result = pca.fit_transform(scaled_data)

    # Plot
    sns.scatterplot(x=pca_result[:, 0], y=pca_result[:, 1], hue=clusters, palette='Set2')
    plt.title(f"K-means clustering (k={k})")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.show()

    return pd.Series(clusters, index=data.index, name='Cluster')


# Chi Squared
def chi(g10, g100, g1000, g5000, g10000):
    results = []
    cluster_sets = {'10': g10, '100': g100, '1000': g1000, '5000': g5000, '10000': g10000}
    pairs = [('10', '100'), ('10', '1000'), ('10', '5000'), ('10', '10000'), ('100', '1000'), ('100', '5000'),
             ('100', '10000'), ('1000', '10000'), ('1000', '5000'), ('5000', '10000')]

    for a, b in pairs:
        gA, gB = cluster_sets[a], cluster_sets[b]
        common = gA.index.intersection(gB.index)
        contingency = pd.crosstab(gA.loc[common], gB.loc[common])
        chi2, p, dof, expected = chi2_contingency(contingency)
        results.append({'Comparison': f'{a} vs {b}', 'Chi²': chi2, 'p-value': p, 'DOF': dof, 'Common Genes': len(common)})

    results_df = pd.DataFrame(results)
    print("\nChi-squared comparison between clusterings:\n")
    print(results_df.to_string(index=False))


g_10 = kmeans("data/with_gene_names.tsv", 8, 10)
g_100 = kmeans("data/with_gene_names.tsv", 8, 100)
g_1000 = kmeans("data/with_gene_names.tsv", 8, 1000)
g_5000 = kmeans("data/with_gene_names.tsv", 8, 5000)
g_10000 = kmeans("data/with_gene_names.tsv", 8, 10000)

chi(g_10, g_100, g_1000, g_5000, g_10000)

"""
Comparison        Chi²  p-value  DOF  Common Genes
    10 vs 100   30.000000 0.091988   21            10
   10 vs 1000   10.000000 0.188573    7            10
   10 vs 5000   10.000000 0.188573    7            10
  10 vs 10000   10.000000 0.188573    7            10
  100 vs 1000   13.657771 0.057614    7           100
  100 vs 5000   13.657771 0.057614    7           100
 100 vs 10000   13.657771 0.057614    7           100
1000 vs 10000 4613.000000 0.000000   49           659
 1000 vs 5000 4613.000000 0.000000   49           659
5000 vs 10000 4613.000000 0.000000   49           659

As the number of genes among the comparisons increased the amount of shared genes also increased. An interesting
observation occurred where The top 1000 genes did not correlate 1:1 to the top 5000, or 10000 genes. This shows
that as the number of genes in clustering increased new information was gained causing the algorithm to reach a
different result.
"""