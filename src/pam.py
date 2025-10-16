import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn_extra.cluster import KMedoids
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import chi2_contingency
from sklearn.metrics import silhouette_score

# Find optimal number of clusters
def optimize_pam(datapath, max_k):
    means = []
    inertias = []

    df = pd.read_csv(datapath, sep="\t")
    expression_data = df.iloc[:, 2:]
    expression_data = expression_data.set_index(expression_data.columns[0])
    
    gene_var = expression_data.var(axis=1)
    genes = gene_var.sort_values(ascending=False).head(5000).index
    selected_genes_data = expression_data.loc[genes]
    
    transposed = selected_genes_data.T
    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0)
    numeric_data.columns = numeric_data.columns.astype(str)
    scaled_data = StandardScaler().fit_transform(numeric_data)
    for k in range(1, max_k):
        pam = KMedoids(n_clusters=k, random_state=42)
        pam.fit(scaled_data)
        means.append(k)
        inertias.append(pam.inertia_ if hasattr(pam, "inertia_") else (pam.transform(scaled_data).min(axis=1).sum()))
    fig, ax = plt.subplots(figsize=(10,5))
    ax.plot(means, inertias, 'o-')
    ax.set_xlabel('Numbers of Clusters')
    ax.set_ylabel('Sum of Distances to Medoids')
    ax.grid(True)
    plt.savefig('../results/Assn3_Q2/PAM/pam_optimization.png', dpi=300, bbox_inches='tight')
    plt.close()


def compare_k_values(datapath, num_of_genes, max_k=10):
    df = pd.read_csv(datapath, sep="\t")
    expression_data = df.iloc[:, 2:]
    expression_data = expression_data.set_index(expression_data.columns[0])
    
    gene_var = expression_data.var(axis=1)
    genes = gene_var.sort_values(ascending=False).head(num_of_genes).index
    selected_genes_data = expression_data.loc[genes]
    
    transposed = selected_genes_data.T
    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0)
    numeric_data.columns = numeric_data.columns.astype(str)
    scaled_data = StandardScaler().fit_transform(numeric_data)
    k_values = range(2, max_k + 1)
    inertias = []
    
    for k in k_values:
        pam = KMedoids(n_clusters=k, random_state=42)
        clusters = pam.fit_predict(scaled_data)
        
        inertia = pam.inertia_ if hasattr(pam, 'inertia_') else pam.transform(scaled_data).min(axis=1).sum()
        inertias.append(inertia)
        print(f"k={k}: inertia={inertia:.2f}")
    plt.figure(figsize=(10, 6))
    plt.plot(k_values, inertias, 'o-')
    plt.xlabel('Number of Clusters (k)')
    plt.ylabel('Inertia')
    plt.title('PAM Clustering: Elbow Method')
    plt.grid(True)
    plt.savefig('../results/Assn3_Q2/PAM/k_value_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    inertia_diffs = np.diff(inertias)
    elbow_idx = np.argmin(inertia_diffs) + 1
    best_k_elbow = k_values[elbow_idx]
    print(f"\nBest k by elbow method: {best_k_elbow}")
    
    return best_k_elbow, inertias

# Return list of labels for all samples.
def pam_for_heatmap(top_5000_df: pd.DataFrame, k: int) -> list[int]:
    pam = KMedoids(n_clusters=k, random_state=42)
    clustering_labels = pam.fit_predict(top_5000_df)
    return clustering_labels

# PAM clustering with PCA visualization
def pam(datapath, k, num_of_genes):
    df = pd.read_csv(datapath, sep="\t")
    expression_data = df.iloc[:, 2:]
    expression_data = expression_data.set_index(expression_data.columns[0])
    
    gene_var = expression_data.var(axis=1)
    genes = gene_var.sort_values(ascending=False).head(num_of_genes).index
    selected_genes_data = expression_data.loc[genes]
    numeric_data = selected_genes_data.apply(pd.to_numeric, errors="coerce").fillna(0)
    scaled_data = StandardScaler().fit_transform(numeric_data.T).T

    pam = KMedoids(n_clusters=k, random_state=42)
    clusters = pam.fit_predict(scaled_data)

    clustered_df = numeric_data.copy()
    clustered_df["Cluster"] = clusters

    pca = PCA(n_components=2, random_state=42)
    pca_result = pca.fit_transform(scaled_data)

    fig, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=pca_result[:, 0], y=pca_result[:, 1], hue=clusters, palette="Set2", ax=ax)
    ax.set_title(f"PAM Gene Clustering (k={k}, {num_of_genes} genes)")
    ax.set_xlabel("PCA 1")
    ax.set_ylabel("PCA 2")
    plt.savefig(f'../results/Assn3_Q2/PAM/pam_{num_of_genes}.png', dpi=300, bbox_inches='tight')
    plt.close()

    return pd.Series(clusters, index=numeric_data.index, name="Cluster")


# Chi-squared test for comparison
def chi(g10, g100, g1000, g5000, g10000):
    results = []
    cluster_sets = {'10': g10, '100': g100, '1000': g1000, '5000': g5000, '10000': g10000}
    pairs = [
        ('10', '100'), ('10', '1000'), ('10', '5000'), ('10', '10000'),
        ('100', '1000'), ('100', '5000'), ('100', '10000'),
        ('1000', '5000'), ('1000', '10000'), ('5000', '10000')
    ]

    for a, b in pairs:
        gA, gB = cluster_sets[a], cluster_sets[b]
        common = gA.index.intersection(gB.index)
        contingency = pd.crosstab(gA.loc[common], gB.loc[common])
        chi2, p, dof, expected = chi2_contingency(contingency)
        results.append({
            'Comparison': f'{a} vs {b}',
            'Chi²': chi2,
            'p-value': p,
            'DOF': dof,
            'Common Genes': len(common)
        })

    results_df = pd.DataFrame(results)
    print("\nChi-squared comparison between PAM clusterings:\n")
    print(results_df.to_string(index=False))
    
    results_df.to_csv('../results/Assn3_Q2/PAM/pam_chi_square_results.csv', index=False)
    print(f"\nChi-square results saved to: ../results/Assn3_Q2/PAM/pam_chi_square_results.csv")


if __name__ == "__main__":
    g_10 = pam("data/with_gene_names.tsv", 4, 10)
    g_100 = pam("data/with_gene_names.tsv", 4, 100)
    g_1000 = pam("data/with_gene_names.tsv", 4, 1000)
    g_5000 = pam("data/with_gene_names.tsv", 4, 5000)
    g_10000 = pam("data/with_gene_names.tsv", 4, 10000)

    chi(g_10, g_100, g_1000, g_5000, g_10000)
