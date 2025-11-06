import numpy as np
import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, roc_auc_score
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt
from gmm import gmm
from gmm import processDF
import seaborn as sns

# Create path if does not exist
def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def _plot_schart(y_test, y_prob, title, label_name, outpath, figsize=(8, 5)):
    """Helper function to create and save S-chart plots."""
    sorted_idx = np.argsort(y_prob)
    y_prob_sorted = y_prob[sorted_idx]
    y_test_sorted = y_test.iloc[sorted_idx] if isinstance(y_test, pd.Series) else y_test[sorted_idx]
    
    plt.figure(figsize=figsize)
    plt.scatter(range(len(y_test_sorted)), y_prob_sorted, c=y_test_sorted, cmap='bwr', alpha=0.6)
    plt.axhline(0.5, color='gray', linestyle='--')
    plt.xlabel("Samples sorted by predicted probability")
    plt.ylabel("Predicted Probability")
    plt.title(title)
    plt.colorbar(label=label_name)
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight", dpi=300)
    plt.close()


def _plot_confusion_matrix(cm, title, outpath, figsize=(6, 5)):
    """Helper function to create and save confusion matrix heatmaps."""
    plt.figure(figsize=figsize)
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
    plt.title(title)
    plt.xlabel("Predicted Cluster")
    plt.ylabel("True Cluster")
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight", dpi=300)
    plt.close()

def svm_classification(datapath, metadata, num_of_genes, kernel='rbf', C=1.0, results_dir="../results/Assn4/SVM"):
    ensure_dir(results_dir)

    expression_df = pd.read_csv(datapath, sep="\t", index_col=0)
    metadata = pd.read_csv(metadata, sep="\t")

    expression_df = expression_df.apply(pd.to_numeric, errors="coerce")
    gene_var = expression_df.var(axis=1)
    genes = gene_var.sort_values(ascending=False).head(num_of_genes).index
    expression_df = expression_df.loc[genes]
    expr_log = expression_df.apply(lambda x: np.log1p(x))
    score = expr_log.mean()
    metadata["reprog_score"] = metadata["refinebio_accession_code"].map(score)

    threshold = np.median(metadata["reprog_score"])
    metadata["label"] = (metadata["reprog_score"] > threshold).astype(int)

    train_data = metadata.query("refinebio_time in ['t0','t24','t48']")

    X = expr_log.loc[:, train_data["refinebio_accession_code"]].T
    y = train_data["label"]

    # Ensure column names are strings to avoid sklearn errors
    X.columns = X.columns.astype(str)
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, train_size=0.8, random_state=42)

    model = SVC(kernel=kernel, C=C, probability=True, random_state=42)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]
    
    # Calculate metrics once
    cm = confusion_matrix(y_test, y_pred)
    auc = roc_auc_score(y_test, y_prob)
    acc = accuracy_score(y_test, y_pred)

    print(f"SVM (kernel={kernel}, C={C})")
    print(f"AUC: {auc:.3f}")
    print(f"Accuracy: {acc:.3f}")
    print("Confusion Matrix:\n", cm)
    print("Classification Report:\n", classification_report(y_test, y_pred))

    # Use helper function for S-chart
    outpath = os.path.join(results_dir, f"SVM_base_{kernel}_C{C}_genes{num_of_genes}.png")
    _plot_schart(y_test, y_prob, 
                 f"SVM S-chart (kernel={kernel}, C={C}): High vs Low Reprogramming Score",
                 "True Label (0=Low, 1=High)", 
                 outpath)
    print(f"Saved plot to: {outpath}\n")


def svm_on_gmm(datapath, num_of_genes, kernel='rbf', C=1.0, n_clusters=8, results_dir="../results/Assn4/SVM"):
    ensure_dir(results_dir)

    # Use processDF to get clusters - it already loads and processes the data
    data_5000 = processDF(num_of_genes)
    # Ensure column names are strings to avoid PCA errors
    data_5000.columns = data_5000.columns.astype(str)
    clusters = pd.Series(gmm(data_5000, n_clusters, num_of_genes), index=data_5000.index)

    # Reuse the processed data from processDF instead of reloading
    # data_5000 already has samples as rows and genes as columns (transposed format)
    # We just need to ensure column names are strings and scale
    X = data_5000.copy()
    X.columns = X.columns.astype(str)
    X_scaled = StandardScaler().fit_transform(X)

    # Align indices - ensure X_scaled and clusters match
    common_idx = X.index.intersection(clusters.index)
    # Filter X_scaled using boolean indexing
    mask = X.index.isin(common_idx)
    X_scaled = X_scaled[mask]
    y = clusters.loc[common_idx].values

    # Binary classification: one cluster vs all others
    # Find the cluster with the most samples for better balance
    cluster_counts = pd.Series(y).value_counts()
    target_cluster = cluster_counts.index[0]  # Use the largest cluster
    y_binary = (y == target_cluster).astype(int)

    try:
        X_train, X_test, y_train, y_test = train_test_split(
            X_scaled, y_binary, train_size=0.8, random_state=42, stratify=y_binary
        )
    except ValueError:
        X_train, X_test, y_train, y_test = train_test_split(
            X_scaled, y_binary, train_size=0.8, random_state=42, stratify=None
        )

    model = SVC(kernel=kernel, C=C, probability=True, class_weight='balanced', random_state=42)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]

    # Calculate metrics once
    cm = confusion_matrix(y_test, y_pred)
    acc = accuracy_score(y_test, y_pred)
    auc = roc_auc_score(y_test, y_prob)
    
    print(f"\nSVM (kernel={kernel}, C={C}) attempting to match GMM clustering")
    print(f"Number of clusters: {n_clusters}")
    print(f"Target cluster: {target_cluster} (largest cluster with {cluster_counts[target_cluster]} samples)")
    print(f"Accuracy: {acc:.3f}")
    print(f"AUC: {auc:.3f}")
    print("Confusion Matrix:\n", cm)
    print("Classification Report:\n", classification_report(y_test, y_pred))

    # Use helper functions for plotting
    cm_outpath = os.path.join(results_dir, f"SVM_vs_GMM_confusion_matrix_{kernel}_C{C}_clusters{n_clusters}_genes{num_of_genes}.png")
    _plot_confusion_matrix(cm, f"SVM vs GMM Clustering (kernel={kernel}, C={C}, clusters={n_clusters})", cm_outpath)
    print(f"Saved confusion matrix to: {cm_outpath}")

    schart_outpath = os.path.join(results_dir, f"SVM_vs_GMM_schart_{kernel}_C{C}_clusters{n_clusters}_genes{num_of_genes}.png")
    _plot_schart(y_test, y_prob, 
                 f"S-chart: Cluster {target_cluster} vs All Other Clusters",
                 f"True Label (0=Other, 1=Cluster {target_cluster})",
                 schart_outpath)
    print(f"Saved S-chart to: {schart_outpath}\n")


# Example usage
if __name__ == "__main__":
    # Run with different numbers of genes
    gene_counts = [10, 100, 1000, 5000]
    
    for num_genes in gene_counts:
        print(f"\n{'='*80}")
        print(f"Running SVM with {num_genes} genes")
        print(f"{'='*80}\n")
        
        # Base classification (Assignment 1 groups)
        svm_classification("data/with_gene_names.tsv", "data/metadata_SRP120552.tsv", 
                          num_genes, kernel='rbf', C=1.0)
        
        # Cluster prediction (Assignment 3)
        svm_on_gmm("data/with_gene_names.tsv", num_genes, kernel='rbf', C=1.0)
    
    print(f"\n{'='*80}")
    print("All runs completed!")
    print(f"{'='*80}\n")

