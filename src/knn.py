import numpy as np
import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, roc_auc_score
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt
from gmm import gmm
from gmm import processDF
import seaborn as sns

#Create path if does not exist
def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def knn_classification(datapath, metadata, num_of_genes, n_neighbors=8, results_dir="../results/Assn4/KNN"):
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

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, train_size=0.8, random_state=42)

    model = KNeighborsClassifier(n_neighbors=n_neighbors)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, y_prob)

    print(f"KNN (k={n_neighbors})")
    print(f"AUC: {auc:.3f}")
    print(f"Accuracy: {accuracy_score(y_test, y_pred)}")
    print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))
    print("Classification Report:\n", classification_report(y_test, y_pred))

    sorted_idx = np.argsort(y_prob)
    y_prob_sorted = y_prob[sorted_idx]
    y_test_sorted = y_test.iloc[sorted_idx] if isinstance(y_test, pd.Series) else y_test[sorted_idx]

    plt.figure(figsize=(8, 5))
    plt.scatter(range(len(y_test_sorted)), y_prob_sorted, c=y_test_sorted, cmap='bwr', alpha=0.6)
    plt.axhline(0.5, color='gray', linestyle='--')
    plt.xlabel("Samples sorted by predicted probability")
    plt.ylabel("Predicted Probability")
    plt.title(f"KNN S-chart (k={n_neighbors}): High vs Low Reprogramming Score")
    plt.colorbar(label="True Label (0=Low, 1=High)")

    outpath = os.path.join(results_dir, f"KNN_base_k{n_neighbors}_genes{num_of_genes}.png")
    plt.savefig(outpath, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved plot to: {outpath}\n")


def knn_on_gmm(datapath, num_of_genes, n_neighbors=8, results_dir="../results/Assn4/KNN"):
    ensure_dir(results_dir)

    data_5000 = processDF(num_of_genes)
    clusters = pd.Series(gmm(data_5000, 8, num_of_genes), index=data_5000.index)

    # Load full dataset
    df = pd.read_csv(datapath, sep="\t", index_col=0)
    if df.columns[0].lower().startswith("ensg") or not pd.api.types.is_numeric_dtype(df.iloc[:, 0]):
        df = df.set_index(df.columns[0])

    # Transpose so samples are rows
    transposed = df.T
    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0)

    # Select top variable genes
    gene_var = numeric_data.var(axis=0)
    genes = gene_var.sort_values(ascending=False).head(num_of_genes).index
    X = numeric_data.loc[:, genes]
    X_scaled = StandardScaler().fit_transform(X)

    common_idx = X.index.intersection(clusters.index)
    X_scaled = pd.DataFrame(X_scaled, index=X.index).loc[common_idx].values
    y = clusters.loc[common_idx].values

    #X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, train_size=0.8, random_state=42, stratify=y)
    try:
        X_train, X_test, y_train, y_test = train_test_split(
            X_scaled, y, train_size=0.8, random_state=42, stratify=y
        )
    except ValueError:
        X_train, X_test, y_train, y_test = train_test_split(
            X_scaled, y, train_size=0.8, random_state=42, stratify=None
        )

    model = KNeighborsClassifier(n_neighbors=n_neighbors)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    acc = accuracy_score(y_test, y_pred)
    print(f"\nKNN (k={n_neighbors}) attempting to match GMM clustering")
    print(f"Accuracy: {acc:.3f}")
    print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))
    print("Classification Report:\n", classification_report(y_test, y_pred))

    plt.figure(figsize=(6, 5))
    sns.heatmap(confusion_matrix(y_test, y_pred), annot=True, fmt="d", cmap="Blues")
    plt.title(f"KNN vs GMM Clustering (k={n_neighbors})")
    plt.xlabel("Predicted Cluster")
    plt.ylabel("True Cluster")

    outpath = os.path.join(results_dir, f"KNN_vs_GMM_k{n_neighbors}_genes{num_of_genes}.png")
    plt.savefig(outpath, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved confusion matrix to: {outpath}\n")



# Example usage


knn_classification("data/with_gene_names.tsv", "data/metadata_SRP120552.tsv", 10, n_neighbors=8)
knn_on_gmm("data/with_gene_names.tsv", 10, n_neighbors=8)

#knn_classification("data/with_gene_names.tsv", "data/metadata_SRP120552.tsv", 100, n_neighbors=8)
#knn_on_gmm("data/with_gene_names.tsv", 100, n_neighbors=8)


#knn_classification("data/with_gene_names.tsv", "data/metadata_SRP120552.tsv",1000, n_neighbors=8)
#knn_on_gmm("data/with_gene_names.tsv", 1000, n_neighbors=8)

#knn_classification("data/with_gene_names.tsv", "data/metadata_SRP120552.tsv", 5000, n_neighbors=8)
#knn_on_gmm("data/with_gene_names.tsv", 5000, n_neighbors=8)

#knn_classification("data/with_gene_names.tsv", "data/metadata_SRP120552.tsv", 10000, n_neighbors=8)
#knn_on_gmm("data/with_gene_names.tsv", 10000, n_neighbors=8)
