import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, roc_auc_score
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt
from kmeans import kmeans


def logical_regression(datapath, metadata, num_of_genes):
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

    train_data = metadata.query("refinebio_time in ['t0','t24', 't48']")

    X = expr_log.loc[:, train_data["refinebio_accession_code"]].T
    y = train_data["label"]
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.8, random_state=42)

    model = LogisticRegression()
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, y_prob)

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
    plt.title("S-chart: High vs Low Reprogramming Score")
    plt.colorbar(label="True Label (0=Low, 1=High)")
    plt.show()


def regression_on_kmeans(datapath, num_of_genes):
    clusters = kmeans(datapath, 8, 5000)
    df = pd.read_csv(datapath, sep="\t", index_col=0)
    if df.columns[0].lower().startswith("ensg") or not pd.api.types.is_numeric_dtype(df.iloc[:, 0]):
        df = df.set_index(df.columns[0])

    transposed = df.T
    numeric_data = transposed.apply(pd.to_numeric, errors="coerce").fillna(0)

    gene_var = numeric_data.var(axis=0)
    genes = gene_var.sort_values(ascending=False).head(num_of_genes).index
    X = numeric_data.loc[:, genes]
    X_scaled = StandardScaler().fit_transform(X)

    target_cluster = 4
    y = (clusters.loc[X.index] == target_cluster).astype(int)

    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, train_size=0.8, random_state=42)

    model = LogisticRegression(max_iter=1000, class_weight="balanced")
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, y_prob)

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
    plt.title(f"S-chart: Cluster {target_cluster} vs All Other Clusters")
    plt.colorbar(label="True Label (0=Other, 1=Cluster 4)")
    plt.show()


logical_regression("data/with_gene_names.tsv", "data/metadata_SRP120552.tsv", 10000)
regression_on_kmeans("data/with_gene_names.tsv", 10000)

"""
A binary regression where cluster 4 is compared against all the other clusters was used. This was done because it 
produced the best outcome and sigmoid curve. Cluster 4 has a large sample size which helped a lot as several minority
clusters exist that have a little as one member causing the regression to give bad results.
"""

""""
AUC / Accuracy:
BASE DATA
10: 1.0 / 0.99
100: 0.90 / 0.84
1000: 0.93 / 0.80
5000: 0.99 / 0.89
10000: 0.99 / 0.88

KMEANS
10: 0.72 / 0.65
100: 0.79 / 0.72
1000: 0.82 / 0.79
5000: 0.84 / 0.80
10000: 0.86 / 0.83
"""
