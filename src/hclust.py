import pandas as pd
from sklearn.cluster import AgglomerativeClustering

def hclust(top_5000_df: pd.DataFrame, k:int) -> list[int]:
    """
    Return list of labels for all samples.
    """
    agg_clustering = AgglomerativeClustering(
                        n_clusters=k, 
                        affinity='euclidean', 
                        linkage='ward',
                    )
    
    clustering_labels = agg_clustering.fit_predict(top_5000_df)

    return clustering_labels