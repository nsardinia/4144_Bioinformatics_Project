import bioinfokit.visuz
from diff_analysis import dif_analysis


def volcanoplot(datapath, metadatapah):

    da = dif_analysis(datapath, metadatapah)

    bioinfokit.visuz.GeneExpression.volcano(df=da['t24'], lfc='log2FoldChange', pv='pvalue', show=True)


volcanoplot("data/SRP120552.tsv", "data/metadata_SRP120552.tsv")