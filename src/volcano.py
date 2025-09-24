import bioinfokit.visuz
from diff_analysiscolor import dif_analysis
import matplotlib.pyplot as plt


def volcanoplot(datapath, metadatapah):

    da = dif_analysis(datapath, metadatapah)

    bioinfokit.visuz.GeneExpression.volcano(df=da['t24'], lfc='log2FoldChange', pv='pvalue', show=True, plotlegend=True,
                                            sign_line=True, theme="dark")
    bioinfokit.visuz.GeneExpression.volcano(df=da['t48'], lfc='log2FoldChange', pv='pvalue', show=True, plotlegend=True,
                                            sign_line=True, theme="dark")
    bioinfokit.visuz.GeneExpression.volcano(df=da['t72'], lfc='log2FoldChange', pv='pvalue', show=True, plotlegend=True,
                                            sign_line=True, theme="dark")


volcanoplot("data/SRP120552.tsv", "data/metadata_SRP120552.tsv")