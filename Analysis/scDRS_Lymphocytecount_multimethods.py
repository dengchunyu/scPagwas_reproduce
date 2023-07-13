import scdrs
from scipy import stats
import pandas as pd
import scanpy as sc
from anndata import AnnData
sc.set_figure_params(dpi=125)
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
import sys


methods=['smultixcan', 'spredixcan' ,'TWAS']
gwass=['eosinophilcount', 'basophilcount','LymphocytePercent','monocytecount' ,'neutrophilcount', 'WhiteBloodCellcount', 'Hemoglobinconcen', 'MeanCorpuscularHemoglobin', 'MeanCorpusVolume']

gwass=['Lymphocytecount3']
for method in methods:
    for gwas in gwass:
        scDRS_file= "/share/pub/dengcy/GWAS_Multiomics/MultiMethods/" + method + "/" + gwas + ".scDRS.multiMethods.genes.csv"
        df_gs=pd.read_csv(scDRS_file, index_col=0)
        df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/" + method + "/scDRS_result/" + gwas +".geneset.gs", sep="\t", index=False)

