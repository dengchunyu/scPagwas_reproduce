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

method = sys.argv[1]
gwas = sys.argv[2]

DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/"
H5AD_FILE = os.path.join(DATA_PATH, "modeldata_NK_addata.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
scdrs.preprocess(adata)


df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/" + method + "/scDRS_result/" + gwas + ".geneset.gs")
df_scores=[]
for x in range(99):
    n=method + '_top' + str(x+1)
    df_score = scdrs.score_cell(
                data=adata,
                gene_list=df_gs[n][0],
                gene_weight=df_gs[n][1],
                ctrl_match_key="mean_var",
                n_ctrl=200,
                weight_opt="vs",
                return_ctrl_raw_score=False,
                return_ctrl_norm_score=True,
                verbose=False)
    df_scores.append(df_score.raw_score)

df_re=pd.concat(objs=df_scores,axis=1,ignore_index=True)
df_re.to_csv("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/" + method + "/scDRS_result/modeldata_NK/"+ gwas + "_scDRSresult.csv", sep=",")
