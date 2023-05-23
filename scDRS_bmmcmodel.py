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


DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/realgroundtruth"
H5AD_FILE = os.path.join(DATA_PATH, "Pagwas_groundtruth_addata.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
#sc.tl.pca(adata, svd_solver="arpack")
#sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
scdrs.preprocess(adata)

methods=['smultixcan', 'spredixcan' ,'TWAS']
#gwass=['eosinophilcount', 'basophilcount','LymphocytePercent','monocytecount' ,'neutrophilcount', 'WhiteBloodCellcount', 'Hemoglobinconcen', 'MeanCorpuscularHemoglobin', 'MeanCorpusVolume']
gwass=['monocytecount' ,'Lymphocytecount3']

for method in methods:
    for gwas in gwass:
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
            df_scores.append(df_score.norm_score)

        df_re=pd.concat(objs=df_scores,axis=1,ignore_index=True)
        df_re.to_csv("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/" + method + "/scDRS_result/bmmcmodel/"+ gwas + "_scDRSresult.csv", sep=",")


df_scores=[]
names=[]
for method in methods:
    for gwas in gwass:
        df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/" + method + "/scDRS_result/" + gwas + ".geneset.gs")
        n=method + '_top' + str(100)
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
        df_scores.append(df_score.norm_score)
        n=n + "_" +gwas
        names.append(n)

df_re=pd.concat(objs=df_scores,axis=1,ignore_index=True)
df_re.columns=names
df_re.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/ThreeMethods_scDRSresult_bmmcmodel.csv", sep=",")

######MAGMA
df_scores=[]
names=[]
for i in range(1,101):
    df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/magma_scRDS_10_1000.csv", index_col=0)
    n= "MAGMAtop" + str(i)
    df_gs = df_gs.loc[[n],:]
    df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/geneset.gs", sep="\t", index=False)
    df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/geneset.gs")
    df_score = scdrs.score_cell(
        data=adata,
        gene_list=df_gs[n][0],
        gene_weight=df_gs[n][1],
        ctrl_match_key="mean_var",
        n_ctrl=20,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False)
    df_scores.append(df_score.norm_score)
    names.append(n)

df_re=pd.concat(objs=df_scores,axis=1,ignore_index=True)
df_re.columns=names
df_re.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/MAGMA_scDRSresult_bmmcmodel.csv", sep=",")


df_scores=[]
names=[]
for i in range(1,101):
    df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Lymphocytecount/magma_scRDS_10_1000.csv", index_col=0)
    n= "MAGMAtop" + str(i)
    df_gs = df_gs.loc[[n],:]
    df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Lymphocytecount/geneset.gs", sep="\t", index=False)
    df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Lymphocytecount/geneset.gs")
    df_score = scdrs.score_cell(
        data=adata,
        gene_list=df_gs[n][0],
        gene_weight=df_gs[n][1],
        ctrl_match_key="mean_var",
        n_ctrl=20,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False)
    df_scores.append(df_score.norm_score)
    names.append(n)

df_re=pd.concat(objs=df_scores,axis=1,ignore_index=True)
df_re.columns=names
df_re.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Lymphocytecount/MAGMA_scDRSresult_bmmcmodel_Lymphocytecount.csv", sep=",")
