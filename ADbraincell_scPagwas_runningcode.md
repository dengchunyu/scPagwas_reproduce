# AD gwas data based on bain single cell

## 1. Single cell data progress

```R
library(Seurat)
library(SeuratData)
library(dplyr)
Single_data = readRDS("/share/pub/qiuf/brain/01-data/AD/GSE138852.rds")
#Single_data =readRDS("E:/RPakage/scPagwas/inst/extdata/GSE138852_ad.rds")
Single_data <- FindVariableFeatures(Single_data,nfeatures = 3000)
Single_data <- NormalizeData(Single_data, normalization.method = "LogNormalize", scale.factor = 10000)
Single_data <- ScaleData(Single_data)

Single_data <- RunPCA(object = Single_data, assay = "RNA", npcs = 50)
Single_data <- RunTSNE(object = Single_data,assay =  "RNA", reduction = "pca",dims = 1:50)
Single_data <- RunUMAP(object = Single_data, assay = "RNA", reduction = "pca",dims = 1:50)
#Idents(Single_data)
saveRDS(Single_data,file = "/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE138852.rds")

Single_data<-Single_data[,Single_data@meta.data$disease.state=="AD"]
saveRDS(Single_data,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE138852_ad.rds")
###########################################
library(Seurat)
library(SeuratData)
library(dplyr)
Single_data = readRDS("/share/pub/qiuf/brain/01-data/AD/supplement_resolution0.2.rds")
Single_data <- FindVariableFeatures(Single_data,nfeatures = 3000)
Single_data <- NormalizeData(Single_data, normalization.method = "LogNormalize", scale.factor = 10000)
Single_data <- ScaleData(Single_data)

Single_data <- RunPCA(object = Single_data, assay = "RNA", npcs = 50)
Single_data <- RunTSNE(object = Single_data,assay =  "RNA", reduction = "pca",dims = 1:50)
Single_data <- RunUMAP(object = Single_data, assay = "RNA", reduction = "pca",dims = 1:50)
save(Single_data,file = "/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE160936.rds")

```

## 2. GSE138852

#### v1.10.0

export OPENBLAS_NUM_THREADS=1

```R
 library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))

setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test")

Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE138852.rds",
                     output.prefix="AD",
                     output.dirs="GSE138852AD_outputv1.10",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     seruat_return=T,
                     celltype=T,
                     ncores = 10)
save(Pagwas,file="Pagwas_seu_GSE138852_v1.10.RData")
load("Pagwas_seu_GSE138852_v1.10.RData")
#删除CellScalepValue，CellScaleqValue
Pagwas@meta.data[,c("CellScalepValue","CellScaleqValue")]<-NULL
 scPagwas_topgenes <- names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation, decreasing = T), ])[1:500]
 correct_pdf<-Get_CorrectBg_p(Single_data=Pagwas,
                                 scPagwas.TRS.Score=Pagwas$scPagwas.TRS.Score1,
                                 iters_singlecell=200,
                                 n_topgenes=500,
                                 scPagwas_topgenes=scPagwas_topgenes
    )
Pagwas$Random_Correct_BG_p <- correct_pdf$pooled_p
Pagwas$Random_Correct_BG_adjp <- correct_pdf$adj_p
Pagwas$Random_Correct_BG_z <- correct_pdf$pooled_z
scPagwas_Visualization(Single_data =Pagwas,
                                   p_thre = 0.05,
                                   output.dirs = "GSE138852AD_outputv1.10",
                                   FigureType = "umap",
                                   width = 7,
                                   height = 7,
                                   title = "",
                                   lowColor = "#000957",
                                   highColor = "#EBE645",
                                   size = 0.5,
                                   do_plot = TRUE)
saveRDS(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/Pagwas_seu_GSE138852_v1.10.rds")
```

## 3. GSE160936

```R
 library(scPagwas)
 suppressMessages(library(Seurat))
  library(parallel)
 suppressMessages(library("dplyr"))
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.16test")
 Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     output.prefix="test",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE160936.rds",
                     ncores=1,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)


```

## 4.magma

```
gwas<-bigreadr::fread2("/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt")

cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
./magma --annotate window=10,10 --snp-loc /share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/magma_input2.txt \
--gene-loc /share/pub/dengcy/Singlecell/COVID19/MAGMA/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded \
--out /share/pub/dengcy/GWAS_Multiomics/ad_test/annotated_10kbup_10_down

./magma --bfile /share/pub/dengcy/Singlecell/COVID19/MAGMA/g1000_eur/g1000_eur \
--pval /share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/magma_input1.txt ncol=3 \
--gene-annot /share/pub/dengcy/GWAS_Multiomics/test/covid19/annotated_10kbup_10_down.genes.annot \
--out /share/pub/dengcy/GWAS_Multiomics/ad_test/annotated_10kbup_10down
```

## 5.scDRS

```
export OPENBLAS_NUM_THREADS=1
```

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test")
a<-readRDS("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE138852.rds")
library(SeuratDisk)
library(Seurat)
DefaultAssay(a) <- "RNA"
SaveH5Seurat(a, "GSE138852_pagwas.h5seurat")
Convert("GSE138852_pagwas.h5seurat", dest="h5ad")

load("/share/pub/dengcy/GWAS_Multiomics/ad_test/Pagwas_seu_GSE138852_v1.10.RData")
magma_genes<-read.table("/share/pub/dengcy/GWAS_Multiomics/ad_test/annotated_10kbup_10down.genes.out",header=T)
magma_genes<-magma_genes[order(magma_genes$P,decreasing=F),]
scPagwas_genes<-names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation,decreasing=T),])



library(org.Hs.eg.db)
library(dplyr)
library(data.table)
#a<-data.frame(gene_id=magma_genes$GENE)
g2s=toTable(org.Hs.egSYMBOL)
colnames(g2s)<-c("GENE", "symbol")
#g2e=toTable(org.Hs.egENSEMBL)
magma_genes=merge(magma_genes,g2s,by="GENE",all.x=T)
magma_genes<-magma_genes[order(magma_genes$P,decreasing=F),]
i<-"AD"
magmatop1000<-magma_genes$symbol
magmatop1000<-magmatop1000[1:1000]
scPagwastop1000<-scPagwas_genes[1:1000]
magmatop500<-magmatop1000[1:500]
scPagwastop500<-scPagwas_genes[1:500]
topgene<-list(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500)
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
save(topgene,file=paste0("scPagwas.magmatopgenes_",i,".RData")) 

magmatop1000<-paste(magmatop1000,collapse=",")
scPagwastop1000<-paste(scPagwastop1000,collapse=",")
magmatop500<-paste(magmatop500,collapse=",")
scPagwastop500<-paste(scPagwastop500,collapse=",")

a<-data.frame(genes=c(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500))
rownames(a)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
write.csv(a,file=paste0("scPagwas.magmatopgenes_scRDS",i,".csv"))

b<-data.frame(index=colnames(Pagwas),const=rep(1,ncol(Pagwas)))
write.table( b,file="AD_cov.cov",sep="\t",row.names=F,quote=F,)


#####python
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
warnings.filterwarnings("ignore")
# load adata

DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/ad_test"
H5AD_FILE = os.path.join(DATA_PATH, "GSE138852_pagwas.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)

sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
scdrs.preprocess(adata)
####################
i="AD"
df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/ad_test/scPagwas.magmatopgenes_scRDS"+i+".csv", index_col=0)
df_gs = df_gs.loc[["scPagwastop1000","magmatop1000","scPagwastop500","magmatop500"],:]
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/ad_test/"+i+".geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/ad_test/"+i+".geneset.gs")

dict_df_score = dict()
for trait in df_gs:
    gene_list, gene_weights = df_gs[trait]
    dict_df_score[trait] = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=200,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False,
    )
    dict_df_score[trait].to_csv("/share/pub/dengcy/GWAS_Multiomics/ad_test/"+i+trait+".df_res.csv", sep=",", index=False)

df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=dict_df_score["magmatop1000"],
    group_cols=["oupSample.cellType"],
)["oupSample.cellType"]
df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/ad_test/"+i+".magmatop1000.scDRS.celltype.csv", sep=",", index=False)

df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=dict_df_score["scPagwastop1000"],
    group_cols=["oupSample.cellType"],
)["oupSample.cellType"]
df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/ad_test/"+i+".scPagwastop1000.scDRS.celltype.csv", sep=",", index=False)

```





