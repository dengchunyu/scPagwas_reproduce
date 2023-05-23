# scDRS for BMMC and PBMC data



## BMMC

**Supplementary Figure S8**

### preprogress data

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
load("Lymphocytecount3_Hema_bmmc_scPagwas_v1.10.0.RData")
library(SeuratDisk)
library(Seurat)
a<-GetAssayData(Pagwas,assay = "RNA",slot="counts")
Pagwas1 <- CreateSeuratObject(counts = a,
                                             project = "scPagwas",
                                             min.cells = 3,
                              meta.data =Pagwas@meta.data,
                                             min.features = 200)
DefaultAssay(Pagwas1) <- "RNA"
SaveH5Seurat(Pagwas1, "bmmc_pagwas.h5seurat")
Convert("bmmc_pagwas.h5seurat", dest="h5ad")
rm(a)
rm(Pagwas1)


traits<-c("monocytecount","MeanCorpusVolume")

for(i in traits){
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
 
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_",i,".RData")) 
magmatop1000<-topgene$magmatop1000
scPagwastop1000<-topgene$scPagwastop1000
magmatop500<-topgene$magmatop500
scPagwastop500<-topgene$scPagwastop500

magmatop1000<-paste(magmatop1000,collapse=",")
scPagwastop1000<-paste(scPagwastop1000,collapse=",")
magmatop500<-paste(magmatop500,collapse=",")
scPagwastop500<-paste(scPagwastop500,collapse=",")

a<-data.frame(genes=c(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500))
rownames(a)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
write.csv(a,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_scRDS",i,".csv"))
}
```

### scDRS

```
export OPENBLAS_NUM_THREADS=1
```

```
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
DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test"
H5AD_FILE = os.path.join(DATA_PATH, "bmmc_pagwas.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
sc.pp.filter_cells(adata, min_genes=0)
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
scdrs.preprocess(adata)
#adata.write('/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/bmmc_pagwas.h5ad')

traits = ['Lymphocytecount3','monocytecount','MeanCorpusVolume']
for i in traits:
 df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_scRDS"+i+".csv", index_col=0)
 df_gs = df_gs.loc[["scPagwastop1000","magmatop1000"],:]
 df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".geneset.gs", sep="\t", index=False)
 df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".geneset.gs")

 df_score = scdrs.score_cell(
        data=adata,
        gene_list=df_gs['scPagwastop1000'][0],
        gene_weight=df_gs['scPagwastop1000'][1],
        ctrl_match_key="mean_var",
        n_ctrl=200,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False,
    )
 df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+"scPagwastop1000.df_res.csv", sep=",")
 df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=df_score,
    group_cols=["celltypes"],
 )["celltypes"]
 df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".scPagwastop1000.scDRS.celltype.csv", sep=",")

 df_score = scdrs.score_cell(
        data=adata,
        gene_list=df_gs['magmatop1000'][0],
        gene_weight=df_gs['magmatop1000'][1],
        ctrl_match_key="mean_var",
        n_ctrl=200,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False,
    )
 df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+"magmatop1000.df_res.csv", sep=",")
 df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=df_score,
    group_cols=["celltypes"],
 )["celltypes"]
df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".magmatop1000.scDRS.celltype.csv", sep=",")
```

### visualize

```
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 library(ggpubr)
 library(ggtext)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
load("Lymphocytecount3_Hema_bmmc_scPagwas_v1.10.0.RData")
fortify.Seurat.umap <- function(x) {
  xy1 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "umap")
  )
  colnames(xy1) <- c("UMAP_1", "UMAP_2")
  xy1$UMAP_1 <- as.numeric(xy1$UMAP_1)
  xy1$UMAP_2 <- as.numeric(xy1$UMAP_2)
  
  return(cbind(xy1, as.data.frame(x@meta.data)))
}


#' fortify.Seurat.tsne
#' @description set data frame to ggplot
#'
#' @param x seruat
#' @export

fortify.Seurat.tsne <- function(x) {
  xy2 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "tsne")
  )
  colnames(xy2) <- c("TSNE_1", "TSNE_2")
  xy2$TSNE_1 <- as.numeric(xy2$TSNE_1)
  xy2$TSNE_2 <- as.numeric(xy2$TSNE_2)
  
  return(cbind(xy2, as.data.frame(x@meta.data)))
}
#############
#read the result for scDRS
a<-read.csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/Lymphocytecount3scPagwastop1000.df_res.csv")

Pagwas$Lymphocytecount3_scPagwas1000_scdrs.score<-a$norm_score
all_fortify_can <- fortify.Seurat.tsne(Pagwas)
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = TSNE_1, y = TSNE_2,color =Lymphocytecount3_scPagwas1000_scdrs.score), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas1000_scdrs")

pdf(file="bmmc_Lymphocytecount3_scPagwas1000_scdrs.score.pdf",width =10, height = 10)
print(p1)
dev.off()

rm(a)
a<-read.csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/Lymphocytecount3magmatop1000.df_res.csv")
Pagwas$Lymphocytecount3_magma1000_scdrs.score<-a$norm_score

all_fortify_can <- fortify.Seurat.tsne(Pagwas)
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = TSNE_1, y = TSNE_2,color =Lymphocytecount3_magma1000_scdrs.score), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas1000_scdrs")

pdf(file="bmmc_Lymphocytecount3_magma1000_scdrs.score.pdf",width =10, height =10)
print(p1)
dev.off()


rm(a)

a<-read.csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/monocytecountscPagwastop1000.df_res.csv")
Pagwas$mono.scPagwas1000_scdrs.score <-a$norm_score
all_fortify_can <- fortify.Seurat.tsne(Pagwas)
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = TSNE_1, y = TSNE_2,color =mono.scPagwas1000_scdrs.score), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas1000_scdrs")

pdf(file="bmmc_monocytecount_scPagwas1000_scdrs.score.pdf",width =10, height =10)
print(p1)
dev.off()


a<-read.csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/monocytecountmagmatop1000.df_res.csv")
Pagwas$mono.magma1000_scdrs.score <-a$norm_score
all_fortify_can <- fortify.Seurat.tsne(Pagwas)
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = TSNE_1, y = TSNE_2,color =mono.magma1000_scdrs.score), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("mono.magma1000_scdrs.score")

pdf(file="bmmc_monocytecount_magma1000_scdrs.score.pdf",width =10, height =10)
print(p1)
dev.off()

a<-read.csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/MeanCorpusVolumescPagwastop1000.df_res.csv")
Pagwas$Corpus.scPagwas1000_scdrs.score <-a$norm_score
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = TSNE_1, y = TSNE_2,color =Corpus.scPagwas1000_scdrs.score), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("Corpus.scPagwas1000_scdrs.score")

pdf(file="bmmc_Corpus.scPagwas1000_scdrs.score.pdf",width =10, height =10)
print(p1)
dev.off()

a<-read.csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/MeanCorpusVolumemagmatop1000.df_res.csv")
Pagwas$Corpus.magma1000_scdrs.score <-a$norm_score
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = TSNE_1, y = TSNE_2,color =Corpus.magma1000_scdrs.score), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("Corpus.magma1000_scdrs.score")

pdf(file="bmmc_Corpus.magma1000_scdrs.score.pdf",width =10, height =10)
print(p1)
dev.off()

###############
load("Lymphocytecount3_Hema_bmmc_scPagwas_v1.10.0.RData")
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = TSNE_1, y = TSNE_2,color =scPagwas.TRS.Score1), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas + Seurat")

pdf(file="bmmc_Lymphocytecount3_scPagwas.score.pdf",width =9, height =9)
print(p1)
dev.off()


load("monocytecount_bmmc_Pagwas1.10.0.RData")
all_fortify_can<-fortify.Seurat.tsne(Pagwas)
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = TSNE_1, y = TSNE_2,color =scPagwas.TRS.Score1), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas + Seurat")

pdf(file="bmmc_monocytecount_scPagwas.score.pdf",width =9, height =9)
print(p1)
dev.off()

load("MeanCorpusVolume_Hema_bmmc_scPagwas_v1.9.1.RData")
all_fortify_can<-fortify.Seurat.tsne(Pagwas)
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = TSNE_1, y = TSNE_2,color =scPagwas.topgenes.Score1), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas + Seurat")

pdf(file="bmmc_MeanCorpusVolume_scPagwas.score.pdf",width =9, height =9)
print(p1)
dev.off()
```



## PBMC

**Supplementary Figure S9**

### Preprogress data

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc")
load("monocytecount_pbmc_scPagwas_singlecell.RData")
library(SeuratDisk)
library(Seurat)
a<-GetAssayData(Pagwas,assay = "RNA",slot="counts")
Pagwas1 <- CreateSeuratObject(counts = a,
                                             project = "scPagwas",
                                             min.cells = 3,
                              meta.data =Pagwas@meta.data,
                                             min.features = 200)
DefaultAssay(Pagwas1) <- "RNA"
SaveH5Seurat(Pagwas1, "pbmc_pagwas.h5seurat")
Convert("pbmc_pagwas.h5seurat", dest="h5ad")
rm(a)
rm(Pagwas1)
head(Pagwas@meta.data)

traits<-c("Lymphocytecount3","monocytecount","MeanCorpusVolume")

for(i in traits[2:3]){
print(i)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/",i,"_pbmc_scPagwas_singlecell.RData"))
scPagwas_genes<-names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation,decreasing=T),])
rm(Pagwas)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_",i,".RData")) 
scPagwastop1000<-scPagwas_genes[1:1000]
scPagwastop500<-scPagwas_genes[1:500]

scPagwastop1000<-paste(scPagwastop1000,collapse=",")
scPagwastop500<-paste(scPagwastop500,collapse=",")
magmatop1000<-paste(topgene$magmatop1000,collapse=",")
magmatop500<-paste(topgene$magmatop500,collapse=",")
a<-data.frame(genes=c(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500))
rownames(a)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
  
write.csv(a,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/PBMCscPagwas.topgenes_scRDS",i,".csv"))
}

```

### scDRS

```python
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
DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc"
H5AD_FILE = os.path.join(DATA_PATH, "pbmc_pagwas.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
sc.pp.filter_cells(adata, min_genes=0)
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
scdrs.preprocess(adata)

traits = ['Lymphocytecount3','monocytecount','MeanCorpusVolume']
for i in traits:
 df_gs = pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/PBMCscPagwas.topgenes_scRDS"+i+".csv", index_col=0)
 df_gs = df_gs.loc[["scPagwastop1000","magmatop1000"],:]
 df_gs = df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/"+i+".geneset.gs", sep="\t", index=False)
 df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/"+i+".geneset.gs")

 df_score = scdrs.score_cell(
        data=adata,
        gene_list=df_gs['scPagwastop1000'][0],
        gene_weight=df_gs['scPagwastop1000'][1],
        ctrl_match_key="mean_var",
        n_ctrl=200,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False,
    )
 df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/"+i+"scPagwastop1000.df_res.csv", sep=",")
 df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=dict_df_score["scPagwastop1000"],
    group_cols=["initial_clustering"],
 )["initial_clustering"]
 df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/"+i+".scPagwastop1000.scDRS.celltype.csv", sep=",")

 df_score = scdrs.score_cell(
        data=adata,
        gene_list=df_gs['magmatop1000'][0],
        gene_weight=df_gs['magmatop1000'][1],
        ctrl_match_key="mean_var",
        n_ctrl=200,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False,
    )
 df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/"+i+"magmatop1000.df_res.csv", sep=",")
 df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=dict_df_score["scPagwastop1000"],
    group_cols=["initial_clustering"],
 )["initial_clustering"]
 df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/"+i+".magmatop1000.scDRS.celltype.csv", sep=",")
    
#!/usr/bin/sh
#PBS -N pbmc.scdrs
#PBS -q workq
#PBS -l nodes=node03
#PBS -l mem=150GB
#PBS -l ncpus=5
#PBS -j oe

export OPENBLAS_NUM_THREADS=1
python /share/pub/dengcy/GWAS_Multiomics/test/covid19/6.py

```

### visualize

```R
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 library(ggpubr)
 library(ggtext)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc")
load("monocytecount_pbmc_scPagwas_singlecell.RData")
fortify.Seurat.umap <- function(x) {
  xy1 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "umap")
  )
  colnames(xy1) <- c("UMAP_1", "UMAP_2")
  xy1$UMAP_1 <- as.numeric(xy1$UMAP_1)
  xy1$UMAP_2 <- as.numeric(xy1$UMAP_2)
  
  return(cbind(xy1, as.data.frame(x@meta.data)))
}


fortify.Seurat.tsne <- function(x) {
  xy2 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "tsne")
  )
  colnames(xy2) <- c("TSNE_1", "TSNE_2")
  xy2$TSNE_1 <- as.numeric(xy2$TSNE_1)
  xy2$TSNE_2 <- as.numeric(xy2$TSNE_2)
  
  return(cbind(xy2, as.data.frame(x@meta.data)))
}
#############
#read the result for scDRS
all_fortify_can <- fortify.Seurat.umap(Pagwas)
rm(Pagwas)
a<-read.csv("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/Lymphocytecount3scPagwastop1000.df_res.csv")
all_fortify_can$Lymph.scPagwas1000_scdrs.score<-a$norm_score
rm(a)

a<-read.csv("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/Lymphocytecount3magmatop1000.df_res.csv")
all_fortify_can$Lymph.magma1000_scdrs.score<-a$norm_score

a<-read.csv("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/monocytecountscPagwastop1000.df_res.csv")
all_fortify_can$mono.scPagwas1000_scdrs.score <-a$norm_score

a<-read.csv("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/monocytecountmagmatop1000.df_res.csv")
all_fortify_can$mono.magma1000_scdrs.score <-a$norm_score
a<-read.csv("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/MeanCorpusVolumescPagwastop1000.df_res.csv")
all_fortify_can$Corpus.scPagwas1000_scdrs.score <-a$norm_score
a<-read.csv("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scDRS/MeanCorpusVolumemagmatop1000.df_res.csv")
all_fortify_can$Corpus.magma1000_scdrs.score <-a$norm_score

p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =scPagwas.topgenes.Score1), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas score")

pdf(file="pbmc_monocyte_scPagwas_score.pdf",width =10, height = 10)
print(p1)
dev.off()


p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =Lymph.scPagwas1000_scdrs.score), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas1000_scdrs")

pdf(file="pbmc_Lymphocytecount3_scPagwas1000_scdrs.zscore.pdf",width =10, height = 10)
print(p1)
dev.off()


p2<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =Lymph.magma1000_scdrs.score), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("magma1000_scdrs.zscore")

pdf(file="pbmc_Lymphocytecount3_magma1000_scdrs.zscore.pdf",width =10, height =10)
print(p2)
dev.off()


p3<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x =UMAP_1, y = UMAP_2,color =mono.scPagwas1000_scdrs.score), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas1000_scdrs")

pdf(file="pbmc_monocytecount_scPagwas1000_scdrs.zscore.pdf",width =10, height =10)
print(p3)
dev.off()


p4<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color = mono.magma1000_scdrs.score), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("magma1000_scdrs.zscore")

pdf(file="pbmc_monocytecount_magma1000_scdrs.zscore.pdf",width = 10, height =10)
print(p4)
dev.off()


p4<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color = Corpus.scPagwas1000_scdrs.score), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas1000_scdrs")

pdf(file="pbmc_Corpus_scPagwas1000_scdrs.zscore.pdf",width = 10, height = 10)
print(p4)
dev.off()

p4<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color = Corpus.magma1000_scdrs.score), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("magma1000_scdrs")

pdf(file="pbmc_Corpus_magma1000_scdrs.zscore.pdf",width = 10, height = 10)
print(p4)
dev.off()

save(all_fortify_can,file="all_fortify_can.RData")


load("Lymphocytecount3_pbmc_scPagwas_singlecell.RData")
all_fortify_can <- fortify.Seurat.umap(Pagwas)
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =scPagwas.topgenes.Score1), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas score")

pdf(file="pbmc_Lymphocytecount3_scPagwas_score.pdf",width =10, height = 10)
print(p1)
dev.off()
```

