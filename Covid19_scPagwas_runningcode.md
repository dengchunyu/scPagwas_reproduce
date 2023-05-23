# Covid19 traits for PBMC scRNA-seq data for scPagwas

### Severe

> covid
> An object of class Seurat 
> 37985 features across 104484 samples within 2 assays 
> Active assay: SCT (17803 features, 0 variable features)
>  1 other assay present: RNA



### moderate

> An object of class Seurat 
> 37985 features across 204508 samples within 2 assays 
> Active assay: SCT (17803 features, 0 variable features)
>  1 other assay present: RNA



### mild

> An object of class Seurat 
> 37985 features across 122686 samples within 2 assays 
> Active assay: SCT (17803 features, 0 variable features)
>  1 other assay present: RNA



### Normal

> An object of class Seurat 
> 37985 features across 37775 samples within 2 assays 
> Active assay: SCT (17803 features, 0 variable features)
>  1 other assay present: RNA

## scPagwas v1.9

The scPagwas and other results for Figure5

### 1. severe

#### Run celltypes scPagwas result

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test5.6")
severe_all<-readRDS("/share2/pub/jiangdp/jiangdp/COVID/data/severe_all.rds")
table(Idents(severe_all))
#Effector CD8+T cell   Memory CD8+T cell    Naive CD8+T cell    Naive CD4+T cell 
#               7245                3127                 766               18525 
# Memory CD4+T cell        CD14+monocyte                  NK        Naive B cell 
#               7034               44577               10922                5230 
#      CD16+monocyte                  DC            Platelet     CD34+Progenitor 
#               3634                1160                1425                 830 
#      Mature B cell 
#                  9

Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/severe_all.rds",
                     singlecell=F,
                      output.prefix="test",
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
Pagwas$Pathway_ld_gwas_data<-NULL
save(Pagwas,file="scpagwas_severe.v1.9.1.RData")

load("scpagwas_severe.v1.9.1.RData")
Pagwas<-scPagwas_main(Pagwas = Pagwas,
                     gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                       assay="RNA",
                     block_annotation = block_annotation,
                     singlecell=T,
                      output.prefix="severe",
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file="scpagwas_severe.v1.9.1.RData")

Bootstrap_P_Barplot(p_results=Pagwas$bootstrap_results$bp_value[-1],
                                p_names=rownames(Pagwas$bootstrap_results)[-1],
                                title = "severe",
                                figurenames = "barplot_severe_kegg.pdf",
                                width = 5,
                                height = 7,
                                do_plot=F)

```

#### Run sub celltypes

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
celltype<-c("Effector CD8+T cell","Memory CD8+T cell","Naive CD8+T cell" ,"Naive CD4+T cell" ,"Memory CD4+T cell","CD14+monocyte","NK","Naive B cell", "CD16+monocyte","DC","Platelet","CD34+Progenitor" )

for(i in celltype){
Pagwas<-scPagwas_main(Pagwas = NULL,
                       gwas_data ="/share/pub/dengcy/GWAS_Multiomics/test/covid19/COVID19_GWAS.txt",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/severe_",i,".rds"),
                      output.prefix="COVID19severe_kegg",
                     output.dirs=i,
                      ncores=5,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file=paste0("scpagwas_severe_",i,".RData"))
}
############cd14monocyte
i<-celltype[6]
for(j in 1:5){
 Pagwas<-scPagwas_main(Pagwas = NULL,
                       gwas_data ="/share/pub/dengcy/GWAS_Multiomics/test/covid19/COVID19_GWAS.txt",
                     add_eqtls="OnlyTSS",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/severe_",i,j,".rds"),
                      output.prefix="COVID19severe_kegg",
                      celltype=F,
                      split_n=2,
                      ncores=10,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file=paste0("scpagwas_severe_",i,j,".RData"))   
}
```

#### Integrate the result：

```R
#severe
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
i<-celltype[6]
Pagwas_l<-lapply(1:5,function(j){
load(paste0("scpagwas_severe_",i,j,".RData"))
 return(Pagwas)   
})

Pagwas_CD14_mono <- merge_pagwas(Pagwas_list = Pagwas_l,n_topgenes = 1000)   
Pagwas_l<-lapply(celltype[-6],function(i){
load(paste0("scpagwas_severe_",i,".RData"))
 return(Pagwas)   
})
Pagwas <- merge_pagwas(Pagwas_list = Pagwas_l,n_topgenes = 1000)   
Pagwas <- merge_pagwas(Pagwas_list = c(Pagwas,Pagwas_CD14_mono),n_topgenes = 1000)  
```

#### Visualization

```R
load("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_severe.v1.91.RData")
Pagwas1 <- FindVariableFeatures(Pagwas1,nfeatures = 5000,assay ="RNA")
#Pagwas1 <- NormalizeData(Pagwas1, normalization.method = "LogNormalize", scale.factor = 10000)
Pagwas1 <- ScaleData(Pagwas1,assay = "RNA")
Pagwas1 <- RunPCA(object = Pagwas1, assay = "RNA", npcs = 50)
Pagwas1 <- RunTSNE(object = Pagwas1,assay =  "RNA", reduction = "pca",dims = 1:50)
Pagwas1 <- RunUMAP(object = Pagwas1, assay = "RNA", reduction = "pca",dims = 1:50)

 scPagwas_Visualization(Single_data=Pagwas1,
                        p_thre = 0.05,
                        FigureType = "umap",
                        width = 5,
                        height =5,
                        lowColor = "white", 
                        highColor = "red",
                        output.dirs="COVID19severe_kegg",
                        size = 0.3,
                        do_plot = F)
save(Pagwas1,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_severe.v1.91.RData")
```

#### MAGMA

```
gwas<-bigreadr::fread2("/share/pub/dengcy/GWAS_Multiomics/test/covid19/COVID19_GWAS.txt")
write.table(gwas,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/COVID19_GWAS.txt",quote=F)

gwas$N <- 
magma_Input1<-gwas[,c("rsid","pvalue","N")]
magma_Input2<-gwas[,c("rsid","chrom","pos","pos")]
write.table(magma_Input2,file=paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/magma_Input2.txt"),sep="\t",row.names=F,quote=F,col.names=F)
colnames(magma_Input1)<-c("SNP","P","N")
write.table(magma_Input1,file=paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/magma_Input1.txt"),sep="\t",row.names=F,quote=F)

cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
./magma --annotate window=10,10 --snp-loc /share/pub/dengcy/GWAS_Multiomics/test/covid19/magma_Input2.txt \
--gene-loc /share/pub/dengcy/Singlecell/COVID19/MAGMA/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded \
--out /share/pub/dengcy/GWAS_Multiomics/test/covid19/annotated_10kbup_10_down

./magma --bfile /share/pub/dengcy/Singlecell/COVID19/MAGMA/g1000_eur/g1000_eur \
--pval /share/pub/dengcy/GWAS_Multiomics/test/covid19/magma_Input1.txt ncol=3 \
--gene-annot /share/pub/dengcy/GWAS_Multiomics/test/covid19/annotated_10kbup_10_down.genes.annot \
--out /share/pub/dengcy/GWAS_Multiomics/test/covid19/annotated_10kbup_10down
```

#### scDRS

1.get pre data

export OPENBLAS_NUM_THREADS=1

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19")
load("scpagwas_severe.v1.91.RData")

library(SeuratDisk)
library(Seurat)
a<-GetAssayData(Pagwas1,assay = "RNA",slot="counts")
Pagwas1 <- CreateSeuratObject(counts = a,
                              project = "scPagwas",
                              min.cells = 3,
                              meta.data =Pagwas1@meta.data,
                              min.features = 200)
DefaultAssay(Pagwas1) <- "RNA"
SaveH5Seurat(Pagwas1, "severe_pagwas.h5seurat")
Convert("severe_pagwas.h5seurat", dest="h5ad")
############GENE FILES
magma_genes<-read.table("/share/pub/dengcy/GWAS_Multiomics/test/covid19/annotated_10kbup_10down.genes.out",header=T)
magma_genes<-magma_genes[order(magma_genes$P,decreasing=F),]
scPagwas_genes<-names(Pagwas1@misc$gene_heritability_correlation[order(Pagwas1@misc$gene_heritability_correlation,decreasing=T),])
i<-"severe"

library(org.Hs.eg.db)
library(dplyr)
library(data.table)
#a<-data.frame(gene_id=magma_genes$GENE)
g2s=toTable(org.Hs.egSYMBOL)
colnames(g2s)<-c("GENE", "symbol")
#g2e=toTable(org.Hs.egENSEMBL)
magma_genes=merge(magma_genes,g2s,by="GENE",all.x=T)
magma_genes<-magma_genes[order(magma_genes$P,decreasing=F),]


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
a<-a[,c(56,1:55)]
write.table( a,file="severe_cov.cov",sep="\t",row.names=F,quote=F,)
```

2.scDRS

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
DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/test/covid19"
H5AD_FILE = os.path.join(DATA_PATH, "severe_pagwas.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
sc.pp.filter_cells(adata, min_genes=0)
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
scdrs.preprocess(adata)

#adata.write('/share/pub/dengcy/GWAS_Multiomics/test/covid19/severe_pagwas.h5ad')
i="severe"
df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scPagwas.magmatopgenes_scRDS"+i+".csv", index_col=0)
df_gs = df_gs.loc[["scPagwastop1000","magmatop1000","scPagwastop500","magmatop500"],:]
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/"+i+".geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/test/covid19/"+i+".geneset.gs")

dict_df_score = dict()
for trait in df_gs:
    gene_list, gene_weights = df_gs[trait]
    dict_df_score[trait] = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=500,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False,
    )
     dict_df_score[trait].to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+trait+".df_res.csv", sep=",", index=False)

df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=dict_df_score["magmatop1000"],
    group_cols=["annotation"],
)["annotation"]
df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+".magmatop1000.scDRS.celltype.csv", sep=",", index=False)

df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=dict_df_score["scPagwastop1000"],
    group_cols=["annotatione"],
)["annotation"]
df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+".scPagwastop1000.scDRS.celltype.csv", sep=",", index=False)



##################
df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scPagwas.magmatopgenes_scRDS"+i+".csv", index_col=0)
df_gs = df_gs.loc[["scPagwastop1000","magmatop1000"],:]
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/"+i+".geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/test/covid19/"+i+".geneset.gs")

gene_list, gene_weights = df_gs["magmatop1000"]
df_score = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=1000,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False)
df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+"magmatop1000.df_res.csv", sep=",")
df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=df_score,
    group_cols=["annotation"],
)["annotation"]
df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+".magmatop1000.scDRS.celltype.csv", sep=",")
```



### 2. mild

#### Run celltypes scPagwas result

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))

 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test5.6")
mild_all<-readRDS("/share2/pub/jiangdp/jiangdp/COVID/data/mild_all.rds")
table(Idents(mild_all))
#Effector CD8+T cell   Memory CD8+T cell    Naive CD8+T cell    Naive CD4+T cell 
#               7245                3127                 766               18525 
# Memory CD4+T cell        CD14+monocyte                  NK        Naive B cell 
#               7034               44577               10922                5230 
#      CD16+monocyte                  DC            Platelet     CD34+Progenitor 
#               3634                1160                1425                 830 
#      Mature B cell 
#                  9
lapply(as.vector(unique(Idents(mild_all))[2:12]),function(x){
  a<-mild_all[,Idents(mild_all)==x] 
  saveRDS(a,file=paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/mild_",x,".rds"))
})

###########
library(scPagwas)
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/mild_all.rds",
                     singlecell=F,
                     output.prefix="mild",
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
Pagwas$Pathway_ld_gwas_data<-NULL
save(Pagwas,file="scpagwas_mild.v1.73.RData")
```

#### Run sub celltypes

```R

library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
celltype<-c("Effector CD8+T cell","Memory CD8+T cell","Naive CD8+T cell" ,"Naive CD4+T cell" ,"Memory CD4+T cell","CD14+monocyte","NK","Naive B cell", "CD16+monocyte","DC","Platelet","CD34+Progenitor" )

Pagwas<-list()
gwas_data <- bigreadr::fread2("/share/pub/dengcy/GWAS_Multiomics/test/covid19/COVID19_GWAS.txt")
Pagwas <- GWAS_summary_input(
    Pagwas = Pagwas,
    gwas_data = gwas_data,
    maf_filter = 0.1
  )
Pagwas$snp_gene_df <- Snp2Gene(snp = Pagwas$gwas_data, refGene = block_annotation, marg = 10000)
names(Pagwas)
i<-"Effector CD8+T cell"
Pagwas_l<-lapply(celltype,function(i){
    Pagwas<-scPagwas_main(Pagwas = Pagwas,
                       assay="RNA",
                          gwas_data =NULL,
                     block_annotation = block_annotation,
                     Single_data = paste0("/share2/pub/jiangdp/jiangdp/COVID/data/mild_",i,".rds"),
                      output.prefix="COVID19mild_kegg",
                     output.dirs=i,
                     ncores=5,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file=paste0("/share2/pub/jiangdp/jiangdp/COVID/data/scpagwas_mild_",i,".RData"))
    return(Pagwas)
})

if(length(SOAR::Objects())>0){
 SOAR::Remove(SOAR::Objects()) 
}
Pagwas_integrate <- merge_pagwas(Pagwas_list = Pagwas_l),n_topgenes = 1000)
save(Pagwas_integrate,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_mild.v1.91.RData")
```

#### Integrate the result

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
setwd("/share2/pub/jiangdp/jiangdp/COVID/data")

celltype <- c("CD14+monocyte","CD16+monocyte","CD34+Progenitor","DC","Effector CD8+T cell" ,"Memory CD4+T cell ","Memory CD8+T cell","Naive B cell","Naive CD4+T cell","Naive CD8+T cell","NK", "Platelet")

Pagwas_l<-lapply(celltype,function(i){
load(paste0("scpagwas_mild_",i,".RData"))
    Pagwas
 return(Pagwas)   
})
Pagwas <- merge_pagwas(Pagwas_list = Pagwas_l,n_topgenes = 1000)
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_mild.v1.91.RData")
```

#### scDRS

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19")
load("scpagwas_mild.v1.10.RData")
library(SeuratDisk)
library(Seurat)
a<-GetAssayData(Pagwas,assay = "RNA",slot="counts")
Pagwas1 <- CreateSeuratObject(counts = a,
                                             project = "scPagwas",
                                             min.cells = 3,
                              meta.data =Pagwas@meta.data,
                                             min.features = 200)
DefaultAssay(Pagwas1) <- "RNA"
SaveH5Seurat(Pagwas1, "mild_pagwas.h5seurat")
Convert("mild_pagwas.h5seurat", dest="h5ad")

############GENE FILES
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19")
load(paste0("scPagwas.magmatopgenes_severe.RData"))
load("scpagwas_mild.v1.10.RData")
scPagwas_genes<-names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation,decreasing=T),])
i<-"mild"

magmatop1000<-topgene$magmatop1000
scPagwastop1000<-scPagwas_genes[1:1000]
magmatop500<-topgene$magmatop500
scPagwastop500<-scPagwas_genes[1:500]
topgene<-list(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500)
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
save(topgene,file=paste0("scPagwas.magmatopgenes_",i,".RData")) 

magmatop1000<-paste(magmatop1000,collapse=",")
scPagwastop1000<-paste(scPagwastop1000,collapse=",")
magmatop500<-paste(magmatop500,collapse=",")
scPagwastop500<-paste(scPagwastop500,collapse=",")

a<-data.frame(genes=c(scPagwastop1000,magmatop1000))
rownames(a)<-c("scPagwastop1000","magmatop1000")
write.csv(a,file=paste0("scPagwas.magmatopgenes_scRDS",i,".csv"))
b<-data.frame(index=colnames(Pagwas),const=rep(1,ncol(Pagwas)))
write.table( b,file=paste0(i,"_cov.cov"),sep="\t",row.names=F,quote=F,)
```

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
DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/test/covid19"
H5AD_FILE = os.path.join(DATA_PATH, "mild_pagwas.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
sc.pp.filter_cells(adata, min_genes=0)
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
scdrs.preprocess(adata)

i="mild"
df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scPagwas.magmatopgenes_scRDS"+i+".csv", index_col=0)
df_gs = df_gs.loc[["scPagwastop1000","magmatop1000"],:]
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/"+i+".geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/test/covid19/"+i+".geneset.gs")

dict_df_score = dict()
for trait in df_gs:
    gene_list, gene_weights = df_gs[trait]
    dict_df_score[trait] = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=500,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False,
    )
     dict_df_score[trait].to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+trait+".df_res.csv", sep=",", index=False)

df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=dict_df_score["magmatop1000"],
    group_cols=["annotation"],
)["annotation"]
df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+".magmatop1000.scDRS.celltype.csv", sep=",", index=False)

df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=dict_df_score["scPagwastop1000"],
    group_cols=["annotation"],
)["annotation"]
df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+".scPagwastop1000.scDRS.celltype.csv", sep=",", index=False)


##########
df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scPagwas.magmatopgenes_scRDS"+i+".csv", index_col=0)
df_gs = df_gs.loc[["scPagwastop1000","magmatop1000"],:]
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/"+i+".geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/test/covid19/"+i+".geneset.gs")

gene_list, gene_weights = df_gs["magmatop1000"]
df_score = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=1000,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False)
df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+"magmatop1000.df_res.csv", sep=",")
df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=df_score,
    group_cols=["annotation"],
)["annotation"]
df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+".magmatop1000.scDRS.celltype.csv", sep=",")
```



### 3. moderate

#### Run celltypes scPagwas result

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/moderate_all.rds",
                      output.prefix="COVID19moderate_kegg",
                     singlecell=F,
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
Pagwas$Pathway_ld_gwas_data<-NULL
save(Pagwas,file="scpagwas_moderate.v1.71.RData")


Bootstrap_P_Barplot(p_results=Pagwas$bootstrap_results$bp_value[-1],
                                p_names=rownames(Pagwas$bootstrap_results)[-1],
                                title = "moderate",
                                figurenames = "barplot_moderate_kegg.pdf",
                                width = 5,
                                height = 7,
                                do_plot=F)

 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test5.6")
moderate_all<-readRDS("/share2/pub/jiangdp/jiangdp/COVID/data/moderate_all.rds")

lapply(as.vector(unique(Idents(moderate_all))),function(x){
  a<-moderate_all[,Idents(moderate_all)==x] 
  saveRDS(a,file=paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/moderate_",x,".rds"))
})
```

#### Run sub celltypes

```R
library(scPagwas)
 suppressMessages(library(Seurat))

celltype<-c("Effector CD8+T cell","Memory CD8+T cell","Naive CD8+T cell" ,"Naive CD4+T cell" ,"Memory CD4+T cell","CD14+monocyte","NK","Naive B cell", "CD16+monocyte","DC","Platelet","CD34+Progenitor" )

Pagwas<-list()
gwas_data <- bigreadr::fread2("/share/pub/dengcy/GWAS_Multiomics/test/covid19/COVID19_GWAS.txt")
Pagwas <- GWAS_summary_input(
    Pagwas = Pagwas,
    gwas_data = gwas_data,
    maf_filter = 0.1
  )
Pagwas$snp_gene_df <- Snp2Gene(snp = Pagwas$gwas_data, refGene = block_annotation, marg = 10000)
rm(gwas_data)
Pagwas_l<-lapply(celltype,function(i){
    Pagwas<-scPagwas_main(Pagwas = Pagwas,
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = paste0("/share2/pub/jiangdp/jiangdp/COVID/data/moderate_",i,".rds"),
                      output.prefix="COVID19moderate_kegg",
                     output.dirs=i,
                     ncores=3,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file=paste0("/share2/pub/jiangdp/jiangdp/COVID/data/scpagwas_moderate_",i,".RData"))
    return(Pagwas)
})

if(length(SOAR::Objects())>0){
 SOAR::Remove(SOAR::Objects()) 
}
Pagwas_integrate <- merge_pagwas(Pagwas_list = Pagwas_l),n_topgenes = 1000)
save(Pagwas_integrate,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_moderate.v1.91.RData")
```

#### scDRS

export OPENBLAS_NUM_THREADS=1

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19")
load("scpagwas_moderate.v1.10.RData")
library(SeuratDisk)
library(Seurat)
a<-GetAssayData(Pagwas,assay = "RNA",slot="counts")
Pagwas1 <- CreateSeuratObject(counts = a,
                                             project = "scPagwas",
                                             min.cells = 3,
                              meta.data =Pagwas@meta.data,
                                             min.features = 200)
DefaultAssay(Pagwas1) <- "RNA"
SaveH5Seurat(Pagwas1, "moderate_pagwas.h5seurat")
Convert("moderate_pagwas.h5seurat", dest="h5ad")
rm(Pagwas1)
rm(a)
gc()
############GENE FILES
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19")
load(paste0("scPagwas.magmatopgenes_severe.RData"))
#load("scpagwas_moderate.v1.10.RData")
scPagwas_genes<-names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation,decreasing=T),])
i<-"moderate"

magmatop1000<-topgene$magmatop1000
scPagwastop1000<-scPagwas_genes[1:1000]
magmatop500<-topgene$magmatop500
scPagwastop500<-scPagwas_genes[1:500]
topgene<-list(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500)
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
save(topgene,file=paste0("scPagwas.magmatopgenes_",i,".RData")) 

magmatop1000<-paste(magmatop1000,collapse=",")
scPagwastop1000<-paste(scPagwastop1000,collapse=",")
magmatop500<-paste(magmatop500,collapse=",")
scPagwastop500<-paste(scPagwastop500,collapse=",")

a<-data.frame(genes=c(scPagwastop1000,magmatop1000))
rownames(a)<-c("scPagwastop1000","magmatop1000")
write.csv(a,file=paste0("scPagwas.magmatopgenes_scRDS",i,".csv"))
b<-data.frame(index=colnames(Pagwas),const=rep(1,ncol(Pagwas)))
write.table( b,file=paste0(i,"_cov.cov"),sep="\t",row.names=F,quote=F,)
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
DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/test/covid19"
H5AD_FILE = os.path.join(DATA_PATH, "moderate_pagwas.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
sc.pp.filter_cells(adata, min_genes=0)
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
scdrs.preprocess(adata)

i="moderate"
df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scPagwas.magmatopgenes_scRDS"+i+".csv", index_col=0)
df_gs = df_gs.loc[["scPagwastop1000","magmatop1000"],:]
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/"+i+".geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/test/covid19/"+i+".geneset.gs")

gene_list, gene_weights = df_gs["magmatop1000"]
df_score = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=1000,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False)
df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+"magmatop1000.df_res.csv", sep=",")

df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=df_score,
    group_cols=["annotation"],
)["annotation"]
df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+".magmatop1000.scDRS.celltype.csv", sep=",")

```



### 4. normal

```R

setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
data(Genes_by_pathway_kegg)
 data(block_annotation)
 data(chrom_ld)
 
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                     add_eqtls="OnlyTSS",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/Normal.rds",
                     singlecell=F,
                      output.prefix="test",
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
Pagwas$Pathway_ld_gwas_data<-NULL
save(Pagwas,file="scpagwas_Normal.v1.71.RData")

###########单纯跑细胞类型
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
data(Genes_by_pathway_kegg)
 data(block_annotation)
 data(chrom_ld)
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")

Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                     add_eqtls="OnlyTSS",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/Normal.rds",
                     singlecell=F,
                      output.prefix="Normal",
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
Pagwas$Pathway_ld_gwas_data<-NULL
save(Pagwas,file="scpagwas_Normal.v1.73.RData")
```

#### scDRS

export OPENBLAS_NUM_THREADS=1

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19")
load("scpagwas_Normal.v1.10.RData")
library(SeuratDisk)
library(Seurat)
a<-GetAssayData(Pagwas1,assay = "RNA",slot="counts")
Pagwas1 <- CreateSeuratObject(counts = a,
                                             project = "scPagwas",
                                             min.cells = 3,
                              meta.data =Pagwas1@meta.data,
                                             min.features = 200)
DefaultAssay(Pagwas1) <- "RNA"
SaveH5Seurat(Pagwas1, "Normal_pagwas.h5seurat")
Convert("Normal_pagwas.h5seurat", dest="h5ad")
rm(a)
rm(Pagwas1)
############GENE FILES
load("scpagwas_Normal.v1.10.RData")
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19")
load(paste0("scPagwas.magmatopgenes_severe.RData"))
#load("scpagwas_Normal.v1.10.RData")
scPagwas_genes<-names(Pagwas@misc$gene_heritability_correlation[order(Pagwas1@misc$gene_heritability_correlation,decreasing=T),])
i<-"Normal"

magmatop1000<-topgene$magmatop1000
scPagwastop1000<-scPagwas_genes[1:1000]
magmatop500<-topgene$magmatop500
scPagwastop500<-scPagwas_genes[1:500]
topgene<-list(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500)
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
save(topgene,file=paste0("scPagwas.magmatopgenes_",i,".RData")) 

magmatop1000<-paste(magmatop1000,collapse=",")
scPagwastop1000<-paste(scPagwastop1000,collapse=",")
magmatop500<-paste(magmatop500,collapse=",")
scPagwastop500<-paste(scPagwastop500,collapse=",")

a<-data.frame(genes=c(scPagwastop1000,magmatop1000))
rownames(a)<-c("scPagwastop1000","magmatop1000")
write.csv(a,file=paste0("scPagwas.magmatopgenes_scRDS",i,".csv"))
b<-data.frame(index=colnames(Pagwas),const=rep(1,ncol(Pagwas)))
write.table( b,file=paste0(i,"_cov.cov"),sep="\t",row.names=F,quote=F,)
```



```R
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
DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/test/covid19"
H5AD_FILE = os.path.join(DATA_PATH, "Normal_pagwas.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
sc.pp.filter_cells(adata, min_genes=0)
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
scdrs.preprocess(adata)

i="Normal"
df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scPagwas.magmatopgenes_scRDS"+i+".csv", index_col=0)
df_gs = df_gs.loc[["scPagwastop1000","magmatop1000"],:]
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/"+i+".geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/test/covid19/"+i+".geneset.gs")

#dict_df_score = dict()
#for trait in df_gs:
gene_list, gene_weights = df_gs["magmatop1000"]
df_score = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=500,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False)
df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+"magmatop1000.df_res.csv", sep=",")

df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=df_score,
    group_cols=["annotation"],
)["annotation"]
df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+".magmatop1000.scDRS.celltype.csv", sep=",")

df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=dict_df_score["scPagwastop1000"],
    group_cols=["annotation"],
)["annotation"]
df_stats.to_csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scDRS/"+i+".scPagwastop1000.scDRS.celltype.csv", sep=",", index=False)
```



## Visualize the result

Figure

```R
library(scPagwas)
suppressMessages(library(Seurat))
suppressMessages(library("dplyr"))
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19")
load("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_severe.v1.91.RData")
Pagwas1 <- FindVariableFeatures(Pagwas1,nfeatures = 5000,assay ="RNA")
Pagwas1 <- NormalizeData(Pagwas1, normalization.method = "LogNormalize", scale.factor = 10000)
Pagwas1 <- ScaleData(Pagwas1,assay = "RNA")
Pagwas1 <- RunPCA(object = Pagwas1, assay = "RNA", npcs = 50)
Pagwas1 <- RunTSNE(object = Pagwas1,assay =  "RNA", reduction = "pca",dims = 1:50)
Pagwas1 <- RunUMAP(object = Pagwas1, assay = "RNA", reduction = "pca",dims = 1:50)
save(Pagwas1，file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_severe.v1.91.RData")
 scPagwas_Visualization(Single_data=Pagwas1,
                                   p_thre = 0.05,
                                   FigureType = "tsne",
                                   width = 7,
                                   height = 7,
                                   lowColor = "#000957", 
                        highColor = "#EBE645",
                        output.dirs="Covid19severe_plot_1.91",
                                   size = 0.5,
                                   do_plot = F)
 scPagwas_Visualization(Single_data=Pagwas1,
                                   p_thre = 0.05,
                                   FigureType = "umap",
                                   width = 7,
                                   height = 7,
                                   lowColor = "#000957", 
                        highColor = "#EBE645",
                        output.dirs="Covid19severe_plot_1.91",
                                   size = 0.5,
                                   do_plot = F)
                                   

library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
load("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_moderate.v1.91.RData")
Pagwas1 <- FindVariableFeatures(Pagwas1,nfeatures = 5000,assay ="RNA")
Pagwas1 <- NormalizeData(Pagwas1, normalization.method = "LogNormalize", scale.factor = 10000)
Pagwas1 <- ScaleData(Pagwas1,assay = "RNA")
Pagwas1 <- RunPCA(object = Pagwas1, assay = "RNA", npcs = 50)
Pagwas1 <- RunTSNE(object = Pagwas1,assay =  "RNA", reduction = "pca",dims = 1:50)
Pagwas1 <- RunUMAP(object = Pagwas1, assay = "RNA", reduction = "pca",dims = 1:50)
save(Pagwas1，file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_moderate.v1.91.RData")
 scPagwas_Visualization(Single_data=Pagwas1,
                                   p_thre = 0.05,
                                   FigureType = "tsne",
                                   width = 7,
                                   height = 7,
                                   lowColor = "#000957", 
                        highColor = "#EBE645",
                        output.dirs="Covid19moderate_plot_1.91",
                                   size = 0.5,
                                   do_plot = F)
 scPagwas_Visualization(Single_data=Pagwas1,
                                   p_thre = 0.05,
                                   FigureType = "umap",
                                   width = 7,
                                   height = 7,
                                   lowColor = "#000957", 
                        highColor = "#EBE645",
                        output.dirs="Covid19moderate_plot_1.91",
                                   size = 0.5,
                                   do_plot = F)



library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
load("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_mild.v1.91.RData")
Pagwas1 <- FindVariableFeatures(Pagwas1,nfeatures = 5000,assay ="RNA")
Pagwas1 <- NormalizeData(Pagwas1, normalization.method = "LogNormalize", scale.factor = 10000)
Pagwas1 <- ScaleData(Pagwas1,assay = "RNA")
Pagwas1 <- RunPCA(object = Pagwas1, assay = "RNA", npcs = 50)
Pagwas1 <- RunTSNE(object = Pagwas1,assay =  "RNA", reduction = "pca",dims = 1:50)
Pagwas1 <- RunUMAP(object = Pagwas1, assay = "RNA", reduction = "pca",dims = 1:50)
save(Pagwas1，file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_mild.v1.91.RData")
 scPagwas_Visualization(Single_data=Pagwas1,
                                   p_thre = 0.05,
                                   FigureType = "tsne",
                                   width = 7,
                                   height = 7,
                                   lowColor = "#000957", 
                        highColor = "#EBE645",
                        output.dirs="Covid19mild_plot_1.91",
                                   size = 0.5,
                                   do_plot = F)
 scPagwas_Visualization(Single_data=Pagwas1,
                                   p_thre = 0.05,
                                   FigureType = "umap",
                                   width = 7,
                                   height = 7,
                                   lowColor = "#000957", 
                        highColor = "#EBE645",
                        output.dirs="Covid19mild_plot_1.91",
                                   size = 0.5,
                                   do_plot = F)


library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
load("/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_Normal.v1.91.RData")
Pagwas1 <- FindVariableFeatures(Pagwas1,nfeatures = 5000,assay ="RNA")
Pagwas1 <- NormalizeData(Pagwas1, normalization.method = "LogNormalize", scale.factor = 10000)
Pagwas1 <- ScaleData(Pagwas1,assay = "RNA")
Pagwas1 <- RunPCA(object = Pagwas1, assay = "RNA", npcs = 50)
Pagwas1 <- RunTSNE(object = Pagwas1,assay =  "RNA", reduction = "pca",dims = 1:50)
Pagwas1 <- RunUMAP(object = Pagwas1, assay = "RNA", reduction = "pca",dims = 1:50)
save(Pagwas1，file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_Normal.v1.91.RData")

 scPagwas_Visualization(Single_data=Pagwas1,
                                   p_thre = 0.05,
                                   FigureType = "tsne",
                                   width = 7,
                                   height = 7,
                                   lowColor = "#000957", 
                        highColor = "#EBE645",
                        output.dirs="Covid19Normal_plot_1.91",
                                   size = 0.5,
                                   do_plot = F)
```



