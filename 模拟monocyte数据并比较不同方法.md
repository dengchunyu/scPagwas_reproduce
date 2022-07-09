# 模拟monocyte数据

数据来自：

Interrogation of human hematopoiesis at single-cell and single-variant resolution

```R
library(scDesign2)
library(copula)    # corKendall
library(Rtsne)
library(plyr)      # mapvalues
library(reshape2)  # melt
library(gridExtra) # arrangeGrob
library(ggpubr)    # as_ggplot
library(cowplot)   # draw_plot_label
library(ggplot2); theme_set(theme_bw());
set.seed(42)
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")
bulk_count<-read.delim("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/16populations_RNAcounts.txt",header=T)
bulk_count<-bulk_count[!duplicated(bulk_count$Genes),]
rownames(bulk_count)<-bulk_count$Genes
bulk_count<-bulk_count[,-1]
bulk_count<-data.matrix(bulk_count)
bulk_count2<-bulk_count[,c(10,10,15,15)]

copula_result <- fit_model_scDesign2(bulk_count2, cell_type_sel=c('Mono','NK'), sim_method = 'copula',ncores = 10)
sim_count_1000 <- simulate_count_scDesign2(copula_result, n_cell_new=1000, sim_method = 'copula',cell_type_prop = c(0.5,0.5))
rownames(sim_count_1000)<-rownames(bulk_count2)
saveRDS(copula_result, file = 'copula_result_population.rds')
saveRDS(sim_count_1000, file = 'sim_count_1000_population.rds')

#sim_count_1000<-readRDS("sim_count_1000.rds")
#构造seruat格式
library(SeuratObject)
library(Seurat)
sim_count_1000<-readRDS("sim_count_1000.rds")
sim_data<-CreateSeuratObject(counts=sim_count_1000,assay = "RNA")
sim_data$celltype<-colnames(sim_count_1000)
sim_data$type<- c(rep(1,500),rep(0,500))
sim_data <- NormalizeData(sim_data, normalization.method = "LogNormalize", scale.factor = 10000)
sim_data <- ScaleData(sim_data)

sim_data<-FindVariableFeatures(object=sim_data)

Idents(sim_data)<-sim_data$celltype
#sim_data<-RunTsne(object=sim_data)
saveRDS(sim_data,file="sim_data.rds")


```

### 利用单细胞数据构造模拟数据

```R
library(scDesign2)
library(copula)    # corKendall
library(Rtsne)
library(plyr)      # mapvalues
library(reshape2)  # melt
library(gridExtra) # arrangeGrob
library(ggpubr)    # as_ggplot
library(cowplot)   # draw_plot_label
library(ggplot2); theme_set(theme_bw())
library(SeuratObject)
library(Seurat)
set.seed(42)

setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")
load("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/phe_NM_Healthy_pbmc.Rdata")
Idents(NM_Healthy_pbmc)<-NM_Healthy_pbmc@meta.data$full_clustering
#GetAssayData(object=NM_Healthy_pbmc[["RNA"]],slot="count")[1:5,1:5]

NM_Healthy_pbmc<-subset(NM_Healthy_pbmc,idents=c("CD14_mono","CD16_mono","NK_16hi"))

N1_pbmc<-subset(NM_Healthy_pbmc,idents=c("CD14_mono","CD16_mono"))
Idents(N1_pbmc)<-N1_pbmc@meta.data$sample_id
N1_pbmc_meanexpr <- AverageExpression(N1_pbmc, group.by = "ident",assays = "RNA",slot = "data")
Mono_pbmc_meanexpr<-N1_pbmc_meanexpr[[1]]
colnames(Mono_pbmc_meanexpr)<-paste0("Mono",1:ncol(Mono_pbmc_meanexpr))
Mono_pbmc_meanexpr<-Mono_pbmc_meanexpr*10000

N1_pbmc<-subset(NM_Healthy_pbmc,idents=c("NK_16hi"))
Idents(N1_pbmc)<-N1_pbmc@meta.data$sample_id
N1_pbmc_meanexpr <- AverageExpression(N1_pbmc, group.by = "ident",assays = "RNA",slot = "data")
NK_pbmc_meanexpr<-N1_pbmc_meanexpr[[1]]
colnames(NK_pbmc_meanexpr)<-paste0("NK",1:ncol(NK_pbmc_meanexpr))
NK_pbmc_meanexpr<-NK_pbmc_meanexpr*10000


bulk_count<-cbind(Mono_pbmc_meanexpr,NK_pbmc_meanexpr)
#bulk_count<-apply(bulk_count,2,as.integer)
rownames(bulk_count)<-rownames(NM_Healthy_pbmc)
colnames(bulk_count)<-c(rep("Mono",23),rep("NK",24))

copula_result <- fit_model_scDesign2(bulk_count, cell_type_sel=c('Mono','NK'), sim_method = 'copula',ncores = 10)

sim_count_1000 <- simulate_count_scDesign2(copula_result, n_cell_new=1000, sim_method = 'copula',cell_type_prop = c(0.5,0.5))

saveRDS(copula_result, file = 'copula_result2.rds')
saveRDS(sim_count_1000, file = 'sim_count2_1000.rds')
rownames(sim_count_1000)<-rownames(bulk_count)
#构造seruat格式


#sim_count_1000<-readRDS("sim_count_1000.rds")
sim_data<-CreateSeuratObject(counts=sim_count_1000,assay = "RNA")
sim_data$celltype<-colnames(sim_count_1000)
sim_data$type<- c(rep(1,500),rep(0,500))
sim_data <- NormalizeData(sim_data, normalization.method = "LogNormalize", scale.factor = 10000)
sim_data <- ScaleData(sim_data)
Idents(sim_data)<-sim_data$celltype
sim_data <- FindVariableFeatures(sim_data,nfeatures = 3000)
sim_data <- RunPCA(object = sim_data, assay = "RNA", npcs = 50)
sim_data <- RunTSNE(object = sim_data,assay =  "RNA", reduction = "pca",dims = 1:50)
sim_data <- RunUMAP(object = sim_data, assay = "RNA", reduction = "pca",dims = 1:50)

#sim_data<-RunTsne(object=sim_data)
saveRDS(sim_data,file="sim_data_NK.rds")
############绘制模拟数据的tsne
pdf(file="modeldata_UMAP.pdf",height=5,width=5)
DimPlot(sim_data,group.by="celltype",reduction="umap",pt.size=0.3,label = TRUE, repel=TRUE,label.size = 4)+ 
umap_theme()+ labs(x="UMAP",y="")+
ggtitle("Model data")+
        scale_colour_manual(name = "celltype", values = c("#ff7f0e","#1f77b4")) +
        theme(aspect.ratio=1)
dev.off()

```

### 计算scPagwas

```
 library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))

setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")
i<-"monocytecount"

Pagwas<-scPagwas_main(Pagwas = NULL,
gwas_data=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"),
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data.rds",
                     output.prefix="modeldata_dc",
                     output.dirs="modeldata_outputv1.9.1",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     seruat_return=T,
                     celltype=F,
                     ncores = 15,
                     split_n=1)
                     
save(Pagwas,file="Pagwas_monocytecount_nk_v1.9.1.RData")
#load("Pagwas_monocytecount_nk_v1.9.1.RData")


 scPagwas_Visualization(Single_data=Pagwas,
                        p_thre = 0.05,
                        FigureType = "umap",
                        width = 5,
                        height =5,
                        lowColor = "white", 
                        highColor = "red",
                        output.dirs="modeldata_outputv1.9.1",
                        size = 0.3,
                        do_plot = F)
  
   library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))

setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")

Pagwas<-scPagwas_main(Pagwas = NULL,
gwas_data="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/Lymphocytecount3_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_NK.rds",
                     output.prefix="modeldata_NK",
                     output.dirs="modeldata_Lymphocytecount_outputv1.9.1",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     seruat_return=T,
                     celltype=F,
                     ncores = 15,
                     split_n=1)
save(Pagwas,file="Pagwas_Lymphocytecount_nk_v1.9.1.RData")

 scPagwas_Visualization(Single_data=Pagwas,
                        p_thre = 0.05,
                        FigureType = "umap",
                        width = 5,
                        height =5,
                        lowColor = "white", 
                        highColor = "red",
                        output.dirs="modeldata_Lymphocytecount_outputv1.9.1",
                        size = 0.3,
                        do_plot = F)


```

#!/usr/bin/sh
#PBS -N roly1
#PBS -q workq
#PBS -l nodes=node05
#PBS -l ncpus=5
#PBS -j oe
####mild

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scPagwas1.8/r1.r



### rolypoly

```
lapply(list.files("/share/pub/dengcy/GWAS_Multiomics/Pkg/rolypoly-master/rolypoly-master/R/"),function(x){
source(paste0("/share/pub/dengcy/GWAS_Multiomics/Pkg/rolypoly-master/rolypoly-master/R/",x))
})

library(Seurat)
library(SeuratObject)
library("dplyr")
library("foreach")
library("ggplot2")
library("glmnet")
library("MASS")
library("Matrix")
library("matrixcalc")

load("/share/pub/dengcy/Singlecell/COVID19/data/covid_ld.RData")
library("data.table")
load("/share/pub/dengcy/Singlecell/COVID19/1.rolypoly_result/geneid_df1.RData")
i<-"monocytecount"
suppressMessages(gwas_data <- bigreadr::fread2(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt")))
 Single_data =readRDS("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_NK.rds")
 Single_mat<-GetAssayData(Single_data,assay="RNA")
 
 geneid_df1<-geneid_df1[!duplicated(geneid_df1$label),]
ld_path <- "/share/pub/dengcy/Singlecell/COVID19/data/LD"
rolypoly_result <- rolypoly_roll(
  gwas_data =as.data.frame(gwas_data),
  block_annotation = geneid_df1,
  block_data =as.data.frame(Single_mat[,1:50]) ,
  ld_folder =ld_path,
  bootstrap_iters = 200
)
save(rolypoly_result,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/rolypoly_modeldata4.RData")
```

#!/usr/bin/sh
#PBS -N roly4
#PBS -q workq
#PBS -l nodes=node04
#PBS -l ncpus=1
#PBS -j oe
####mild

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scPagwas1.8/5.r

整合rolypoly结果：

```
rolypoly_score<-unlist(lapply(1:20,function(i){
load(paste0("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/rolypoly_modeldata",i,".RData"))
return(rolypoly_result$block_heritability_contribution)
}))
```

#### 整合结果到scPagwas：

```R
 
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_nk_v1.9.1.RData") 
#Pagwas$rolypoly_score<-rolypoly_score
rolypoly_score<-unlist(lapply(1:20,function(i){
load(paste0("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/rolypoly_modeldata",i,".RData"))
return(rolypoly_result$block_heritability_contribution)
}))

rolypoly_p<-unlist(lapply(1:20,function(i){
load(paste0("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/rolypoly_modeldata",i,".RData"))
return(rolypoly_result$bootstrap_results$bp_value[-1])
}))


Pagwas$rolypoly_score<-rolypoly_score
Pagwas$rolypoly_p<-rolypoly_p
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_nk_v1.9.1.RData")


```

## 和其他方法得结果比较

### 1.计算model数据的magma-scDRS

```R
library(SeuratDisk)
library(Seurat)

#sim_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_NK.rds")
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_nk_v1.9.1.RData")

DefaultAssay(Pagwas) <- "RNA"
SaveH5Seurat(Pagwas, "modeldata_NK_addata.h5seurat")
Convert("modeldata_NK_addata.h5seurat", dest="h5ad")
#write.csv(Single_data@meta.data,file="modeldata_lymph.metadata.csv")
#scDRS
i<-"monocytecount"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_magma_genes.RData"))

scPagwas_genes<-names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation,decreasing=T),])

magmatop1000<-intersect(magma_genes$symbol[1:1000],rownames(Pagwas))
scPagwastop1000<-scPagwas_genes[1:1000]
 magmatop500<-intersect(magma_genes$symbol[1:500],rownames(Pagwas))
scPagwastop500<-scPagwas_genes[1:500]
topgene<-list(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500)
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")

a<-data.frame(genes=c(paste(scPagwastop1000,collapse=","),paste(magmatop1000,collapse=","),paste(scPagwastop500,collapse=","),paste(magmatop500,collapse=",")))
rownames(a)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
write.csv(a,file=paste0("scPagwas.Model.magmatopgenes_scRDS",i,".csv"))


#rownames(a)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
a<-c("magmatop1000","NA",magmatop1000)
a<-matrix(a,nrow=1)
write.table(a,file=paste0("scPagwas.magma.Model.",i,".gmt"),append = T,col.names = F,row.names = F,quote = F,sep = "\t")
a<-c("scPagwastop1000","NA",scPagwastop1000)
a<-matrix(a,nrow=1)
write.table(a,file=paste0("scPagwas.magma.Model.",i,".gmt"),append = T,col.names = F,row.names = F,quote = F,sep = "\t")
a<-c("scPagwastop500","NA",scPagwastop500)
a<-matrix(a,nrow=1)
write.table(a,file=paste0("scPagwas.magma.Model.",i,".gmt"),append = T,col.names = F,row.names = F,quote = F,sep = "\t")
a<-c("magmatop500","NA",magmatop500)
a<-matrix(a,nrow=1)
write.table(a,file=paste0("scPagwas.magma.Model.",i,".gmt"),append = T,col.names = F,row.names = F,quote = F,sep = "\t")

################################计算scdrs
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
DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/"
H5AD_FILE = os.path.join(DATA_PATH, "modeldata_NK_addata.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
scdrs.preprocess(adata)
i="monocytecount"
df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scPagwas.Model.magmatopgenes_scRDS"+i+".csv", index_col=0)
df_gs = df_gs.loc[["scPagwastop1000","magmatop1000","scPagwastop500","magmatop500"],:]
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/"+i+".geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/monocytecount.geneset.gs")

gene_list  = df_gs['magmatop1000'][0]
gene_weight  = df_gs['magmatop1000'][1]
df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=20)
df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/monocytecount.1000magma.df_modelnk.csv", sep=",", index=False)

gene_list  = df_gs['magmatop500'][0]
gene_weight  = df_gs['magmatop500'][1]
df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=20)
df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/monocytecount.500magma.df_modelnk.csv", sep=",", index=False)

gene_list  = df_gs['scPagwastop1000'][0]
gene_weight  = df_gs['scPagwastop1000'][1]
df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=20)
df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/monocytecount.1000scPagwas.df_modelnk.csv", sep=",", index=False)

gene_list  = df_gs['scPagwastop500'][0]
gene_weight  = df_gs['scPagwastop500'][1]
df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=20)
df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/monocytecount.500scPagwas.df_modelnk.csv", sep=",", index=False)

```

### 2.整合结果到seruat格式中

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")
#load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_nk_v1.9.1.RData")

for(i in c(10)){
    j<-100*i
scPagwas_topgenes <- names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation, decreasing = T), ])[1:j]
Pagwas <- Seurat::AddModuleScore(Pagwas, assay = "RNA", list(scPagwas_topgenes), name = paste0("scPagwas.",j,"topgenes.Score"))
}

i<-"monocytecount"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_magma_genes.RData"))

for(i in c(1:10)){
    j<-100*i
magmatop1000<-intersect(magma_genes$symbol[1:j],rownames(Pagwas))
Pagwas <- Seurat::AddModuleScore(Pagwas, assay = "RNA", list(magmatop1000), name = paste0("magma.",j,"seruat.Score"))
}
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_nk_v1.9.1.RData")
library(scPagwas)
suppressMessages(library(Seurat))
library(org.Hs.eg.db)
library(dplyr)
 library(ggplot2)
library(ggtext)

library('ComplexHeatmap')
library(circlize)
magma1000_scDRS_re<-read.csv("monocytecount.1000magma.df_modelnk.csv")
Pagwas@meta.data$magma1000_scDRS_score<- rev(magma1000_scDRS_re$raw_score)
magma500_scDRS_re<-read.csv("monocytecount.500magma.df_modelnk.csv")
Pagwas@meta.data$magma500_scDRS_score<-rev(magma500_scDRS_re$raw_score)

scPagwas.1000_scDRS_re<-read.csv("monocytecount.1000scPagwas.df_modelnk.csv")
Pagwas@meta.data$scPagwas.1000_scDRS_score<- rev(scPagwas.1000_scDRS_re$raw_score)
scPagwas.500_scDRS_re<-read.csv("monocytecount.500scPagwas.df_modelnk.csv")
Pagwas@meta.data$scPagwas.500_scDRS_score<- rev(scPagwas.500_scDRS_re$raw_score)
save(Pagwas,file="Pagwas_monocytecount_nk_v1.9.1.RData")

```

### 3.UMAP得分映射

包括scPagwas和magma得得分映射

```R
###############umap图
library(Seurat)
library(SeuratObject)
library("dplyr")
library("foreach")
library("ggplot2")
library(reshape2)  # melt
library(gridExtra) # arrangeGrob
library(ggpubr)    # as_ggplot
library(cowplot)   # draw_plot_label
library(ggtext)
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")
load("Pagwas_monocytecount_nk_v1.9.1.RData")

#Pagwas$scPagwas.1000_scDRS_score<-rev(Pagwas$scPagwas.1000_scDRS_score)

all_fortify_can <- fortify.Seurat.umap(Pagwas)
all_fortify_can$scPagwas.1000_scDRS_score <- rev(Pagwas$scPagwas.1000_scDRS_score)

all_fortify_can$rolypoly_logp<- -log10(all_fortify_can$rolypoly_p)
    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =rolypoly_logp), size = 0.3, alpha = 1) +
        umap_theme() +
         scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("rolypoly log10pvalue")

pdf(file="./modeldata_outputv1.9.1/Umap.monocytecount_nk.rolypoly_logpe.pdf",width = 6, height = 6)
print(p1)
dev.off()

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =magma500_scDRS_score), size = 0.3, alpha = 1) +
        umap_theme() +
         scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("magma500_scDRS_score")

pdf(file="./modeldata_outputv1.9.1/Umap.monocytecount_nk.magma500_scDRS_score.pdf",width = 6, height = 6)
print(p1)
dev.off()

#all_fortify_can <- fortify.Seurat.umap(Pagwas)
    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =magma1000_scDRS_score), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("magma1000_scDRS_score")

pdf(file="./modeldata_outputv1.9.1/Umap.modelmonocytecount_nk.magma1000_scDRS_score.pdf",width = 6, height = 6)
print(p1)
dev.off()


    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas.topgenes.Score1), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("scPagwas.topgenes.Score1")

pdf(file="./modeldata_outputv1.9.1/Umap.modelmonocytecount_nk.scPagwas.topgenes.Score1.pdf",width = 6, height = 6)
print(p1)
dev.off()

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas.1000_scDRS_score), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("scPagwas.1000_scDRS")

pdf(file="./modeldata_outputv1.9.1/Umap.modelmonocytecount_nk.scPagwas.1000_scDRS_score.pdf",width = 6, height = 6)
print(p1)
dev.off()

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas.500topgenes.Score1), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("scPagwas.500topgenes.Score")

pdf(file="./modeldata_outputv1.9.1/Umap.modelmonocytecount_nk.scPagwas.500topgenes.Score.pdf",width = 6, height = 6)
print(p1)
dev.off()

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =magma.1000seruat.Score1), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("magma.1000seruat.Score")

pdf(file="./modeldata_outputv1.9.1/Umap.modelmonocytecount_nk.magma.1000seruat.Score.pdf",width = 6, height = 6)
print(p1)
dev.off()
```

### 4.计算其他得分：auccell， vision

```R
####################

library(testSctpa)
library(dplyr)
library(GSVA)
library(VISION)
library(AUCell)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")

counts = load_counts()
se_oj = CreateSeuratObject(counts)
se_oj = cal_PAS(seurat_object = se_oj,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = 'log',
              species = 'mouse', 
              pathway='kegg')

DefaultAssay(object = Pagwas) <- "RNA"

auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
gmt_file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scPagwas.magma.Model.monocytecount.gmt")

auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
Pagwas$aucell_magmatop1000<-auccell_df$magmatop1000 
Pagwas$aucell_scPagwastop1000<-auccell_df$scPagwastop1000 
Pagwas$aucell_scPagwastop500<-auccell_df$scPagwastop500 
Pagwas$aucell_magmatop500<-auccell_df$magmatop500
save(Pagwas,file="Pagwas_monocytecount_nk_v1.9.1.RData")

Visiondf = cal_PAS(seurat_object = Pagwas,
                       tool = 'Vision',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
gmt_file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scPagwas.magma.Model.monocytecount.gmt")
Visiondf<- as.data.frame(t(GetAssayData(Visiondf, slot="data", assay="PAS")))
Pagwas$Vision_magmatop1000<-Visiondf$magmatop1000 
Pagwas$Vision_scPagwastop1000<-Visiondf$scPagwastop1000 
Pagwas$Vision_scPagwastop500<-Visiondf$scPagwastop500 
Pagwas$Vision_magmatop500<-Visiondf$magmatop500
save(Pagwas,file="Pagwas_monocytecount_nk_v1.9.1.RData")

gsva_df<-gsva(data.matrix(Pagwas@assays$RNA@data), gset.idx.list=topgene,method="gsva")
gsva_df<-as.data.frame(t(gsva_df))
Pagwas$gsva_magmatop1000<-gsva_df$magmatop1000 
Pagwas$gsva_scPagwastop1000<-gsva_df$scPagwastop1000 
Pagwas$gsva_scPagwastop500<-gsva_df$scPagwastop500 
Pagwas$gsva_magmatop500<-gsva_df$magmatop500
save(Pagwas,file="Pagwas_monocytecount_nk_v1.9.1.RData")


```

### 5.AUC结果画图

```R
####################auc画图
library(circlize)
library(dplyr)
require(pROC)
require(ggplot2)
##不同数量的top基因的打分

load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_nk_v1.9.1.RData")
all_fortify_can <- fortify.Seurat.umap(Pagwas)

auc_list1<-lapply(all_fortify_can[,c(29:32,12,33:37)],function(x){
  roc(predictor=x,response=all_fortify_can$type)  
})
names(auc_list1)<-colnames(all_fortify_can)[c(29:32,12,33:37)]
 auc_l1<- unlist(lapply(1:length(auc_list1),function(x) round(as.numeric(auc_list1[[x]]["auc"]),3)))
 names(auc_list1)<- paste0(names(auc_list1),"(AUC=",auc_l1,")")
 pdf("AUC_modeldata.pagwas.100.1000.pdf",height = 6)
 ggroc(auc_list1, linetype = 2, size = 1,alpha=0.8)+
    ggtitle("ROC curve for different number of gene") + 
    theme_classic()+ggsci::scale_color_lancet()+
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
                 color="grey", linetype="dashed")
    #theme(panel.background=element_rect(fill="white",colour="blue"))
  dev.off()


#,"scPagwas.1000_scDRS_score","aucell_magmatop1000","aucell_scPagwastop1000","Vision_magmatop1000","Vision_scPagwastop1000","gsva_magmatop1000", "gsva_scPagwastop1000"

all_fortify_can$type[501:503]<-1
auc_list1<-lapply(all_fortify_can[,c("scPagwas.lm.score","scPagwas.topgenes.Score1","scPagwas.500topgenes.Score1","magma1000_scDRS_score","magma500_scDRS_score","rolypoly_logp")],function(x){
  roc(predictor=x,response=all_fortify_can$type)  
})
names(auc_list1)<-c("scPagwas.gPAS.score","scPagwas.top1000genes.Score","scPagwas.top500geness.Score","scDRS.top1000gene.score","scDRS.top500gene.score","rolypoly.logp")

 auc_l1<- unlist(lapply(1:length(auc_list1),function(x) round(as.numeric(auc_list1[[x]]["auc"]),3)))
 names(auc_list1)<- paste0(names(auc_list1),"(AUC=",auc_l1,")")
 pdf("AUC_modeldata.pagwas.methods.pdf",height = 5,width=8)
 ggroc(auc_list1, linetype = 2, size = 1,alpha=0.8)+
    ggtitle("ROC curve for different methods") + 
    theme_classic()+ggsci::scale_color_lancet()+
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
                 color="grey", linetype="dashed")
    #theme(panel.background=element_rect(fill="white",colour="blue"))
  dev.off()


auc_list1<-lapply(all_fortify_can[,c("scPagwas.topgenes.Score1","scPagwas.1000_scDRS_score","aucell_scPagwastop1000","Vision_scPagwastop1000","gsva_scPagwastop1000")],function(x){
  roc(predictor=x,response=all_fortify_can$type)  
})

names(auc_list1)<-c("scPagwas.topgenes.Score1","scPagwas.1000_scDRS_score","aucell_scPagwastop1000","Vision_scPagwastop1000","gsva_scPagwastop1000")
 auc_l1<- unlist(lapply(1:length(auc_list1),function(x) round(as.numeric(auc_list1[[x]]["auc"]),3)))
 names(auc_list1)<- c("scPagwas.seruat.Score","scPagwas.scDRS.score","scPagwas.aucell.score","scPagwas.Vision.score","scPagwas.gsva.score")
                        
 names(auc_list1)<- paste0(names(auc_list1),"(AUC=",auc_l1,")")
 pdf("AUC_modeldata.pagwas.methods_supplement1.pdf",height = 5,width=8)
 ggroc(auc_list1, linetype = 2, size = 1,alpha=0.8)+
    ggtitle("ROC curve for different methods") + 
    theme_classic()+ggsci::scale_color_lancet()+
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
                 color="grey", linetype="dashed")
    #theme(panel.background=element_rect(fill="white",colour="blue"))
  dev.off()

###################################
auc_list1<-lapply(all_fortify_can[,29:37],function(x){
  roc(predictor=x,response=all_fortify_can$type)  
})

names(auc_list1)<-colnames(all_fortify_can)[29:37]
 auc_l1<- unlist(lapply(1:length(auc_list1),function(x) round(as.numeric(auc_list1[[x]]["auc"]),3)))

                        
 names(auc_list1)<- paste0(names(auc_list1),"(AUC=",auc_l1,")")
 pdf("AUC_modeldata.pagwas.methods_supplement2.pdf",height = 5,width=8)
 ggroc(auc_list1, linetype = 2, size = 1,alpha=0.8)+
    ggtitle("ROC curve for different methods") + 
    theme_classic()+ggsci::scale_color_lancet()+
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
                 color="grey", linetype="dashed")
    #theme(panel.background=element_rect(fill="white",colour="blue"))
  dev.off()
#############################100genes auc
 
all_fortify_can$type[501:503]<-1
auc_list1<-lapply(all_fortify_can[,c("scPagwas.100topgenes.Score1","magma.100seruat.Score1")],function(x){
  roc(predictor=x,response=all_fortify_can$type)  
})
names(auc_list1)<-c("scPagwas.100.topgenes","magma.100.topgenes")

 auc_l1<- unlist(lapply(1:length(auc_list1),function(x) round(as.numeric(auc_list1[[x]]["auc"]),3)))
 names(auc_list1)<- paste0(names(auc_list1),"(AUC=",auc_l1,")")
 pdf("AUC_modeldata.100genescompare.pdf",height = 5,width=7)
 ggroc(auc_list1, linetype = 2, size = 1,alpha=0.8)+
    ggtitle("") + 
    theme_classic()+ggsci::scale_color_lancet()+
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
                 color="grey", linetype="dashed")
    #theme(panel.background=element_rect(fill="white",colour="blue"))
  dev.off()
#############################500genes auc
 
all_fortify_can$type[501:503]<-1
auc_list1<-lapply(all_fortify_can[,c("scPagwas.500topgenes.Score1","magma500_scDRS_score")],function(x){
  roc(predictor=x,response=all_fortify_can$type)  
})
names(auc_list1)<-c("scPagwas.500.topgenes","magma.500.topgenes")

 auc_l1<- unlist(lapply(1:length(auc_list1),function(x) round(as.numeric(auc_list1[[x]]["auc"]),3)))
 names(auc_list1)<- paste0(names(auc_list1),"(AUC=",auc_l1,")")
 pdf("AUC_modeldata.500genescompare.pdf",height = 5,width=7)
 ggroc(auc_list1, linetype = 2, size = 1,alpha=0.8)+
    ggtitle("") + 
    theme_classic()+ggsci::scale_color_lancet()+
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
                 color="grey", linetype="dashed")
    #theme(panel.background=element_rect(fill="white",colour="blue"))
  dev.off()
```

cd8细胞的模拟数据画图：

```

setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocyte_modeldata_cd8_v1.9.1.RData")
library(scPagwas)
 suppressMessages(library(Seurat))
library(org.Hs.eg.db)
library(dplyr)


library('ComplexHeatmap')
library(circlize)
magma1000_scDRS_re<-read.csv("monocytecount.1000magma.df_modelcd8.csv")
Pagwas@meta.data$magma1000_scDRS_score<-magma1000_scDRS_re$raw_score
magma500_scDRS_re<-read.csv("monocytecount.500magma.df_modelcd8.csv")
Pagwas@meta.data$magma500_scDRS_score<-magma500_scDRS_re$raw_score

meta.data=Pagwas@meta.data;
filecell="monocytecount";
#compare_type=c("sclm_allsnpscore");
celltypes= c("Mono","NK")
colors_celltypes= c("#ff7f0e","#1f77b4");
    
df<-meta.data

df<-df[order(df$magma500_scDRS_score,decreasing=T),]
pdf(paste0("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Figure_modelcd8_magma500_scDRS_score_",filecell,".rankplot.pdf"),height =7,width =2)
print(Heatmap(data.matrix(df$celltype),
        col =colors_celltypes,
        cluster_columns = F,
        cluster_rows = F,
        border=T,
        row_gap = unit(0.25, "mm"),
        show_parent_dend_line=T,
        name = "magma rank"
))

dev.off()

df<-df[order(df$scPagwas.lm.score,decreasing=T),]
pdf(paste0("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Figure_modelcd8_scPagwas.lm.score_",filecell,".rankplot.pdf"),height =7,width =2)
print(Heatmap(data.matrix(df$celltype),
        col =colors_celltypes,
        cluster_columns = F,
        cluster_rows = F,
        border=T,
        row_gap = unit(0.25, "mm"),
        show_parent_dend_line=T,
        name = "scPagwas rank"
))

dev.off()

#"#EE8572","#ECB390","#E5F4E7","#E5F4E7"

df<-df[order(df$rolypoly_score,decreasing=T),]
pdf(paste0("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Figure_rolypoly_score_",filecell,".rankplot.pdf"),height =7,width =2)
print(Heatmap(data.matrix(df$celltype),
        col =colors_celltypes,
        cluster_columns = F,
        cluster_rows = F,
        border=T,
        row_gap = unit(0.25, "mm"),
        show_parent_dend_line=T,
        name = "scPagwas rank"
))

dev.off()


```

## 子函数

```
umap_theme <- function() {
  theme_grey() %+replace%
    theme(
      panel.background = element_rect(fill = "white", colour = "black", size = 2),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_blank()
    )
}


fortify.Seurat.umap <- function(x) {
  xy1 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "umap"))
  colnames(xy1) <- c("UMAP_1", "UMAP_2")
  xy1$UMAP_1 <- as.numeric(xy1$UMAP_1)
  xy1$UMAP_2 <- as.numeric(xy1$UMAP_2)

  return(cbind(xy1, as.data.frame(x@meta.data)))
}


#' fortify.Seurat.tsne
#' @description set data frame to ggplot
#' @param x seruat
#' @export
#' @return

fortify.Seurat.tsne <- function(x) {
  xy2 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "tsne"))
  colnames(xy2) <- c("TSNE_1", "TSNE_2")
  xy2$TSNE_1 <- as.numeric(xy2$TSNE_1)
  xy2$TSNE_2 <- as.numeric(xy2$TSNE_2)

  return(cbind(xy2, as.data.frame(x@meta.data)))
}

```

