# AD 单细胞计算，基于之前的计算结果

## 1.GSE138852 scPagwas

```R
#导入之前的计算结果
> brain1
An object of class Seurat 
10850 features across 13214 samples within 1 assay 
Active assay: RNA (10850 features, 0 variable features)
> table(Idents(brain1))

  oligo    unID   astro     OPC  neuron    endo      mg doublet 
   7432     925    2171    1078     656      98     449     405 

setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test")
#load("/share/pub/qiuf/brain/01-data/AD/scPagwas_AD.RData")
library(dplyr)
library(irlba)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(stringr)
library(parallel)
library(glmnet)
library(GenomicRanges)
library(IRanges)
library(utils)
library(ggplot2)
library(ggpubr)
library(bigstatsr)
library(RMTstat)
library(gridExtra)
library(data.table)
library(bigmemory)
library(biganalytics)
library(bigreadr)

lapply(list.files("/share/pub/dengcy/GWAS_Multiomics/pagwas/R3.25/")[-1],function(x){
source(paste0("/share/pub/dengcy/GWAS_Multiomics/pagwas/R3.25/",x))
})
 
 setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test")
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/Genes_by_pathway_kegg.RData")
load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/block_annotation.RData")
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/chrom_ld.RData")
 
library(Seurat)

pagwas<-scPagwas_main(Pagwas = NULL,
                      gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data = "/share/pub/qiuf/brain/01-data/AD/GSE138852.rds",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
pagwas<-Celltype_heritability_contributions(pagwas,iters = 200)
saveRDS(pagwas,"/share/pub/dengcy/GWAS_Multiomics/ad_test/pagwas_ad.rds")

Bootstrap_P_Barplot(Pagwas=pagwas,
                    figurenames="pagwas_ad.pdf",
                    width = 5,
                    height = 7,
                    do_plot=T,
                    title = "ad scPagwas")

saveRDS(pagwas,"pagwas_ad.rds")

pagwas_ad<-readRDS("pagwas_ad.rds")
pagwas_ad<-Singlecell_heritability_contributions(Pagwas=pagwas_ad)
saveRDS(pagwas,"pagwas_ad.rds")
############################################



```

## 2.GSE138852 亚集数据跑scPagwas(结果删除)

```R
library(scPagwas)
library(Seurat)

#2 Cell Subtypes
data_AD <-readRDS("/share/pub/qiuf/brain/01-data/AD/GSE138852.rds")
data_AD <- NormalizeData(data_AD, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(data_AD)<-data_AD$oupSample.cellType
AD_final <- subset(data_AD,idents=c("oligo","astro","OPC","neuron","mg"))
AD_final <- AD_final[,AD_final$disease.state=="AD"]

> AD_final
An object of class Seurat 
10850 features across 5727 samples within 1 assay 
Active assay: RNA (10850 features, 0 variable features)


Idents(AD_final)<-AD_final$oupSample.cellType
saveRDS(AD_final,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/AD_5cellsubset.rds")

##############
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test")
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)

pagwas<-scPagwas_main(Pagwas = NULL,
                      gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data = "/share/pub/dengcy/GWAS_Multiomics/ad_test/AD_5cellsubset.rds",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(pagwas,file="pagwas_ad_5cellsubset.RData")

load("pagwas_ad_5cellsubset.RData")
pagwas<-Celltype_heritability_contributions(pagwas,iters = 200)
save(pagwas,file="pagwas_ad_5cellsubset.RData")

Bootstrap_P_Barplot(Pagwas=pagwas,
                    figurenames="pagwas_ad_5cellsubset.pdf",
                    width = 5,
                    height = 7,
                    do_plot=T,
                    title = "ad scPagwas 5cellsubset")
pagwas<-Singlecell_heritability_contributions(Pagwas=pagwas,split_n=5)
save(pagwas,file="pagwas_ad_5cellsubset.RData")


```

## 3.GSE160936验证数据跑scPagwas

/share/pub/qiuf/brain/01-data/AD/supplement_resolution0.2.rds

```R
library(scPagwas)
library(Seurat)
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)

setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test")
pagwas<-scPagwas_main(Pagwas = NULL,
                      gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data = "/share/pub/qiuf/brain/01-data/AD/supplement_resolution0.2.rds",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
pagwas<-Celltype_heritability_contributions(pagwas,iters = 200)
saveRDS(pagwas,file="pagwas_ad_GSE160936.rds")

Bootstrap_P_Barplot(Pagwas=pagwas,
                    figurenames="pagwas_ad_GSE160936.pdf",
                    width = 5,
                    height = 7,
                    do_plot=T,
                    title = "ad scPagwas GSE160936")
lapply(list.files("/home/guofj/scPagwas/R/"),function(x){
source(paste0("/home/guofj/scPagwas/R/",x))
})
setwd("/home/guofj/scPagwas/ADtest")
pagwas<-readRDS("pagwas_ad_GSE160936.rds")
#library(scPagwas) 
pagwas<-Singlecell_heritability_contributions(Pagwas=pagwas,split_n=5)
saveRDS(pagwas,file="pagwas_ad_GSE160936.rds")

#s2.sh
```

<img src="E:\OneDrive\GWAS_Multiomics\scriptmd\image-20220421103451904.png" width="40%" />



可视化测试：

```R
pagwas_ad<-readRDS("pagwas_ad.rds")
pdf("pagwas_ad_barplot4.20.pdf",width=5)
Bootstrap_P_Barplot(Pagwas=pagwas_ad,
                    figurenames = NULL,
                    width = 5,
                    height = 7,
                    do_plot=T,
                    title = "ad scPagwas")
dev.off()
pdf("AD_Bootstrap_estimate_Plot4.20.pdf",width=5)
Bootstrap_estimate_Plot(Pagwas=pagwas_ad,
                        #figurenames = "AD_Bootstrap_estimate_Plot4.20.PDF",

                        do_plot=T)
dev.off()

ad_genes<-sort(pagwas_ad$gene_heritability_correlation[,1],decreasing=T)[1:500]
```

通路可视化：

```R
BiocManager::install(c( "impute", "preprocessCore"))
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
install.packages("WGCNA", repos=site)

suppressMessages(require("WGCNA"))
suppressMessages(require("patchwork"))
suppressMessages(require("tidygraph"))
suppressMessages(require("ggraph"))
suppressMessages(require("igraph"))
#check the objects
pdf("AD_pathway_contribution_network4.20.pdf")
plot_pathway_contribution_network(
                  mat_datExpr=pagwas_ad$pca_cell_df,
                  vec_pathwaycontribution=pagwas_ad$Pathway_block_heritability,
                  vec_pathways_highlight=names(sort(pagwas_ad$Pathway_block_heritability,decreasing = T)[1:5]),
                  n_max_pathways=20,
                  igraph_algorithm = "drl",
                  fontface_labels="bold.italic",
                  color_edge = "#9D9D9D",
                  fontSize_label_lg=4,
                  fontSize_legend_lg=4,
                  fontSize_legend_xlg=4,
                  edge_thickness = 1,
                  do_plot=T
                  
  )
dev.off()

```

通路展示

```R
library(ggplot2)
library(grDevices)
library(stats)
library(FactoMineR)
library(scales)
library(reshape2)
library(ggdendro)
library(grImport2)
library(gridExtra)
library(grid)
library(sisal)
#load("E:/OneDrive/GWAS_Multiomics/dox/FlexDotPlot-master/FlexDotPlot-master/data/PBMC3K_example_data.rda")

Pagwas$pca_scCell_mat<- apply(Pagwas$pca_scCell_mat,2, function(x) (x - min(x)) / (max(x) - min(x)))

pca_scCell_mat<-Pagwas$pca_scCell_mat[colnames(Pagwas$Pathway_sclm_results),]

Pagwas$Pathway_sclm_mat<-Pagwas$Pathway_sclm_results * t(pca_scCell_mat)

rownames(Pagwas$Pathway_sclm_mat)<-colnames(Pagwas$pca_scCell_mat)
Pagwas$Celltype_anno$annotation<-as.factor(Pagwas$Celltype_anno$annotation)

celltype_pathway_scproportion<-tapply(Pagwas$Celltype_anno$cellnames,Pagwas$Celltype_anno$annotation,function(x){
  colMeans(Pagwas$Pathway_sclm_mat[x,])
})

celltype_pathway_scproportion<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),celltype_pathway_scproportion)

celltype_pathway_scproportion<-as.data.frame(apply(celltype_pathway_scproportion,2, function(x) (x - min(x)) / (max(x) - min(x))))

rownames(celltype_pathway_scproportion)<-levels(Pagwas$Celltype_anno$annotation)

Pathway_sclm_mat<-t(as(Pagwas$Pathway_sclm_mat,"matrix"))


ls_scCell<-CreateSeuratObject(
  Pathway_sclm_mat,
  project = "CreateSeuratObject",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = Pagwas$Celltype_anno
)

Idents(ls_scCell)<-Pagwas$Celltype_anno$annotation
ls_scCell_markerpa <- FindAllMarkers(object = ls_scCell)

pa_cellytpes<-lapply(unique(ls_scCell_markerpa$cluster),function(x){
  ls_scCell_markerpa[ls_scCell_markerpa$cluster==x,"gene"][1:3]
})

pca_scCell<-CreateSeuratObject(
  pca_scCell_mat,
  project = "CreateSeuratObject",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = Pagwas$Celltype_anno
)

Idents(pca_scCell)<-Pagwas$Celltype_anno$annotation
pca_scCell <- FindVariableFeatures(object = pca_scCell,nfeatures = dim(pca_scCell)[1]*0.1)
pca_scCell_markerpa <- FindAllMarkers(object = pca_scCell,min.pct = 0,logfc.threshold = 0,return.thresh = 0.1)

#pas<-VariableFeatures(object = pca_scCell)
pas<-VariableFeatures(object = pca_scCell)
#pas<-VariableFeatures(object = pca_scCell)

a1<-celltype_pathway_scproportion[,pas]
a2<-as.data.frame(t(Pagwas$pca_cell_df))[,pas]

a1$celltype<-rownames(a1)
a2$celltype<-rownames(a2)

gg_pa<-reshape2::melt(a1,id.vars="celltype",variable.name = "pathway",value.name = "heritability_proporion")

gg_pa2<-reshape2::melt(a2,id.vars="celltype",variable.name = "pathway",value.name = "pca")
gg_pa<-merge(gg_pa,gg_pa2)
gg_pa<-gg_pa[,c("pathway","celltype","heritability_proporion","pca")]
plot_sc_pa_contri_dot(data.to.plot = gg_pa, 
         size_var = "heritability_proporion", 
         col_var="heritability_proporion",
         shape.scale = 8,
         dend_x_var = c("heritability_proporion"),
         dist_method="euclidean", 
         hclust_method="ward.D"
         )

```

## 4.28重新跑AD数据

### GSE138852

```R
library(scPagwas)
library(Seurat)
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test/428test")
pagwas<-scPagwas_main(Pagwas = NULL,
                      gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     add_eqtls="OnlyTSS",
                      output.prefix="ADKEGG",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data = "/share/pub/qiuf/brain/01-data/AD/GSE138852.rds",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
saveRDS(pagwas,"GSE138852pagwas_ad_kegg.rds")

###############reactome
library(scPagwas)
library(Seurat)
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test/428test")
pagwas<-scPagwas_main(Pagwas = NULL,
                      gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     add_eqtls="OnlyTSS",
                      output.prefix="ADreactome",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data = "/share/pub/qiuf/brain/01-data/AD/GSE138852.rds",
                     Pathway_list=genes.by.reactome.pathway,
                     chrom_ld = chrom_ld)
saveRDS(pagwas,"GSE138852pagwas_ad_reactome.rds")

```

### GSE160936

```R
library(scPagwas)
library(Seurat)
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
#setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test/428test")
Pagwas<-readRDS("/home/guofj/scPagwas/ADtest/pagwas_ad_GSE160936.rds")
#pagwas<-scPagwas_main(Pagwas = NULL,
#                      gwas_data ="/share/pub/qiuf/brain/01-#data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
#                     add_eqtls="OnlyTSS",
#                      output.prefix="GSE160936_ADKEGG",
 #                    block_annotation = block_annotation,
#                     assay="RNA",
#                     Single_data ="/share/pub/qiuf/brain/01-#data/AD/supplement_resolution0.2.rds",
#                     Pathway_list=Genes_by_pathway_kegg,
#                     chrom_ld = chrom_ld)
#saveRDS(pagwas,"GSE160936pagwas_ad_kegg.rds")
setwd("/home/guofj/scPagwas/ADtest")
Pagwas<-scPagwas_main(Pagwas = Pagwas,
                     add_eqtls="OnlyTSS",
                      gwas_data ="/home/guofj/scPagwas/ADtest/ad_Pagwas_gwas.txt",
                     output.prefix="GSE160936_ADKEGG",
                     block_annotation = block_annotation,
                     assay="RNA",
                     split_n=6,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
saveRDS(Pagwas,"GSE160936pagwas_ad_kegg.rds")



###############reactome
library(scPagwas)
library(Seurat)
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test/428test")
load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
pagwas<-scPagwas_main(Pagwas = NULL,
                      gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     add_eqtls="OnlyTSS",
                      output.prefix="GSE160936_ADreactome",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data = "/share/pub/qiuf/brain/01-data/AD/supplement_resolution0.2.rds",
                      split_n=5,
                     Pathway_list=genes.by.reactome.pathway,
                     chrom_ld = chrom_ld)
saveRDS(pagwas,"GSE160936pagwas_ad_reactome.rds")
```

## 5.6号重新跑AD v1.7

1.GSE138852 scPagwas

```
 setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test")
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/Genes_by_pathway_kegg.RData")
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/block_annotation.RData")
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/chrom_ld.RData")
library(Seurat)
library(scPagwas)
Single_data = readRDS("/share/pub/qiuf/brain/01-data/AD/GSE138852.rds")
pagwas<-scPagwas_main(Pagwas = NULL,
                      gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data = "/share/pub/qiuf/brain/01-data/AD/GSE138852.rds",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
saveRDS(pagwas,"/share/pub/dengcy/GWAS_Multiomics/ad_test/pagwas_ad.rds")

Bootstrap_P_Barplot(Pagwas=pagwas,
                    figurenames="pagwas_ad_GSE138852.pdf",
                    width = 5,
                    height = 7,
                    do_plot=T,
                    title = "GSE138852 ad scPagwas")
saveRDS(pagwas,"pagwas_ad_GSE138852.rds")


```

### AD疾病GSE138852亚组样本作为软件测试数据

由于测试数据的gwas不能太大，这里选择prune后的数据作为分析数据。

An object of class Seurat 
10850 features across 6673 samples within 1 assay 
Active assay: RNA (10850 features, 0 variable features)

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



#### 1.prune gwas

```
 library(scPagwas)
 library(Seurat)
 library(parallel)
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
 ##############真实例子：
 Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="E:/RPakage/scPagwas/inst/extdata/AD_prune_gwas_data.txt",
                     output.prefix="adPrunetest",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data ="E:/RPakage/scPagwas/inst/extdata/GSE138852_ad.rds",
                     split_n=1,
                     ncores=5,
                     nfeatures =NULL,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     SimpleResult=F)
save(Pagwas,file="E:/OneDrive/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_Prune_adsubset_kegg.RData")
```

<img src="E:\OneDrive\GWAS_Multiomics\scriptmd\figures\image-20220506204517047.png" width="60%" />

#### 2.全gwas

```R
export OPENBLAS_NUM_THREADS=1
############AD疾病亚组样本
 library(scPagwas)
 suppressMessages(library(Seurat))
 library(parallel)
 suppressMessages(library("dplyr"))
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)

 Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     output.prefix="test",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE138852_ad.rds",
                     split_n=1,
                     ncores=20,
                     nfeatures =NULL,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     SimpleResult=F)
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_adsubset_kegg.RData")
#!/usr/bin/sh
#PBS -N ad_test4
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=20
#PBS -j oe
####mild

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/4.r

load("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_adsubset_kegg.RData")
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test")
Bootstrap_P_Barplot(p_results=Pagwas$bootstrap_results$bp_value[-1],
                                p_names=rownames(Pagwas$bootstrap_results)[-1],
                                title = "GSE138852",
                                figurenames = "barplot_GSE138852_adsubset_kegg.pdf",
                                width = 5,
                                height = 7,
                                do_plot=F)
```



<img src="E:\OneDrive\GWAS_Multiomics\scriptmd\figures\image-20220507092158597.png" width="30%" />

#### 3.全gwas+reactome

```R
 library(scPagwas)
 suppressMessages(library(Seurat))
 library(parallel)
 suppressMessages(library("dplyr"))
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)

 Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     output.prefix="test",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE138852_ad.rds",
                     split_n=1,
                     ncores=5,
                     nfeatures =NULL,
                     Pathway_list=genes.by.reactome.pathway,
                     chrom_ld = chrom_ld,
                     SimpleResult=F)
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_adsubset_reactome.RData")

#!/usr/bin/sh
#PBS -N ad_test2
#PBS -q workq
#PBS -l nodes=node04
#PBS -l ncpus=5
#PBS -j oe
####mild

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/2.r

 library(scPagwas)
 suppressMessages(library(Seurat))
 library(parallel)
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test")
load("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_adsubset_reactome.RData")
Bootstrap_P_Barplot(p_results=Pagwas$bootstrap_results$bp_value[-1],
                                p_names=rownames(Pagwas$bootstrap_results)[-1],
                                title = "GSE138852 reactome",
                                figurenames = "barplot_GSE138852_adsubset_reactome.pdf",
                                width = 5,
                                height = 7,
                                do_plot=F)

```

<img src="E:\OneDrive\GWAS_Multiomics\scriptmd\figures\image-20220511092703572.png" width="30%" />

### GSE138852所有样本

#### kegg通路结果

```R
#export OPENBLAS_NUM_THREADS=1
############AD疾病亚组样本
library(scPagwas)
suppressMessages(library(Seurat))
library(parallel)
suppressMessages(library("dplyr"))
data(Genes_by_pathway_kegg)
#gene annotation files.
data(block_annotation)
#LD data
data(chrom_ld)

Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     output.prefix="test",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE138852.rds",
                     split_n=1,
                     ncores=5,
                     nfeatures =NULL,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     SimpleResult=F)
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_kegg.RData")

#!/usr/bin/sh
#PBS -N ad_test1
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=5
#PBS -j oe
####mild

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/1.r


 library(scPagwas)
 suppressMessages(library(Seurat))
 library(parallel)
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test")
load("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_kegg.RData")
Bootstrap_P_Barplot(p_results=Pagwas$bootstrap_results$bp_value[-1],
                                p_names=rownames(Pagwas$bootstrap_results)[-1],
                                title = "GSE138852 kegg all samples",
                                figurenames = "barplot_GSE138852_kegg.pdf",
                                width = 5,
                                height = 7,
                                do_plot=F)

```

<img src="E:\OneDrive\GWAS_Multiomics\scriptmd\figures\image-20220511092343675.png" width="30%" />

#### 对结果进行整合，符合v1.7.3版本得结果格式：

```R
 suppressMessages(library(Seurat))
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test")
load("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_kegg.RData")
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE138852.rds")

 scPagwas_pathway = SeuratObject::CreateAssayObject(data = t(data.matrix(Pagwas$Pathway_single_results))
                                                    
  scPagwas_pca = SeuratObject::CreateAssayObject(data = Pagwas$pca_scCell_mat)
  scPagwas_lm = SeuratObject::CreateAssayObject(data = t(data.matrix(Pagwas$Pathway_sclm_results)))
                                                    
Single_data<-Single_data[,colnames(Pagwas$data_mat)]
                                                    
  Single_data[['scPagwasPaHeritability']] = scPagwas_pathway
  Single_data[['scPagwasPaPca']] = scPagwas_pca
  Single_data[['scPagwaslmHeritability']] = scPagwas_lm

  scPagwas_topgenes<-names(Pagwas$gene_heritability_correlation[order(Pagwas$gene_heritability_correlation,decreasing=T),])[1:1000]
  Single_data<-AddModuleScore(Single_data,assay="RNA",list(scPagwas_topgenes),name=c("scPagwas.topgenes.Score"))
  Single_data$scPagwas.lm.score<-Pagwas$scPagwas_score[rownames(Pagwas$Celltype_anno)]
  Single_data$Cells.lm.rankPvalue<-Pagwas$CellsrankPvalue[rownames(Pagwas$Celltype_anno),"pValueHigh"]
  
    Pagwas$CellsrankPvalue$adj.p<-p.adjust(Pagwas$CellsrankPvalue$pValueHigh,
                                             method ="bonferroni" )
                                             
  Single_data$Cells.lm.adjp<-Pagwas$CellsrankPvalue[rownames(Pagwas$Celltype_anno),"adj.p"]
  
Pagwas[c("VariableFeatures","merge_scexpr","snp_gene_df","data_mat","rawPathway_list",
"scPagwas_score","CellsrankPvalue","Celltype_anno","pca_cell_df","Pathway_lm_results",
"pca_scCell_mat","Pathway_sclm_results","Pathway_single_results")]<-NULL
class(Pagwas)<-"list"
  Single_data@misc<-Pagwas
save(Single_data,file="Pagwas_seu_GSE138852_kegg.RData")

pdf(file = "GSE138852_TSNE.pdf", height = 7, width = 7)                                   
DimPlot(Single_data,group.by="oupSample.cellType",reduction="tsne",pt.size=0.5,label = TRUE, repel=TRUE,label.size = 4)+ 
umap_theme()+ labs(x="TSNE",y="")+
ggtitle("GSE138852")+
        scale_colour_manual(name = "celltypes", values = color_scanpy_viridis28) +
        theme(aspect.ratio=1)
                                                            
dev.off()

##################
Pagwas_fortify <- fortify.Seurat.tsne(Single_data)

library(ggpubr)
library(ggplot2)
library(ggtext)                               
plot2 <- ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = TSNE_1, y = TSNE_2,color =scPagwas.topgenes.Score1), size = 0.5, alpha = 1) +
      umap_theme() +
      scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D",
                             midpoint = median(Pagwas_fortify$scPagwas.topgenes.Score1))+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      ggtitle("GSE138852 scPagwas.topgenes.Score")
                                          
      pdf(file = "GSE138852.topgenes.Score_.tsne.pdf", height = 7, width = 7)
      print(plot2)
      dev.off()
                                                    
                                                    
    Pagwas_fortify$p_thre<-rep("non",nrow(Pagwas_fortify))
    Pagwas_fortify$p_thre[which(Pagwas_fortify$Cells.lm.adjp<0.05)]<-"significant"
    plot3 <-  ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = TSNE_1, y = TSNE_2,color =p_thre), size =0.5, alpha = 0.7) +
      umap_theme() +
      scale_color_manual(values = c("non" = "#CDD0CB", "significant" = "#E45826"))+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      ggtitle("significant cells")

  plot3 <- ggplot() +
   geom_point(data = Pagwas_fortify[order(Pagwas_fortify$Cells.lm.adjp <0.01),],
              aes(x = TSNE_1, y = TSNE_2), size = 0.5, alpha = 0.8, color = "gray") +
   umap_theme() +
   theme(aspect.ratio=1) +
   theme(legend.text=element_markdown(size=14),
         legend.title=element_text(size=14)) +
   guides(colour = guide_legend(override.aes = list(size=3))) +
   geom_point(data = Pagwas_fortify[(Pagwas_fortify$Cells.lm.adjp <0.01),], aes(x = TSNE_1, y = TSNE_2), size = 0.5, alpha = 0.8, color = "#E45826")+
   theme(aspect.ratio=1,
         legend.text=element_markdown(size=16),
         legend.title=element_text(size=16))
                                           
      pdf(file = "GSE138852.adjp_tsne.pdf", height = 7, width = 7)
      print(plot3)
      dev.off()

     plot4 <-  ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = TSNE_1, y = TSNE_2,color =scPagwas.lm.score), size =0.5, alpha = 1) +
      umap_theme() +
      scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D")+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      ggtitle("GSE138852")
   
      pdf(file = "GSE138852.lm_tsne.pdf", height = 7, width = 7)
      print(plot4)
      dev.off()
```

#### 重新计算cellrankpvalue

```R
library(scPagwas)
suppressMessages(library(Seurat))
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test")
#load("Pagwas_seu_GSE138852_kegg.RData")
load("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_kegg.RData")
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R524/sub_functions.R")

scPagwas_genes<-names(Pagwas$gene_heritability_correlation[order(Pagwas$gene_heritability_correlation,decreasing=T),])[1:1000]
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE138852.rds")
data.mat<-GetAssayData(Single_data,assay = "RNA")
data.mat<-data.mat[scPagwas_genes,]

p_cell <- rankPvalue(t(data.matrix(data.mat)),
                     pValueMethod = "scale",
                    calculateQvalue = T)

Single_data$pValueScale<-p_cell$pValueHighScale
Single_data$qValueScale<-p_cell$qValueHighScale

Pagwas_fortify <- fortify.Seurat.tsne(Single_data)

library(ggpubr)
library(ggplot2)
library(ggtext) 

Pagwas_fortify$pValueScale <0.05

plot3 <- ggplot() +
   geom_point(data = Pagwas_fortify[Pagwas_fortify$qValueScale >0.05,],
              aes(x = TSNE_1, y = TSNE_2), size = 0.5, alpha = 0.8, color = "gray") +
   umap_theme() +
   theme(aspect.ratio=1) +
   theme(legend.text=element_markdown(size=14),
         legend.title=element_text(size=14)) +
   guides(colour = guide_legend(override.aes = list(size=3))) +
   geom_point(data = Pagwas_fortify[(Pagwas_fortify$qValueScale <0.05),], aes(x = TSNE_1, y = TSNE_2), size = 0.5, alpha = 0.8, color = "#E45826")+
   theme(aspect.ratio=1,
         legend.text=element_markdown(size=16),
         legend.title=element_text(size=16))
                                           
      pdf(file = "GSE138852.cellp2_tsne.pdf", height = 7, width = 7)
      print(plot3)
      dev.off()

load("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_seu_GSE138852_kegg.RData")
Single_data$pValueScale<-p_cell$pValueHighScale
Single_data$qValueScale<-p_cell$qValueHighScale
save(Single_data,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_seu_GSE138852_kegg.RData")
```

#### 重新计算v1.9.0

export OPENBLAS_NUM_THREADS=1

```R
 library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
#
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test/531test")
suppressMessages(library("dplyr"))

Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE138852.rds",
                     output.prefix="test",
                     output.dirs="GSE138852AD_outputv1.9",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     seruat_return=T,
                     celltype=F,
                     ncores = 2,
                     split_n=1)
save(Pagwas,file="Pagwas_seu_GSE138852_v1.9.RData")



```

#!/usr/bin/sh
#PBS -N ad_test1
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=5
#PBS -j oe
####mild

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/ad_test/531test/1.r



#### reactome

```R
library(scPagwas)
 suppressMessages(library(Seurat))
  library(parallel)
 suppressMessages(library("dplyr"))
load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/reduce_genes.by.reactome.pathway.RData")
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test")
 Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     output.prefix="test",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE138852.rds",
                     split_n=1,
                     ncores=10,
                     nfeatures =NULL,
                     Pathway_list=reduce_genes.by.reactome.pathway,
                     chrom_ld = chrom_ld,
                     SimpleResult=F)
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_reactome.RData")

Bootstrap_P_Barplot(p_results=Pagwas$bootstrap_results$bp_value[-1],
                                p_names=rownames(Pagwas$bootstrap_results)[-1],
                                title = "GSE138852 reactome all samples",
                                figurenames = "barplot_GSE138852_reactome.pdf",
                                width = 5,
                                height = 7,
                                do_plot=F)
```

2test 

```R
library(scPagwas)
 suppressMessages(library(Seurat))
  library(parallel)
 suppressMessages(library("dplyr"))
 #data(Genes_by_pathway_kegg)
load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)

 Pagwas<-scPagwas_main(
     Pagwas = NULL;
                     gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt";
                     output.prefix="GSE138852";
                     add_eqtls="OnlyTSS";
                     block_annotation = block_annotation;
                     assay="RNA";
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE138852.rds";
                     split_n=1;
                     marg=5000;
                     ncores=20;
                     nfeatures =NULL;
                     Pathway_list=genes.by.reactome.pathway;
                     chrom_ld = chrom_ld;
                     SimpleResult=T;
 )


save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_reactome2.RData")
load("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_reactome2.RData")
Bootstrap_P_Barplot(p_results=Pagwas$bootstrap_results$bp_value[-1],
                                p_names=rownames(Pagwas$bootstrap_results)[-1],
                                title = "GSE138852 reactome2",
                                figurenames = "barplot_GSE138852_reactome2.pdf",
                                width = 5,
                                height = 7,
                                do_plot=F)
```



### GSE160936



head(single_cell@meta.data)

```
single_cell<-readRDS("/share/pub/qiuf/brain/01-data/AD/supplement_resolution0.2.rds")
table(Idents(single_cell))
Effector CD8+T cell   Memory CD8+T cell    Naive CD8+T cell    Naive CD4+T cell 
               7245                3127                 766               18525 
 Memory CD4+T cell        CD14+monocyte                  NK        Naive B cell 
               7034               44577               10922                5230 
      CD16+monocyte                  DC            Platelet     CD34+Progenitor 
               3634                1160                1425                 830 
      Mature B cell 
                  9
lapply(as.vector(unique(Idents(severe_all))[2:12]),function(x){
  a<-severe_all[,Idents(severe_all)==x] 
  saveRDS(a,file=paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/severe_",x,".rds"))
})

```



#### kegg

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
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE160936.rds",
                     singlecell=F,
                     split_n=4,
                     ncores=10,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     SimpleResult=T)
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.16test/Pagwas_GSE160936_kegg_celltypes.RData")
Bootstrap_P_Barplot(p_results=Pagwas$bootstrap_results$bp_value[-1],
                                p_names=rownames(Pagwas$bootstrap_results)[-1],
                                title = "GSE160936 kegg",
                                figurenames = "barplot_GSE160936_kegg.pdf",
                                width = 5,
                                height = 7,
                                do_plot=F)
#!/usr/bin/sh
#PBS -N ad_test1
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=5
#PBS -j oe
####mild

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/ad_test/5.16test/1.r

* Start to link gwas and pathway block annotations for 319 pathways!
  |==================================================                    |  72%
 *** caught bus error ***
address 0x7f8dbfb4a000, cause 'non-existent physical address'

Traceback:
 1: replaceMat(x$address_rw, i, j, value)
 2: replace_matrix(x, transform_ind(i, n), transform_ind(j, m), value)
 3: `[<-`(`*tmp*`, , value = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,/var/spool/pbs/mom_priv/jobs/277968.mgt.SC: line 11: 369955 Bus error               (core dumped) Rscript /share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/1.r
```



#### prune gwas

```R
library(scPagwas)
 suppressMessages(library(Seurat))
  library(parallel)
 suppressMessages(library("dplyr"))
setwd("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test")
load("Pagwas_GSE160936_prune_kegg.RData")
    Pagwas[c("Pathway_ld_gwas_data","VariableFeatures","merge_scexpr",
             "rawPathway_list",
             "snp_gene_df")]<-NULL
Bootstrap_P_Barplot(p_results=Pagwas$bootstrap_results$bp_value[-1],
                                p_names=rownames(Pagwas$bootstrap_results)[-1],
                                title = "GSE160936 kegg prunegwas",
                                figurenames = "barplot_GSE160936_kegg_prunegwas.pdf",
                                width = 5,
                                height = 7,
                                do_plot=F)

```

<img src="E:\OneDrive\GWAS_Multiomics\scriptmd\figures\image-20220511090331451.png" width="30%" />

#### reactome

```R
 library(scPagwas)
 suppressMessages(library(Seurat))
  library(parallel)
 suppressMessages(library("dplyr"))
 #data(Genes_by_pathway_kegg)
load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)

 Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     output.prefix="GSE160936",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE160936.rds",
                     split_n=3,
                       singlecell=F,
                     ncores=10,
                     Pathway_list=genes.by.reactome.pathway,
                     chrom_ld = chrom_ld,
                     SimpleResult=F)
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE160936_reactome_celltypes.RData")

Bootstrap_P_Barplot(p_results=Pagwas$bootstrap_results$bp_value[-1],
                                p_names=rownames(Pagwas$bootstrap_results)[-1],
                                title = "GSE160936 reactome",
                                figurenames = "/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/barplot_GSE160936_reactome.pdf",
                                width = 5,
                                height = 7,
                                do_plot=F)
#!/usr/bin/sh
#PBS -N ad_test2
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=2
#PBS -j oe
####mild

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/ad_test/5.16test/2.r
```







