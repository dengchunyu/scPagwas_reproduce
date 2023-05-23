HCA single celll data



## 从GEO中下载

地址：E:/OneDrive/GWAS_Multiomics/HCLdata/GSE134355_RAW.tar

```R

library("stringr") 
library("readr")
library(Seurat)

setwd("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata")
files<-list.files("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata/HCAadult")
##利用edgeR包读入一下数据（主动选择部分数据）因为数据量太大
Seurat_list<-lapply(files,function(x){
    count<-read_table(x)
    gene<-unlist(count[,1])
    count<-as(as(count[,-1],"matrix"),"dgCMatrix")
    rownames(count)<-gene
    a<- strsplit(x,split = "_",fixed=T)[[1]][2]
    a<-str_replace_all(a,"-","_")
    Seurat_object <- CreateSeuratObject(
               counts = count, 
               min.cells = 3, 
               min.features = 200)
    Idents(Seurat_object)<-rep(a,ncol(Seurat_object))
    return(Seurat_object)
                })
    
Seurat_data<-merge(Seurat_list[[1]],Seurat_list[[2]])
for(i in 3:length(files)){
    Seurat_data<-merge(Seurat_data,Seurat_list[[i]])
}
rm(Seurat_list)
Seurat_data <- NormalizeData(Seurat_data, normalization.method = "LogNormalize", scale.factor = 10000)
saveRDS(Seurat_data,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata/HCA_tissues.rds")
```

scPagwas

```
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 #Input pathway gene list, you can construct with youself.

 ####先计算平均表达'
  suppressMessages(library(Seurat))
 HCA_tissues<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata/HCA_tissues.rds")
 ##合并部分细胞
 Idents(HCA_tissues)<-unlist(lapply(as.vector(Idents(HCA_tissues)),function(x){
  return(substr(x,1,nchar(x)-1))
 }))

b<-as.vector(Idents(HCA_tissues))
b[which(b=="Adult_Kidney4_")]<-"Adult_Kidney"
b[which(b=="Adult_Liver1_")]<-"Adult_Liver"
b[which(b=="Adult_Liver4_")]<-"Adult_Liver"
b[which(b=="Adult_Lung3_")]<-"Adult_Lung"
b[which(b=="Adult_Peripheral_Blood3_")]<-"Adult_Peripheral_Blood"
b[which(b=="Adult_Peripheral_Blood4_")]<-"Adult_Peripheral_Blood"
b[which(b=="Adult_Stomach3_")]<-"Adult_Stomach"
b[which(b=="Adult_Transverse_Colon2_")]<-"Adult_Transverse_Colon"
table(b)
Idents(HCA_tissues)<-b
saveRDS(HCA_tissues,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata/HCA_tissues.rds")

#############
 suppressMessages(library(Seurat))
  library(scPagwas)
  
 HCA_tissues<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata/HCA_tissues.rds")
HCA_tissues_scexpr <- Seurat::AverageExpression(HCA_tissues,assays="RNA")[["RNA"]]
save(HCA_tissues_scexpr,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCA_tissues_scexpr.RData")
rm(HCA_tissues_scexpr)

Pagwas <- list();
class(Pagwas) <- 'Pagwas'
Pagwas <- Single_data_input(Pagwas=Pagwas,
                                assay="RNA",
Single_data=/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata/HCA_tissues.rds,
                                Pathway_list=Genes_by_pathway_kegg)
  
  Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,n.cores=1,
                                 Pathway_list=Genes_by_pathway_kegg
                                 )
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Pagwas_HCA_tissues.RData")
```

BMMC

```
setwd("E:/OneDrive/SingleCell/data/PBMCscATAC-seq")
scRNA_Healthy_Hema<-readRDS("scRNA-Healthy-Hematopoiesis-191120.rds")
table(scRNA_Healthy_Hema@colData$BioClassification)

library(scRNAseq)
library(SingleCellExperiment)
library(Seurat)
library("stringr") 
counts <- assay(scRNA_Healthy_Hema, "counts")
Seu_Healthy_Hema <- CreateSeuratObject(
               counts = counts, 
    meta.data=as.data.frame(colData(scRNA_Healthy_Hema)),
               min.cells = 3, 
               min.features = 200)
##查看分类
table(scRNA_Healthy_Hema@colData$BioClassification)
        01_HSC 02_Early.Eryth  03_Late.Eryth  04_Early.Baso    05_CMP.LMPP 
          1425           1653            446            111           2260 
      06_CLP.1         07_GMP    08_GMP.Neut         09_pDC         10_cDC 
           903           2097           1050            544            325 
11_CD14.Mono.1 12_CD14.Mono.2   13_CD16.Mono         14_Unk       15_CLP.2 
          1800           4222            292            520            377 
      16_Pre.B           17_B      18_Plasma       19_CD8.N      20_CD4.N1 
           710           1711             62           1521           2470 
     21_CD4.N2       22_CD4.M      23_CD8.EM      24_CD8.CM          25_NK 
          2364           3539            796           2080           2143 
        26_Unk 
           161 

Idents(Seu_Healthy_Hema)<-scRNA_Healthy_Hema@colData$BioClassification
Seu_Healthy_Hema<-ScaleData(Seu_Healthy_Hema)
saveRDS(Seu_Healthy_Hema,file="Seu_Healthy_Hema.rds")
Seu_Healthy<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Healthy_Hema.rds")
length(table(Idents(Seu_Healthy)))
```

pbmc

```
# 自己安装  mojaveazure/seurat-disk 这个GitHub包：
#remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
library(patchwork)
#～～～～～开始读数据～～～～～
##h5ad是python的Scanpy读取文件格式，需要转换
#～～～～读取adipose～～～～
#Convert('/share/pub/dengcy/GWAS_Multiomics/singlecelldata/COVID19_Healthy.h5ad', "h5ad",overwrite = TRUE,assay = "RNA")

NM_Healthy_pbmc <- LoadH5Seurat("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/COVID19_Healthy.h5seurat")

## Normalizing the data
NM_Healthy_pbmc <- NormalizeData(NM_Healthy_pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

## Identify the 2000 most highly variable genes
NM_Healthy_pbmc <- FindVariableFeatures(NM_Healthy_pbmc, selection.method = "vst", nfeatures = 2000)

## In addition we scale the data
all.genes <- rownames(NM_Healthy_pbmc)
NM_Healthy_pbmc <- ScaleData(NM_Healthy_pbmc, features = all.genes)

NM_Healthy_pbmc <- RunPCA(NM_Healthy_pbmc, features = VariableFeatures(object = NM_Healthy_pbmc), verbose = FALSE)
NM_Healthy_pbmc <- FindNeighbors(NM_Healthy_pbmc, dims = 1:10, verbose = FALSE)
NM_Healthy_pbmc <- FindClusters(NM_Healthy_pbmc, resolution = 0.5, verbose = FALSE)
NM_Healthy_pbmc <- RunUMAP(scRNA, dims = 1:10, umap.method = "uwot", metric = "cosine")
saveRDS(NM_Healthy_pbmc,file="NM_Healthy_pbmc.rds")
table(NM_Healthy_pbmc$seurat_clusters)
phe=NM_Healthy_pbmc@meta.data
save(NM_Healthy_pbmc,file = 'phe_NM_Healthy_pbmc.Rdata')
```

