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

#### v1.9.0

export OPENBLAS_NUM_THREADS=1

```R
 library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))

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
                     celltype=T,
                     ncores = 2)
save(Pagwas,file="Pagwas_seu_GSE138852_v1.9.RData")
```

## 2. GSE160936

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









