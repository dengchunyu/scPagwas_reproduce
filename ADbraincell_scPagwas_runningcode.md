# AD gwas data based on bain single cell

## Single cell data

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

## GSE138852

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

## GSE160936

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
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE160936.rds",
                     singlecell=F,
                     ncores=10,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ldT)
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.16test/Pagwas_GSE160936_kegg_celltypes.RData")
Bootstrap_P_Barplot(p_results=Pagwas$bootstrap_results$bp_value[-1],
                                p_names=rownames(Pagwas$bootstrap_results)[-1],
                                title = "GSE160936 kegg",
                                figurenames = "barplot_GSE160936_kegg.pdf",
                                width = 5,
                                height = 7,
                                do_plot=F)

```

#### reactome

```R
 library(scPagwas)
 suppressMessages(library(Seurat))
  library(parallel)
 suppressMessages(library("dplyr"))
load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")

 Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share/pub/qiuf/brain/01-data/GWAS/00_IEU_GWAS/02_neurodegenerative_disorders/AD/Pagwas.txt",
                     output.prefix="GSE160936",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/GSE160936.rds",
                       singlecell=F,
                     ncores=10,
                     Pathway_list=genes.by.reactome.pathway,
                     chrom_ld = chrom_ldF)
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







