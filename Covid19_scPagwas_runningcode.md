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

save rds files

```R
load("/share/pub/dengcy/GWAS_Multiomics/test/covid19/Single_data_severe.RData")
saveRDS(Single_data_severe,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/severe_all.rds")

load("/share/pub/dengcy/GWAS_Multiomics/test/covid19/Single_data_moderate.RData")
saveRDS(Single_data_moderate,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/moderate_all.rds")
load("/share/pub/dengcy/GWAS_Multiomics/test/covid19/Single_data_mild.RData")
saveRDS(Single_data_mild,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/mild_all.rds")
load("/share/pub/dengcy/GWAS_Multiomics/test/covid19/Single_data_normal.RData")
saveRDS(Single_data_normal,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/normal_all.rds")
```

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
                     Single_data = paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/mild_",i,".rds"),
                      output.prefix="COVID19severe_kegg",
                     output.dirs=i,
                      ncores=5,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file=paste0("scpagwas_mild_",i,".RData"))
}

```

#### Integrate the result

```R
Pagwas_l<-lapply(celltype,function(i){
load(paste0("scpagwas_mild_",i,".RData"))
 return(Pagwas)   
})
Pagwas <- merge_pagwas(Pagwas_list = Pagwas_l,n_topgenes = 1000)
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
celltype<-c("Effector CD8+T cell","Memory CD8+T cell","Naive CD8+T cell" ,"Naive CD4+T cell" ,"Memory CD4+T cell","CD14+monocyte","NK","Naive B cell", "CD16+monocyte","DC","Platelet","CD34+Progenitor" )

for(i in celltype){
Pagwas<-scPagwas_main(Pagwas = NULL,
                       gwas_data ="/share/pub/dengcy/GWAS_Multiomics/test/covid19/COVID19_GWAS.txt",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/moderate_",i,".rds"),
                      output.prefix="COVID19moderate_kegg",
                     output.dirs=i,
                      ncores=5,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file=paste0("scpagwas_moderate_",i,".RData"))
}
```

#### Integrate the result

```R
Pagwas_l<-lapply(celltype,function(i){
load(paste0("scpagwas_moderate_",i,".RData"))
 return(Pagwas)   
})
Pagwas <- merge_pagwas(Pagwas_list = Pagwas_l,n_topgenes = 1000)
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



## Visualize the result

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



