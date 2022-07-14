# Covid19 traits for PBMC scRNA-seq data for scPagwas

### Severe

> covid
> An object of class Seurat 
> 37985 features across 104484 samples within 2 assays 
> Active assay: SCT (17803 features, 0 variable features)
>  1 other assay present: RNA



### moderate

An object of class Seurat 
37985 features across 204508 samples within 2 assays 
Active assay: SCT (17803 features, 0 variable features)
 1 other assay present: RNA



### mild

An object of class Seurat 
37985 features across 122686 samples within 2 assays 
Active assay: SCT (17803 features, 0 variable features)
 1 other assay present: RNA



### Normal

An object of class Seurat 
37985 features across 37775 samples within 2 assays 
Active assay: SCT (17803 features, 0 variable features)
 1 other assay present: RNA

## scPagwas v1.9

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

### severe

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

sub celltypes

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
celltype<-c("Effector CD8+T cell","Memory CD8+T cell","Naive CD8+T cell" ,"Naive CD4+T cell" ,"Memory CD4+T cell","CD14+monocyte","NK","Naive B cell", "CD16+monocyte","DC","Platelet","CD34+Progenitor" )
i<-celltype[1]

#i<-"Naive CD8+T cell"
# setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test5.6")
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
save(Pagwas,file=paste0("scpagwas_severe_",i".RData"))#i<-"Naive CD8+T cell"
# setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test5.6")


############cd14monocyte
celltype<-c("Effector CD8+T cell","Memory CD8+T cell","Naive CD8+T cell" ,"Naive CD4+T cell" ,"Memory CD4+T cell","CD14+monocyte","NK","Naive B cell", "CD16+monocyte","DC","Platelet","CD34+Progenitor" )
i<-celltype[6]
j<-1

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
```

#### Integrate the result：

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
i<-celltype[1]
load("scpagwas_severe_",i,j,".RData")
Pagwas_integrate <- merge_pagwas(Pagwas_list = list(Pagwas1,Pagwas2),n_topgenes = 1000)   


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

### mild

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
data(Genes_by_pathway_kegg)
 data(block_annotation)
 data(chrom_ld)
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

setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
data(Genes_by_pathway_kegg)
 data(block_annotation)
 data(chrom_ld)
 
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/mild_all.rds",
                     singlecell=F,
                      output.prefix="test",
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
Pagwas$Pathway_ld_gwas_data<-NULL
save(Pagwas,file="scpagwas_mild.v1.71.RData")


###########
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
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/mild_all.rds",
                     singlecell=F,
                      output.prefix="mild",
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
Pagwas$Pathway_ld_gwas_data<-NULL
save(Pagwas,file="scpagwas_mild.v1.73.RData")
```

#### based on v1.7.3 to v1.9.1

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
celltype<-c("Effector CD8+T cell","Memory CD8+T cell","Naive CD8+T cell" ,"Naive CD4+T cell" ,"Memory CD4+T cell ","CD14+monocyte","NK","Naive B cell", "CD16+monocyte","DC","Platelet","CD34+Progenitor" )
i<-celltype[1]
load(paste0("/share2/pub/jiangdp/jiangdp/COVID/data/scpagwas_mild_",i,".RData"))
Single_data<-readRDS(paste0("/share2/pub/jiangdp/jiangdp/COVID/data/mild_",i,".rds"))
Pagwas1<-rerun_Pagwas(Single_data=Single_data,Pagwas=Pagwas,assay="RNA",n_topgenes=1000)

for(i in celltype[2:12]){  
    print(i)
load(paste0("/share2/pub/jiangdp/jiangdp/COVID/data/scpagwas_mild_",i,".RData"))
Single_data<-readRDS(paste0("/share2/pub/jiangdp/jiangdp/COVID/data/mild_",i,".rds"))
Pagwas<-rerun_Pagwas(Single_data=Single_data,Pagwas=Pagwas,assay="RNA",n_topgenes=1000)
Pagwas1<-merge(Pagwas1,y= Pagwas, project = "merged", merge.data = TRUE)   
}
save(Pagwas1,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_mild.v1.91.RData")
###############subfunction
rerun_Pagwas<-function(Single_data=Single_data,Pagwas,assay="RNA",n_topgenes=1000){
     class(Pagwas)<-"list"  
    scPagwas_topgenes <- names(Pagwas$gene_heritability_correlation[order(Pagwas$gene_heritability_correlation, decreasing = T), ])[1:n_topgenes]
    Single_data <- Seurat::AddModuleScore(Single_data, assay = assay, list(scPagwas_topgenes), name = c("scPagwas.topgenes.Score"))
      Pagwas[c(
      "VariableFeatures", "merge_scexpr","snp_gene_df",
      "rawPathway_list","data_mat"
    )] <- NULL

    message("* Get rankPvalue for each single cell")
    CellScalepValue <- rankPvalue(datS = t(data.matrix(GetAssayData(Single_data, assay = assay)[scPagwas_topgenes, ])), pValueMethod = "scale")
     Pagwas[c("snp_gene_df",
               "Pathway_sclm_results","CellScalepValue",
               "scPagwas.topgenes.Score"
               )] <- NULL
      #Pagwas[c()] <- NULL
      scPagwas_pathway <- SeuratObject::CreateAssayObject(data = t(data.matrix(Pagwas$Pathway_single_results)))
      scPagwas_pca <- SeuratObject::CreateAssayObject(data = Pagwas$pca_scCell_mat)

      Single_data[["scPagwasPaHeritability"]] <- scPagwas_pathway
      Single_data[["scPagwasPaPca"]] <- scPagwas_pca

      Single_data$scPagwas.lm.score <- Pagwas$scPagwas_score[rownames(Pagwas$Celltype_anno)]
      Single_data$CellScalepValue <- CellScalepValue[rownames(Pagwas$Celltype_anno), "pValueHighScale"]
      Single_data$CellScaleqValue <- CellScalepValue[rownames(Pagwas$Celltype_anno), "qValueHighScale"]
      Pagwas[ c("scPagwas.score","Celltype_anno",
                "Pathway_single_results","pca_scCell_mat")] <- NULL
   
      Single_data@misc<-Pagwas
    return(Single_data)
}

```

Integrate the result：

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
load("scpagwas_severe.v1.71.RData")
names(Pagwas)
pahtways <- rownames(Pagwas$pca_scCell_mat)

celltype<-c("Effector CD8+T cell","Memory CD8+T cell","Naive CD8+T cell" ,"Naive CD4+T cell" ,"Memory CD4+T cell","CD14+monocyte1","CD14+monocyte2","CD14+monocyte3","CD14+monocyte4","CD14+monocyte5","NK","Naive B cell", "CD16+monocyte","DC","Platelet","CD34+Progenitor" )

Pathway_single_results<-list()
Pathway_sclm_results<-list()
cellids<-list()

for(i in celltype){

load(paste0("scpagwas_mild_",i,".RData"))
x<-as.data.frame(data.matrix(Pagwas$Pathway_single_results))
a<- setdiff(pahtways,colnames(x))
x[,a]<-0  
Pathway_single_results[[i]]<-x[,pahtways]
cellids[[i]]<-colnames(Pagwas$pca_scCell_mat)
    
y<-as.data.frame(data.matrix(Pagwas$Pathway_sclm_results))
a<- setdiff(pahtways,colnames(y))
y[,a]<-0  
Pathway_sclm_results[[i]]<-y[,pahtways] 
}

cellids<-unlist(cellids)

Pathway_single_results<-data.matrix(bigreadr::rbind_df(Pathway_single_results))

rownames(Pathway_single_results)<-cellids
colnames(Pathway_single_results)<-pahtways

Pathway_sclm_results<-data.matrix(bigreadr::rbind_df(Pathway_sclm_results))
#,"dgCMatrix")
rownames(Pathway_sclm_results)<-cellids
colnames(Pathway_sclm_results)<-pahtways
save(Pathway_sclm_results,Pathway_single_results,file="pathway_mat_severe.RData")


load("scpagwas_severe.v1.71.RData")

Pagwas$Pathway_sclm_results<-Pathway_sclm_results[colnames(Pagwas$data_mat),]
Pagwas$Pathway_single_results<-Pathway_single_results[colnames(Pagwas$data_mat),]
dim(Pagwas$Pathway_sclm_results)
dim(Pagwas$Pathway_single_results)
save(Pagwas,file="scpagwas_severe.v1.71.RData")

library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
  Pagwas<-scPagwas_perform_score(Pagwas=Pagwas,
                                 remove_outlier=TRUE)
  Pagwas <- scGet_gene_heritability_correlation(Pagwas=Pagwas)
save(Pagwas,file="scpagwas_severe.v1.71.RData")


```

### moderate

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
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/moderate_all.rds",
                      output.prefix="test",
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


###########
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
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/moderate_all.rds",
                     singlecell=F,
                      output.prefix="moderate",
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
Pagwas$Pathway_ld_gwas_data<-NULL
save(Pagwas,file="scpagwas_moderate.v1.73.RData")
```

#### based on v1.7.3 to v1.9.1

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
celltype<-c("Effector CD8+T cell","Memory CD8+T cell","Naive CD8+T cell" ,"Naive CD4+T cell" ,"Memory CD4+T cell ","CD14+monocyte","NK","Naive B cell", "CD16+monocyte","DC","Platelet","CD34+Progenitor" )
i<-celltype[1]
load(paste0("/share2/pub/jiangdp/jiangdp/COVID/data/scpagwas_moderate_",i,".RData"))
Single_data<-readRDS(paste0("/share2/pub/jiangdp/jiangdp/COVID/data/moderate_",i,".rds"))
Pagwas1<-rerun_Pagwas(Single_data=Single_data,Pagwas=Pagwas,assay="RNA",n_topgenes=1000)
#rm(Pagwas)
for(i in celltype[2:12]){  
print(i)
load(paste0("/share2/pub/jiangdp/jiangdp/COVID/data/scpagwas_moderate_",i,".RData"))
Single_data<-readRDS(paste0("/share2/pub/jiangdp/jiangdp/COVID/data/moderate_",i,".rds"))
Pagwas<-rerun_Pagwas(Single_data=Single_data,Pagwas=Pagwas,assay="RNA",n_topgenes=1000)
Pagwas1<-merge(Pagwas1,y= Pagwas, project = "merged", merge.data = TRUE)  
    #rm(Pagwas)
}
save(Pagwas1,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_moderate.v1.91.RData")
```

### normal

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

#### based on v1.7.3 to v1.9.1

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
celltype<-c("Effector CD8+T cell","Memory CD8+T cell","Naive CD8+T cell" ,"Naive CD4+T cell" ,"Memory CD4+T cell ","CD14+monocyte","NK","Naive B cell", "CD16+monocyte","DC","Platelet","CD34+Progenitor" )
i<-celltype[1]
load(paste0("/share2/pub/jiangdp/jiangdp/COVID/data/scpagwas_Normal_",i,".RData"))
Single_data<-readRDS(paste0("/share2/pub/jiangdp/jiangdp/COVID/data/Normall_",i,".rds"))
Pagwas1<-rerun_Pagwas(Single_data=Single_data,Pagwas=Pagwas,assay="RNA",n_topgenes=1000)

for(i in celltype[2:12]){  
    print(i)
load(paste0("/share2/pub/jiangdp/jiangdp/COVID/data/scpagwas_Normal_",i,".RData"))
Single_data<-readRDS(paste0("/share2/pub/jiangdp/jiangdp/COVID/data/Normal_",i,".rds"))
Pagwas<-rerun_Pagwas(Single_data=Single_data,Pagwas=Pagwas,assay="RNA",n_topgenes=1000)
Pagwas1<-merge(Pagwas1,y= Pagwas, project = "merged", merge.data = TRUE)   
}
save(Pagwas1,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_Normal.v1.91.RData")
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



## New 64w covid data：

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD datas
 data(chrom_ld)
#/share2/pub/zhenggw/zhenggw/COVID_combined/covid_64w_raw.rds
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/")
a<-readRDS("/share2/pub/zhenggw/zhenggw/COVID/covid_64w_new.rds")

ci<-c("Severe","Asymptomatic","Mild","Healthy","Critical","Moderate","LPS_10hours","LPS_90mins")
i<-ci[8]
a1<-a[,a$Status_on_day_collection_summary==i]
Idents(a1)<-a1$full_clustering
rm(a)
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = a1,
                      output.prefix=paste0(i,"2"),
                     singlecell=F,
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/covid2_",i,"_Pagwas.RData"))
Bootstrap_P_Barplot(p_results=Pagwas$bootstrap_results$bp_value[-1],
                                p_names=rownames(Pagwas$bootstrap_results)[-1],
                                title = i,
                                figurenames = paste0("barplot_",i,"_kegg2.pdf"),
                                width = 5,
                                height = 7,
                                do_plot=F)
```

