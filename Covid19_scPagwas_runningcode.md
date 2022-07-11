# 新冠数据四个阶段跑scPagwas

## Severe

> covid
> An object of class Seurat 
> 37985 features across 104484 samples within 2 assays 
> Active assay: SCT (17803 features, 0 variable features)
>  1 other assay present: RNA

```R
 library(scPagwas)
 #1.2.0
 suppressMessages(library(Seurat))
 suppressWarnings(library(SOAR))
 suppressMessages(library("dplyr"))
 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test3.23")
 Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                     add_eqtls="OnlyTSS",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/severe_all.rds",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file="scpagwas_severe.RData")
Pagwas<-Singlecell_heritability_contributions(Pagwas,bignumber=100000)
save(Pagwas,file="scpagwas_severe.RData")


##########
#!/usr/bin/sh
#PBS -N severe
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
####mild
source activate R403
Rscript /share/pub/dengcy/GWAS_Multiomics/test/covid19/test3.23/Severe.r

####mayunlong
source activate R403
Rscript /home/guofj/scPagwas/covid19/1.r

library(scPagwas)
setwd("/home/guofj/scPagwas/covid19")
load("scpagwas_severe.RData")
Pagwas<-Singlecell_heritability_contributions(Pagwas=Pagwas,split_n=5)
save(Pagwas,file="scpagwas_severe.RData")

```

## moderate

An object of class Seurat 
37985 features across 204508 samples within 2 assays 
Active assay: SCT (17803 features, 0 variable features)
 1 other assay present: RNA

```R
 library(scPagwas)
 #1.2.0
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test3.23")
 Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                     add_eqtls="OnlyTSS",
                         assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/moderate_all.rds",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)    

 save(Pagwas,file="scpagwas_moderate.RData")
Pagwas<-Singlecell_heritability_contributions(Pagwas)
  save(Pagwas,file="scpagwas_moderate.RData")

  ##########
#!/usr/bin/sh
#PBS -N read_ad
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
####mild
source activate R403
Rscript /share/pub/dengcy/GWAS_Multiomics/test/covid19/test3.23/moderate.r

```

## mild

An object of class Seurat 
37985 features across 122686 samples within 2 assays 
Active assay: SCT (17803 features, 0 variable features)
 1 other assay present: RNA

```R
 library(scPagwas)
 suppressMessages(library(Seurat))
 suppressWarnings(library(SOAR))
 suppressMessages(library("dplyr"))
 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
#lapply(list.files("/share/pub/dengcy/GWAS_Multiomics/pagwas/R3.25/")[2:17],function(x){
#source(paste0("/share/pub/dengcy/GWAS_Multiomics/pagwas/R3.25/",x))
#})
# load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/Genes_by_pathway_kegg.RData")
 #load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/block_annotation.RData")
## load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/chrom_ld.RData")

#source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R3.25/Single_data_input.R")
 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test3.23")
 #Pagwas<-scPagwas_main(Pagwas = NULL,
 #                    gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
  #                   add_eqtls="OnlyTSS",
  #                       assay="RNA",
  #                   block_annotation = block_annotation,
  #                   Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/mild_all.rds",
 #                    Pathway_list=Genes_by_pathway_kegg,
 #                    chrom_ld = chrom_ld)
load("scpagwas_mild.RData")
Pagwas<-Singlecell_heritability_contributions(Pagwas)
save(Pagwas,file="scpagwas_mild.RData")

  ###################
#!/usr/bin/sh
#PBS -N mild
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
####mild
source activate R403
Rscript /share/pub/dengcy/GWAS_Multiomics/test/covid19/test3.23/mild.r

```

## Normal

An object of class Seurat 
37985 features across 37775 samples within 2 assays 
Active assay: SCT (17803 features, 0 variable features)
 1 other assay present: RNA

```R
library(scPagwas)
 #1.2.0
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test3.23")
 Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/Normal.rds",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
  save(Pagwas,file="scpagwas_Normal.RData")
Pagwas<-Singlecell_heritability_contributions(Pagwas,bignumber=100000)
  save(Pagwas,file="scpagwas_Normal.RData")
  
#!/usr/bin/sh
#PBS -N Normal
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
####mild
source activate R403
Rscript /share/pub/dengcy/GWAS_Multiomics/test/covid19/test3.23/Normal.r
```

## 重新计算单细胞数据v1.6.1

### severe

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
load("/share/pub/dengcy/GWAS_Multiomics/test/covid19/Single_data_severe.RData")
 saveRDS(Single_data_severe,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/severe_all.rds")
rm(Single_data_severe)
 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test428")
Pagwas<-scPagwas_main(Pagwas = NULL,
                       gwas_data ="/share/pub/dengcy/GWAS_Multiomics/test/covid19/COVID19_GWAS.txt",
                     add_eqtls="OnlyTSS",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = "/share/pub/dengcy/GWAS_Multiomics/test/covid19/severe_all.rds",
                      output.prefix="COVID19severe_kegg",
                      split_n=5,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file="scpagwas_severe.RData")
```

### moderate

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test428")
pagwas<-scPagwas_main(Pagwas = NULL,
                       gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                     add_eqtls="OnlyTSS",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/moderate_all.rds",
                      output.prefix="COVID19moderate_kegg",
                      split_n=5,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file="scpagwas_moderate.RData")
```

### mild

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test428")
pagwas<-scPagwas_main(Pagwas = NULL,
                       gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                     add_eqtls="OnlyTSS",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/mild_all.rds",
                      output.prefix="COVID19mild_kegg",
                      split_n=5,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file="scpagwas_mild.RData")
```

### normal

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test428")
pagwas<-scPagwas_main(Pagwas = NULL,
                       gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                     add_eqtls="OnlyTSS",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/Normal.rds",
                      output.prefix="COVID19Normal_kegg",
                      split_n=5,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file="scpagwas_Normal.RData")

```



## 重新计算单细胞数据v1.9

```
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
data(Genes_by_pathway_kegg)
 data(block_annotation)
 data(chrom_ld)

 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test5.6")

severe_all<-readRDS("/share2/pub/jiangdp/jiangdp/COVID/data/severe_all.rds")
table(Idents(severe_all))
Effector CD8+T cell   Memory CD8+T cell    Naive CD8+T cell    Naive CD4+T cell 
               7245                3127                 766               18525 
 Memory CD4+T cell        CD14+monocyte                  NK        Naive B cell 
               7034               44577               10922                5230 
      CD16+monocyte                  DC            Platelet     CD34+Progenitor 
               3634                1160                1425                 830 
      Mature B cell 
                  9

 
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

#!/usr/bin/sh
#PBS -N covid2
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=5
#PBS -j oe
####mild

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/test/covid19/test5.6/2.r

* Start to link gwas and pathway block annotations for pathways!
  |==================================================================    |  95%Error: Inconsistency between size of backingfile and dimensions.
Execution halted
```

分细胞类型

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
data(Genes_by_pathway_kegg)
 data(block_annotation)
 data(chrom_ld)
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


#!/usr/bin/sh
#PBS -N 12
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=2
#PBS -j oe
####mild

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16/12.r

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

#### severe数据分细胞类型跑之后整合：

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
#load("scpagwas_severe.v1.9.1.RData")
#names(Pagwas)
#pahtways <- rownames(Pagwas$pca_scCell_mat)

celltype<-c("Effector CD8+T cell","Memory CD8+T cell","Naive CD8+T cell" ,"Naive CD4+T cell" ,"Memory CD4+T cell","CD14+monocyte1","CD14+monocyte2","CD14+monocyte3","CD14+monocyte4","CD14+monocyte5","NK","Naive B cell", "CD16+monocyte","DC","Platelet","CD34+Progenitor" )
i<-celltype[1]
load(paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16/scpagwas_severe_",i,".RData"))
Pagwas@misc<-list()
Pagwas1<-Pagwas
 
for(i in celltype[2:16]){
  load(paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16/scpagwas_severe_",i,".RData"))
    Pagwas@misc<-list()
    Pagwas1<-merge(Pagwas1,y= Pagwas, project = "merged", merge.data = TRUE)
}

save(Pagwas1,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16/Pagwas_merge.RData")

load("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16/Pagwas_merge.RData")
    scPagwas.lm.score <- data.matrix(Pagwas1$scPagwas.lm.score)
    sparse_cor <- corSparse(
      X = t(as_matrix(GetAssayData(Pagwas1, slot = "data", assay = "RNA"))),
      Y = scPagwas.lm.score
    )
 rownames(sparse_cor) <- rownames(GetAssayData(Pagwas1, slot = "data", assay = "RNA"))
  colnames(sparse_cor) <- "gene_heritability_correlation"
  sparse_cor[is.nan(sparse_cor)] <- 0
  Pagwas1@misc$gene_heritability_correlation <- sparse_cor

    scPagwas_topgenes <- names(sparse_cor[order(sparse_cor, decreasing = T), ])[1:1000]
    Pagwas1 <- Seurat::AddModuleScore(Pagwas1, assay =  "RNA", list(scPagwas_topgenes), name = c("scPagwas.topgenes.Score"))
    CellScalepValue <- rankPvalue(datS = t(data.matrix(GetAssayData(Pagwas1, assay = "RNA")[scPagwas_topgenes, ])), pValueMethod = "scale")
    Pagwas1$CellScalepValue <- CellScalepValue[,"pValueHighScale"]
      Pagwas1$CellScaleqValue <- CellScalepValue[,"qValueHighScale"]
 save(Pagwas1,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/scpagwas_severe.v1.91.RData")

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
###################################子函数
corSparse <- function(X, Y) {
  # X <-as_matrix(X)
  n <- nrow(X)
  muX <- colMeans(X)

  stopifnot(nrow(X) == nrow(Y))

  muY <- colMeans(Y)
  covmat <- (as.matrix(crossprod(X, Y)) - n * tcrossprod(muX, muY))
  sdvecX <- sqrt((colSums(X^2) - n * muX^2))
  sdvecY <- sqrt((colSums(Y^2) - n * muY^2))
  cormat <- covmat / tcrossprod(sdvecX, sdvecY)
  # cormat[is.nan(cormat),1]<-0
  return(cormat)
}

as_matrix <- function(mat) {
  tmp <- matrix(data = 0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  row_pos <- mat@i + 1
  col_pos <- findInterval(seq(mat@x) - 1, mat@p[-1]) + 1
  val <- mat@x
  for (i in seq_along(val)) {
    tmp[row_pos[i], col_pos[i]] <- val[i]
  }
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
```

#### naive cd8t细胞的排秩比较：

```R
library('ComplexHeatmap')
library(circlize)
 suppressMessages(library(Seurat))
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
##导入scPagwas结果
load("scpagwas_severe.v1.71.RData")
scPagwas_genes<-names(Pagwas$gene_heritability_correlation[order(Pagwas$gene_heritability_correlation,decreasing=T),])
save(scPagwas_genes,file="severe.scPagwas_rankgenes.RData")

 Pagwas[c("VariableFeatures","merge_scexpr","pca_scCell_mat","Pathway_sclm_results",
           "Pathway_single_results","snp_gene_df",
           "data_mat","rawPathway_list")]<-NULL

Single_data<-readRDS("/share2/pub/jiangdp/jiangdp/COVID/data/severe_all.rds")

Celltype_anno<-Pagwas$Celltype_anno
Celltype_anno$scPagwas_score<-Pagwas$scPagwas_score[rownames(Pagwas$Celltype_anno)]

##计算magma得分
covid19magma<-read.csv("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16/covid19magma_genes.csv")
a<-intersect(covid19magma$Gene_name[1:1000],rownames(Single_data))
Single_data<-AddModuleScore(Single_data,assay="RNA",list(a),name=c("magma.topgenes.Score"))

#################

#magmatop1000<-intersect(magma_genes$symbol[1:1000],rownames(Single_data))
scPagwastop1000<-scPagwas_genes[1:1000]
# magmatop500<-intersect(magma_genes$symbol[1:500],rownames(Single_data))
scPagwastop500<-scPagwas_genes[1:500]

Single_data<-AddModuleScore(Single_data,assay="RNA",list(scPagwastop1000),name=c("scPagwas.topgenes.Score"))
Single_data$celltypes<-Idents(Single_data)

GroundtruthRankplot(meta.data=Single_data@meta.data,filecell="Naive CD8+T",
                    celltypes= c("Naive CD8+T cell","DC"),
         colors_celltypes=c("#E5F4E7","#EE8572"))
#################函数
GroundtruthRankplot<-function(meta.data,
filecell="monocytecount",
celltypes= c("12_CD14.Mono.2","13_CD16.Mono","24_CD8.CM","25_NK"), 
colors_celltypes=c("#EE8572","#ECB390","#E5F4E7","#E5F4E7")){

    df<-meta.data[meta.data$celltypes %in% celltypes,]
df<-df[order(df$scPagwas.topgenes.Score,decreasing=T),]
pdf(paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16/Figure_1000scPagwas_",filecell,".rankplot.pdf"),
    height =7,width =2)
print(Heatmap(data.matrix(df$celltypes),
        col =colors_celltypes,
        #left_annotation = ha, 
        cluster_columns = F,
        cluster_rows = F,
        #color_space="HLS",
        border=T,
        row_gap = unit(0.25, "mm"),
        show_parent_dend_line=T,
        name = "scPagwas rank"
        #row_order =order(rdf$types),
        #show_row_names=T,
        #show_column_names=T
        #row_split=rdf$phenotypes
))
dev.off()

df<-df[order(df$magma.topgenes.Score1,decreasing=T),]
pdf(paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16/Figure_1000magma_",filecell,".rankplot.pdf"),
    height =7,width =2)
print(Heatmap(data.matrix(df$celltypes),
        col =colors_celltypes,
        #left_annotation = ha, 
        cluster_columns = F,
        cluster_rows = F,
        #color_space="HLS",
        border=T,
        row_gap = unit(0.25, "mm"),
        show_parent_dend_line=T,
        name = "scPagwas rank"
        #row_order =order(rdf$types),
        #show_row_names=T,
        #show_column_names=T
        #row_split=rdf$phenotypes
))
dev.off()          
    
}

```

![image-20220521154921679](E:\OneDrive\GWAS_Multiomics\scriptmd\figures\image-20220521154921679.png)

#### 比例图：

```R
library(ggpubr)
library(ggplot2)

filecell="Naive CD8+T"
celltypes= c("Naive CD8+T cell","DC")
colors_celltypes=c("#E5F4E7","#EE8572")
meta.data<-Single_data@meta.data
df<-meta.data[meta.data$celltypes %in% celltypes,] 
df1<-df[order(df$scPagwas.topgenes.Score1,decreasing=T),]

n<-nrow(df1)
b<-rep(1,n)
b[1:(0.25*n)]<-1
b[(0.25*n+1):(0.5*n)]<-2
b[(0.5*n+1):(0.75*n)]<-3
b[(0.75*n+1):n]<-4
df1$rg<-b
df1$celltype<-rep("Correct",nrow(df1))
df1$celltype[df1$celltypes %in% "DC"]<-"non"

print(table(data.frame(df1$celltype,df1$rg)))
t1<-table(data.frame(df1$celltype,df1$rg))
percent1<-t1[,1]/sum(t1[,1])
#  Correct        non 
#0.97722567 0.02277433 

df <- data.frame(
  group = names(percent1),
  value = percent1)

p1<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4") )

percent2<-t1[,2]/sum(t1[,2])
# Correct       non 
#0.6029106 0.3970894
df <- data.frame(
  group = names(percent2),
  value = percent2)

p2<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4"))

percent3<-t1[,3]/sum(t1[,3])
#    Correct         non 
#0.008316008 0.991683992 

df <- data.frame(
  group = names(percent3),
  value = percent3)

p3<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4"))

percent4<-t1[,4]/sum(t1[,4])
#Correct     non 
#      0       1

df <- data.frame(
  group = names(percent4),
  value = percent4)

p4<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4"))
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
pdf(paste0(filecell,"percent.1000.pdf"),width = 10,height = 10)
ggpubr::ggarrange(p1,p2,p3,p4,nrow = 1,ncol = 4)
dev.off()

############magma
filecell="Naive CD8+T"
celltypes= c("Naive CD8+T cell","DC")
colors_celltypes=c("#E5F4E7","#EE8572")
meta.data<-Single_data@meta.data
df<-meta.data[meta.data$celltypes %in% celltypes,] 
df1<-df[order(df$magma.topgenes.Score1,decreasing=T),]

n<-nrow(df1)
b<-rep(1,n)
b[1:(0.25*n)]<-1
b[(0.25*n+1):(0.5*n)]<-2
b[(0.5*n+1):(0.75*n)]<-3
b[(0.75*n+1):n]<-4
df1$rg<-b
df1$celltype<-rep("Correct",nrow(df1))
df1$celltype[df1$celltypes %in% "DC"]<-"non"

print(table(data.frame(df1$celltype,df1$rg)))
t1<-table(data.frame(df1$celltype,df1$rg))
percent1<-t1[,1]/sum(t1[,1])
# Correct       non 
#0.3581781 0.6418219

df <- data.frame(
  group = names(percent1),
  value = percent1)

p1<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4") )

percent2<-t1[,2]/sum(t1[,2])
#  Correct       non 
#0.3679834 0.6320166
df <- data.frame(
  group = names(percent2),
  value = percent2)

p2<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4"))

percent3<-t1[,3]/sum(t1[,3])
#  Correct       non 
#0.4220374 0.5779626 

df <- data.frame(
  group = names(percent3),
  value = percent3)

p3<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4"))

percent4<-t1[,4]/sum(t1[,4])
#  Correct       non 
#0.4428274 0.5571726

df <- data.frame(
  group = names(percent4),
  value = percent4)

p4<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4"))
setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/testsplit5.16")
pdf(paste0(filecell,"percent.magma.1000.pdf"),width = 10,height = 10)
ggpubr::ggarrange(p1,p2,p3,p4,nrow = 1,ncol = 4)
dev.off()
```

![image-20220521154853718](E:\OneDrive\GWAS_Multiomics\scriptmd\figures\image-20220521154853718.png)

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
Effector CD8+T cell   Memory CD8+T cell    Naive CD8+T cell    Naive CD4+T cell 
               7245                3127                 766               18525 
 Memory CD4+T cell        CD14+monocyte                  NK        Naive B cell 
               7034               44577               10922                5230 
      CD16+monocyte                  DC            Platelet     CD34+Progenitor 
               3634                1160                1425                 830 
      Mature B cell 
                  9
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

#!/usr/bin/sh
#PBS -N mild
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=5
#PBS -j oe
####mild

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/test/covid19/test5.6/3.r

/var/spool/pbs/mom_priv/jobs/278227.mgt.SC: line 11:  4895 Bus error               (core dumped) Rscript /share/pub/dengcy/GWAS_Multiomics/test/covid19/test5.6/3.r

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
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/mild_all.rds",
                     singlecell=F,
                      output.prefix="mild",
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
Pagwas$Pathway_ld_gwas_data<-NULL
save(Pagwas,file="scpagwas_mild.v1.73.RData")
```

#### 基于1.7.3版本得结果整合1.9.1结果

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
###############子函数
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

整合：

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

#!/usr/bin/sh
#PBS -N moderate
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=5
#PBS -j oe
####mild

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/test/covid19/test5.6/4.r

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
                     Single_data = "/share2/pub/jiangdp/jiangdp/COVID/data/moderate_all.rds",
                     singlecell=F,
                      output.prefix="moderate",
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
Pagwas$Pathway_ld_gwas_data<-NULL
save(Pagwas,file="scpagwas_moderate.v1.73.RData")
```

#### 基于1.7.3版本得结果整合1.9.1结果

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

#!/usr/bin/sh
#PBS -N normal
#PBS -q workq
#PBS -l nodes=node04
#PBS -l ncpus=5
#PBS -j oe
####mild

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/test/covid19/test5.6/5.r


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

#### 基于1.7.3版本得结果整合1.9.1结果

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

## 结果数据降维计算，可视化

#!/usr/bin/sh
#PBS -N severe
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/test/covid19/1.r



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



## 新的验证数据计算结果：

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

