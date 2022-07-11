

# gene_snpscPagwas数据分布比较

1，对gene，gene-snp，pathway-snp画分布图

![image-20220422101632219](E:\OneDrive\GWAS_Multiomics\scriptmd\image-20220422101632219.png)

单细胞处理综述：

方差稳定：**基因表达水平可以有很大的差异，但是基因的平均表达于基因的变异密切相关**，这种效应被称为均值方差关系。

<img src="E:\OneDrive\GWAS_Multiomics\scriptmd\image-20220422095745112.png" width="40%" />

方差稳定可以调整数据，消除基因表达量对基因方差的影响。这一步骤确保下游分析集中于生物学上最相关的基因（即数据集中在特定细胞类型中表达的基因），而不是简单的集中于表达量最高的基因。例如，方差稳定可能通过使具有低平均表达水平的基因（如转录因子），尽管这些基因在细胞内的总体表达水平较低，但它们在揭示细胞命运，**寻找发育祖细胞相关亚群方面**可能很重要。

一种简单的稳定变异的方法是归一化的count数进行对数变换，它可以减少高表达和地表达基因之间的差异。可以用来消除基因平均表达对基因方差影响的流程包括Seurat,Pagoda2和 SCANPY。

总的来说，虽然方差稳定对于单细胞的分析并不是严格必需的，但是**调整数据集以适应平均基因表达的大范围变化可以增强生物学相关基因对下游分析的贡献**。这种方法**消除了在所有细胞中都有丰富的表达但表达水平相似的基因的影响，例如管家基因，这种基因同时对细胞异质性的研究没有帮助**。在方差稳定，之后，识别和选择高度可变的基因可以提高下游分析中细胞类型的分辨率，特别是在被分析的细胞类型相当相似的情况下。这一可选的处理步骤包括在调整平均基因表达差异后选择残差最大的基因。



## 数据处理

```R
library("data.table")
library(Seurat)
#################数据处理
setwd("/share/pub/dengcy/GWAS_Multiomics/Distribution_test")
i<-"RBCcount"
prune_gwas<-read.table(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"),header=T)
gwas_df<-prune_gwas#[sample(1:nrow(prune_gwas),100000),]
save(gwas_df,file="gwas_df.RData")
write.table(gwas_df,file="gwas_df.txt",quote=F,row.names=F)

Single_data=readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_pbmc.rds")
Single_data<-subset(Single_data,idents=c("RBC","HSC"))
saveRDS(Single_data,file="dis_test_scd.rds")

B_cellnames <- rownames(Single_data@meta.data)[Single_data@meta.data$initial_clustering=="B_cell"]
save(B_cellnames,file="B_cellnames.RData")

load("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_data_scexpr.RData")
merge_scexpr=NM_Healthy_data_scexpr[,c("RBC","HSC")]
save(merge_scexpr,file="merge_scexpr.RData")

merge_scexpr<- AverageExpression(Single_data)
save(merge_scexpr,file="merge_scexpr.RData")

Single_data<-readRDS("dis_test_scd.rds")
single_scexpr<- Single_data[['RNA']]@data
single_scexpr<-as.data.frame(data.matrix(single_scexpr))
save(single_scexpr,file="single_scexpr.RData")

load("/share/pub/dengcy/Singlecell/COVID19/1.rolypoly_result/geneid_df1.RData")
 geneid_df1<-geneid_df1[!duplicated(geneid_df1$label),]
save(geneid_df1,file="geneid_df1.RData")
```

## rolypoly

### 1.处理rolypoly的数据分布

```R
library("rolypoly")
setwd("/share/pub/dengcy/GWAS_Multiomics/Distribution_test")


ld_path <- "/share/pub/dengcy/Singlecell/COVID19/data/LD"
load("gwas_df.RData")
load("merge_scexpr.RData")
load("geneid_df1.RData")
##开始计算细胞类型的数据分布
rolypoly_result <- rolypoly_roll(
  gwas_data = gwas_df,
  block_annotation = geneid_df1,
  block_data =merge_scexpr,
  ld_folder =ld_path,
  bootstrap_iters = 100
)

vector_xy<-vectorize_rolypoly(rolypoly_result$data)
celltypes_gene_snp<-vector_xy$x
save(celltypes_gene_snp,file="celltypes_gene_snp.RData")

#!/usr/bin/sh
#PBS -N 1
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
####mild
source activate R3.6
Rscript /share/pub/dengcy/GWAS_Multiomics/Distribution_test/1.r
######单细胞的snp数据分布

library("rolypoly")
setwd("/share/pub/dengcy/GWAS_Multiomics/Distribution_test")
load("geneid_df1.RData")
load("/share/pub/dengcy/Singlecell/COVID19/1.rolypoly_result/geneid_df1.RData")
ld_path <- "/share/pub/dengcy/Singlecell/COVID19/data/LD"
load("gwas_df.RData")

load("single_scexpr.RData")
#single_scexpr<-as.data.frame(data.matrix(single_scexpr))

rolypoly_result <- rolypoly_roll(
  gwas_data = gwas_df,
  block_annotation = geneid_df1,
  block_data = single_scexpr[,1:2],
  ld_folder =ld_path,
  bootstrap_iters = 100
)

vector_xy<-vectorize_rolypoly(rolypoly_result$data)
singlecell_gene_snp<-vector_xy$x
save(singlecell_gene_snp,file="singlecell_gene_snp.RData")

#!/usr/bin/sh
#PBS -N 2
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
####mild
source activate R3.6
Rscript /share/pub/dengcy/GWAS_Multiomics/Distribution_test/2.r
```

画图：

```R

setwd("E:/OneDrive/GWAS_Multiomics/Distribution_test")
load("singlecell_gene_snp.RData")
singlecell_gene_snp<-as.data.frame(singlecell_gene_snp)
celltypes_gene_snp<-as.data.frame(celltypes_gene_snp)
snps<-intersect(rownames(singlecell_gene_snp),rownames(celltypes_gene_snp))
singlecell_gene_snp<-singlecell_gene_snp[snps,]
colnames(singlecell_gene_snp)<-c("cell1","cell2")
celltypes_gene_snp<-celltypes_gene_snp[snps,]
colnames(celltypes_gene_snp)<-c("cell1","cell2")
rolypoly_gg<-data.frame(snp=rownames(singlecell_gene_snp),
                        singlecell=singlecell_gene_snp$cell1,
                        celltype=celltypes_gene_snp$cell1)
library(reshape2)
rolypoly_gg<-reshape2::melt(rolypoly_gg,id.vars="snp")

pdf("snp_gene_distribution.pdf",width = 4,height = 6)
ggplot(rolypoly_gg, aes(x=value,color=variable))+
  geom_density()+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  theme(legend.position="top")
dev.off()


```

## scPagwas的数据分布

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/Distribution_test") 
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)

 Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/RBCcount_prune_gwas_data.txt",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data = "/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_pbmc.rds",
                     nfeatures =NULL,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R3.25/Celltype_heritability_contributions.R")
Pathway_ld_gwas_data <- link_pwpca_block(Pagwas)
celltype_path_snp <- xy2vector(Pathway_ld_gwas_data)
Pagwas$celltype_path_snp <- celltype_path_snp$x
save(Pagwas,file="Pagwas_celltype_path_snp.RData")

#####################单细胞计算
setwd("/share/pub/dengcy/GWAS_Multiomics/Distribution_test") 
library(scPagwas)
load("Pagwas_celltype_path_snp.RData")

#singlecell_path_snp <- xy2vector(Pathway_ld_gwas_data)
#singlecell_path_snp <- singlecell_path_snp$x
#save(singlecell_path_snp,file="singlecell_path_snp.RData")

  paths<-names(Pagwas$Pathway_ld_gwas_data)
  Pathway_x_results <- lapply(Pagwas$Pathway_ld_gwas_data, function(pa_block) {

    pathway <- unique(pa_block$block_info$pathway)
    x <- matrix(Pagwas$pca_scCell_mat[pathway, ],nrow = 1)
    rownames(x)<-pathway

    if (pa_block$n_snps == 0) {
      pa_block$include_in_inference <- F
      pa_block$x <- NULL # to make sure we totally replace previous stuffs
      return(pa_block)
    }

    mg <- intersect(Pagwas$rawPathway_list[[pathway]],rownames(Pagwas$data_mat))
    if (length(mg) == 1) {
      x2<-matrix(Pagwas$data_mat[mg, ],nrow=1)
      x2<- x2/(x2+0.0001)
      rownames(x2)<-mg
    }else{
      x2 <-  biganalytics::apply(Pagwas$data_mat[mg, ],2,function(ge){
          if (sum(ge) == 0) {
            return(rep(0,length(ge)))
          }else{
            return(ge / sum(ge))
          }
      })
     rownames(x2)<-mg
    }
      x2 <- as(x2,"dgCMatrix")

    if (pa_block$n_snps > 1) {
      x2 <-x2[pa_block$snps$label, ]
      pa_block$n_snps <- nrow(pa_block$snps)

      x <- x[rep(1, pa_block$n_snps), ]
      rownames(x) <- pa_block$snps$rsid
      rownames(Pagwas$snp_gene_df) <- Pagwas$snp_gene_df$rsid
      x <- x * Pagwas$snp_gene_df[pa_block$snps$rsid, "slope"]
      x2 <- x2 * x
    } else {
      x2 <- matrix(x2[pa_block$snps$label, ], nrow = 1)
      rownames(x2) <- pa_block$snps$label
      pa_block$n_snps <- nrow(pa_block$snps)

      x <- matrix(x[rep(1, pa_block$n_snps), ], nrow = 1)
      rownames(x) <- pa_block$snps$rsid

      rownames(Pagwas$snp_gene_df) <- Pagwas$snp_gene_df$rsid
      x <- matrix(as.numeric(x) * as.numeric(Pagwas$snp_gene_df[pa_block$snps$rsid, "slope"]), nrow = 1)
      x2 <- matrix(as.numeric(x2) * as.numeric(x), nrow = 1)
      x2 <- as(x2,"dgCMatrix")
    }
    pa_block$x<- as(pa_block$ld_matrix_squared %*% x2,"dgCMatrix")
    pa_block$include_in_inference <- T
    noise_per_snp <- pa_block$snps$se**2

    if (!is.null(pa_block$x)) {
      if (pa_block$n_snps > 2) {

        na_elements <- is.na(pa_block$y) | apply(pa_block$x, 1, function(x) {
          any(is.na(x))
        }) | is.na(noise_per_snp)

        results <- pa_block$x[!na_elements,]
      }else{
        results<-NULL
      }

    } else {
      results<-NULL
    }
    return(results)
  })
save(Pathway_x_results,file="Pathway_x_results.RData")

```



```R
library(scPagwas)
load("Pagwas_celltype_path_snp.RData")

#singlecell_path_snp <- xy2vector(Pathway_ld_gwas_data)
#singlecell_path_snp <- singlecell_path_snp$x
#save(singlecell_path_snp,file="singlecell_path_snp.RData")

  paths<-names(Pagwas$Pathway_ld_gwas_data)
  Pathway_genex_results <- lapply(Pagwas$Pathway_ld_gwas_data, function(pa_block) {

    pathway <- unique(pa_block$block_info$pathway)
    #rownames(x)<-pathway
    if (pa_block$n_snps == 0) {
      pa_block$include_in_inference <- F
      pa_block$x <- NULL # to make sure we totally replace previous stuffs
      return(pa_block)
    }

    mg <- intersect(Pagwas$rawPathway_list[[pathway]],rownames(Pagwas$data_mat))
    if (length(mg) == 1) {
      x2<-matrix(Pagwas$data_mat[mg, ],nrow=1)
      rownames(x2)<-mg
    }else{
      x2 <- Pagwas$data_mat[mg, ]
     rownames(x2)<-mg
    }
      x2 <- as(x2,"dgCMatrix")

    if (pa_block$n_snps > 1) {
      x2 <-x2[pa_block$snps$label, ]
      pa_block$n_snps <- nrow(pa_block$snps)

      } else {
      x2 <- matrix(x2[pa_block$snps$label, ], nrow = 1)
      rownames(x2) <- pa_block$snps$label
      pa_block$n_snps <- nrow(pa_block$snps)
      x2 <- as(x2,"dgCMatrix")
    }
    pa_block$x<- as(pa_block$ld_matrix_squared %*% x2,"dgCMatrix")
    pa_block$include_in_inference <- T
    noise_per_snp <- pa_block$snps$se**2

    if (!is.null(pa_block$x)) {
      if (pa_block$n_snps > 2) {

        na_elements <- is.na(pa_block$y) | apply(pa_block$x, 1, function(x) {
          any(is.na(x))
        }) | is.na(noise_per_snp)

        results <- pa_block$x[!na_elements,]
      }else{
        results<-NULL
      }

    } else {
      results<-NULL
    }
    return(results)
  })

save(Pathway_genex_results,file="Pathway_genex_results.RData")


```

画图

```R
setwd("E:/OneDrive/GWAS_Multiomics/Distribution_test")
#load("singlecell_path_snp.RData")
#load("celltype_path_snp.RData")
singlecell_path_snp<-as.data.frame(singlecell_path_snp)
celltype_path_snp<-as.data.frame(celltype_path_snp)

snps<-intersect(rownames(singlecell_path_snp),rownames(celltype_path_snp))
singlecell_path_snp<-singlecell_path_snp[snps,]
colnames(singlecell_path_snp)<-c("cell1","cell2")
celltype_path_snp<-celltype_path_snp[snps,]
colnames(celltype_path_snp)<-c("cell1","cell2")

pagwas_gg<-data.frame(snp=rownames(singlecell_path_snp),
                        singlecell=singlecell_path_snp$cell1,
                        celltype=celltype_path_snp$cell1)
library(reshape2)
pagwas_gg<-reshape2::melt(pagwas_gg,id.vars="snp")

pdf("Pathway_gene_distribution.pdf",width = 4,height = 6)
ggplot(pagwas_gg, aes(x=value,color=variable))+
  geom_density()+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  theme(legend.position="top")
dev.off()

############################
load("Pathway_genex_results.RData")
load("Pathway_x_results.RData")
#Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),snp_data_list)
x_Pa<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),lapply(Pathway_x_results,function(x){
    as.data.frame(x[,1:100])
}))

x_ge<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2), lapply(Pathway_genex_results,function(x){
    as.data.frame(x[,1:100])
}))
#snp<-intersect(names(x_ge),names(x_Pa))
#gg_df<-data.frame(snp=names(x_ge),gene_snp=x_ge,pathway_snp=x_Pa)
save(x_Pa,file="x_Pa.RData")
save(x_ge,file="x_ge.RData")

load("gene_pathway_snp_singlegg.RData")
library(reshape2)
library(ggplot2)
gg_df<-x_Pa
gg_df$snp<-rownames(gg_df)
pagwas_gg_df1<-reshape2::melt(gg_df,id.vars="snp")

pdf("Pathway_snp_singledistribution.pdf",width =5,height = 5)
ggplot(pagwas_gg_df1, aes(x=value,color=variable))+
  geom_density()+scale_y_continuous(limits = c(0, 60))+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  labs(title="SNP-gene-Pathway")+
  theme(legend.position="none")
dev.off()

gg_df<-x_ge
gg_df$snp<-rownames(gg_df)
pagwas_gg_df2<-reshape2::melt(gg_df,id.vars="snp")

pdf("Gene_snp_singledistribution.pdf",width =5,height = 5)
ggplot(pagwas_gg_df2, aes(x=value,color=variable))+
  geom_density()+scale_y_continuous(limits = c(0, 6))+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  labs(title="SNP-Gene")+
  theme(legend.position="none")
dev.off()

```

## 方差分布

选择B细胞亚群作为方差计算的数据

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/Distribution_test") 
load("Pathway_genex_results.RData")
load("Pathway_x_results.RData")
load("B_cellnames.RData")
#Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),snp_data_list)
B_x_Pa<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),lapply(Pathway_x_results,function(x){
    x<-as.data.frame(data.matrix(x))
    x<-x[,B_cellnames]
}))

B_x_ge<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2), lapply(Pathway_genex_results,function(x){
     x<-as.data.frame(data.matrix(x))
    x<-x[,B_cellnames]
}))

save(B_x_Pa,file="B_x_Pa.RData")
save(B_x_ge,file="B_x_ge.RData")
var_mean_df<-data.frame(
    snp=rownames(B_x_Pa),
 var_x_Pa=unlist(apply(B_x_Pa,1,var)),
 mean_x_Pa=unlist(apply(B_x_Pa,1,mean)),
 var_x_ge=unlist(apply(B_x_ge,1,var)),
 mean_x_ge=unlist(apply(B_x_ge,1,mean)))
 save(var_mean_df,file="var_mean_df.RData")
```

画图：

```R
setwd("E:/OneDrive/GWAS_Multiomics/Distribution_test")
load("var_mean_df.RData")
var_mean_df$mean_x_Pa<-  log10(abs(var_mean_df$mean_x_Pa))
var_mean_df$mean_x_ge<-  log10(abs(var_mean_df$mean_x_ge))
# blues9
var_mean_df$mean_x_ge[which(var_mean_df$mean_x_ge<=-5)]<- -5
var_mean_df$mean_x_Pa[which(var_mean_df$mean_x_Pa<=-6)]<- -6

pdf("smoothScatter_pathway_snps.pdf")
smoothScatter(var_mean_df$mean_x_Pa,var_mean_df$var_x_Pa, 
              colramp = colorRampPalette(c("white", "#B91646")),
              transformation = function(x) x^.3,
              col = "black",
              ylim=c(0,3),
              xlab = "snp-pathway log10(x)", 
              ylab = "Variance")
lines(lowess(var_mean_df$mean_x_Pa, var_mean_df$var_x_Pa),col='#105652',lwd=1,lty=2)

dev.off()

pdf("smoothScatter_genes_snp.pdf")
smoothScatter(var_mean_df$mean_x_ge,var_mean_df$var_x_ge,
              colramp = colorRampPalette(c("white", "#B91646")),
              col = "black",
              transformation = function(x) x^.3,
              ylim=c(0,3),
              xlab = "snp-genes log10(x)", 
              ylab = "Variance")
lines(lowess(var_mean_df$mean_x_ge, var_mean_df$var_x_ge),col='#105652',lwd=1.5,lty=2)
dev.off()
```



## 通路PCA和单细胞表达分布

```R
#load("E:/OneDrive/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_Prune_adsubset_kegg.RData")
load("E:/OneDrive/GWAS_Multiomics/compare/Hema_test2/Lymphocytecount3_Hema_bmmc_scPagwas_v1.9.1.RData")
pca_scCell_mat <- GetAssayData(Pagwas,assay = "scPagwasPaPca")
median<- apply(pca_scCell_mat,1,median)
#gg_pca<-as.data.frame(pca_scCell_mat[,sample(1:ncol(pca_scCell_mat),100)])
gg_pca<-data.frame(median=median,pathway=rownames(pca_scCell_mat))
#gg_pca$pathway<-rownames(gg_pca)
#gg_pca$median<-median



library(reshape2)
#gg_pca_mel<-reshape2::melt(gg_pca,id.vars="pathway")

pdf("E:/OneDrive/GWAS_Multiomics/Distribution_test/Pathway_gene_distribution.pdf",width = 4,height = 6)
ggplot(gg_pca, aes(x=median))+
  geom_density()+
  #scale_color_brewer()+
  theme_classic()+
labs(title="pathway")+
 theme(legend.position="none")
dev.off()

gg_scgene<-as.data.frame(Pagwas$data_mat[,1:10])
gg_scgene$gene<-rownames(gg_scgene)
#library(reshape2)
gg_scgene<-reshape2::melt(gg_scgene,id.vars="gene")

pdf("E:/OneDrive/GWAS_Multiomics/Distribution_test/singlegene_distribution.pdf",width = 4,height = 6)
ggplot(gg_scgene, aes(x=value,color=variable))+
  geom_density()+
  scale_color_brewer()+
  theme_classic()+
labs(title="gene")+
 theme(legend.position="none")
dev.off()
```

## 基因和通路得方差分布

```R
library("ggplot2")
library(ggpubr)
library(scPagwas)
suppressMessages(library(Seurat))

setwd("/share/pub/dengcy/GWAS_Multiomics/Distribution_test")
#load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas_v1.9.1.RData")
#load("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/monocytecount_pbmc_scPagwas_singlecell.RData")
#load("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_kegg.RData")
#load("/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/scPgwas/scpagwas_brain_Obesity.RData")
load("/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/scPgwas/scpagwas_HCA_Migraine.RData")
#data_mat <- GetAssayData(Pagwas,assay = "RNA")
#data_mat <- as_matrix(data_mat)
#data_mat2<-GetAssayData(Pagwas,assay = "scPagwasPaPca")
Pagwas$Pathway_ld_gwas_data<-NULL
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

###另一种情况
data_mat <- Pagwas$data_mat
#data_mat <- as_matrix(data_mat)
data_mat2 <-Pagwas$pca_scCell_mat
paths<-rownames(data_mat2)
#############需要构造一个pathway得分得矩阵，该矩阵得行是基因，并且是标化之后得值

x_list<-lapply(paths,function(i){
    gene<- Genes_by_pathway_kegg[[i]]
    x2<-data_mat[intersect(gene,rownames(data_mat)),]
    x2<- as_matrix(x2)
     x2 <- apply(x2, 2, function(x) (x - min(x)) / (max(x) - min(x)))
    return(x2)
})
data_mat2<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),x_list)

x_list<-lapply(paths,function(i){
    gene<- Genes_by_pathway_kegg[[i]]
    x2<-data_mat[intersect(gene,rownames(data_mat)),]
    return(x2)
})
data_mat<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),x_list)
              
pa_var_mean_df<-data.frame(
 pathway=rownames(data_mat2),
 var_x_Pa=unlist(apply(data_mat2,1,var)),
 mean_x_Pa=unlist(apply(data_mat2,1,mean)))

 gene_sparsity<-apply(data.matrix(1:nrow(data_mat)),1,function(x){
  a<-data_mat[x,]
  return(sum(a==0)/length(a))
})
 gene_var_mean_df<-data.frame(  
     gene=rownames(data_mat),
     var_x_ge=unlist(apply(data.matrix(1:nrow(data_mat)),1,function(x){
  a<-data_mat[x,]
  return(var(a))
})),
 mean_x_ge=unlist(apply(data.matrix(1:nrow(data_mat)),1,function(x){
  a<-data_mat[x,]
  return(mean(a))
}))


#setwd("E:/OneDrive/GWAS_Multiomics/Distribution_test")
#load("var_mean_df.RData")
pa_var_mean_df$mean_x_Pa<-  log10(abs(pa_var_mean_df$mean_x_Pa))
gene_var_mean_df$mean_x_ge<-  log10(abs(gene_var_mean_df$mean_x_ge))
# blues9
#gene_var_mean_df$mean_x_ge[which(var_mean_df$mean_x_ge<=-5)]<- -5
#pa_var_mean_df$mean_x_Pa[which(var_mean_df$mean_x_Pa<=-6)]<- -6

pdf("HCA_smoothScatter_pathway.pdf")
smoothScatter(pa_var_mean_df$mean_x_Pa,pa_var_mean_df$var_x_Pa, 
              colramp = colorRampPalette(c("white", "#B20600")),
              transformation =  function(x) x^.2,
              col = "black",
              ylim=c(0,3),
              xlab = "snp-pathway log10(x)", 
              ylab = "Variance")
lines(lowess(pa_var_mean_df$mean_x_Pa,pa_var_mean_df$var_x_Pa),col='#B20600',lwd=1,lty=2)
dev.off()

pdf("HCA_smoothScatter_genes.pdf")
smoothScatter(gene_var_mean_df$mean_x_ge,gene_var_mean_df$var_x_ge,
              colramp = colorRampPalette(c("white", "#398AB9")),
              col = "black",
              transformation =  function(x) x^.4,
              ylim=c(0,3),
              xlab = "snp-genes log10(x)", 
              ylab = "Variance")
lines(lowess(gene_var_mean_df$mean_x_ge, gene_var_mean_df$var_x_ge),col='#243A73',lwd=1.5,lty=1)
dev.off()
```

# 单细胞数据与通路数据，零值0膨胀比例图

```R
library("ggplot2")
library(ggpubr)
library(scPagwas)
suppressMessages(library(Seurat))

setwd("/share/pub/dengcy/GWAS_Multiomics/Distribution_test")
#load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas_v1.9.1.RData")

#load("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/monocytecount_pbmc_scPagwas_singlecell.RData")
#load("/share/pub/dengcy/GWAS_Multiomics/ad_test/5.6test/Pagwas_GSE138852_kegg.RData")
#load("/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/scPgwas/scpagwas_brain_Obesity.RData")
load("/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/scPgwas/scpagwas_HCA_Migraine.RData")
Pagwas$Pathway_ld_gwas_data<-NULL
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

#data_mat <- GetAssayData(Pagwas,assay = "RNA")
#data_mat <- as_matrix(data_mat)

data_mat <- Pagwas$data_mat
data_mat2 <-Pagwas$pca_scCell_mat

rm(Pagwas)
#data_mat <- as_matrix(data_mat)
gene_sparsity<-apply(data.matrix(1:nrow(data_mat)),1,function(x){
  a<-data_mat[x,]
  return(sum(a==0)/length(a))
})

cellsingle_sparsity<-apply(data.matrix(1:ncol(data_mat),nrow=1),2,function(x){
  a<-data_mat[,x]  
  return(sum(a==0)/length(a))
})
gene_sparsity2<-round(gene_sparsity,3)
cellsingle_sparsity2<-round(cellsingle_sparsity,3)
#rm(data_mat)

# data_mat2<-GetAssayData(Pagwas,assay = "scPagwasPaPca")
Pathway_sparsity<-apply(data_mat2,1,function(x){
  x<-round(x,2)
  return(sum(x==0)/length(x))
})

cellPathway_sparsity<-apply(data_mat2,2,function(x){
  x<-round(x,2)
  return(sum(x==0)/length(x))
})
Pathway_sparsity2<-round(Pathway_sparsity,4)
cellPathway_sparsity2<-round(cellPathway_sparsity,4)

df<- data.frame(gene_sparsity2)
#pdf("gene_sparsity_distribution.pdf",width = 4,height = 4)
p1<-ggplot(df, aes(x=gene_sparsity2))+
  geom_density(fill="#3A5BA0",color="#3A5BA0",alpha=0.8)+
  scale_color_brewer()+
  theme_classic()+
  labs(title="",x="Sparsity of gene(%)")+
  theme(legend.position="none")
#dev.off()


df2<- data.frame(Pathway_sparsity2)

#pdf("Pathway_sparsity_distribution.pdf",width = 4,height = 4)
p2<-ggplot(df2, aes(x=Pathway_sparsity2))+
  geom_density(fill="#FF5B00",color="#FF5B00",alpha=0.8)+
  scale_color_brewer()+
  theme_classic()+
  labs(title="",x="Sparsity of pathway(%)")+
  theme(legend.position="none")
#dev.off()


df3<- data.frame(cellsingle_sparsity2)
#pdf("cellsingle_sparsity_distribution.pdf",width = 4,height = 4)
p3<-ggplot(df3, aes(x=cellsingle_sparsity2))+
  geom_density(fill="#3A5BA0",color="#3A5BA0",alpha=0.8)+
  scale_color_brewer()+
  theme_classic()+
  labs(title="",x="Sparsity of cell in gene profile(%)")+
  theme(legend.position="none")
#dev.off()


df4<- data.frame(cellPathway_sparsity2)

#pdf("cellPathway_sparsity_distribution.pdf",width = 4,height =4)
p4<-ggplot(df4, aes(x=cellPathway_sparsity2))+
  geom_density(fill="#FF5B00",color="#FF5B00",alpha=0.8)+
  scale_color_brewer()+
  theme_classic()+
  labs(title="",x="Sparsity of cell in pathway profile(%)")+
  theme(legend.position="none")

  median<- apply(data_mat2,1,median)
#gg_pca<-as.data.frame(pca_scCell_mat[,sample(1:ncol(pca_scCell_mat),100)])
gg_pca<-data.frame(median=median,pathway=rownames(data_mat2))
#gg_pca$pathway<-rownames(gg_pca)

#pdf("Pathway_gene_distribution.pdf",width = 4,height = 6)
p5<-ggplot(gg_pca, aes(x=median))+
  geom_density(fill="#FF5B00",color="#FF5B00",alpha=0.8)+
  theme_classic()+
labs(title="pathway",x="Median expression")+
 theme(legend.position="none")
#dev.off()

median2<- apply(data_mat,1,median)
gg_scgene<-data.frame(median=median2,gene=rownames(data_mat))
#library(reshape2)
#gg_scgene<-reshape2::melt(gg_scgene,id.vars="gene")

#pdf("singlegene_distribution.pdf",width = 4,height = 6)
p6<-ggplot(gg_scgene, aes(x=median))+
  geom_density(fill="#3A5BA0",color="#3A5BA0",alpha=0.8)+
  theme_classic()+
labs(title="gene",x="Median expression")+
 theme(legend.position="none")
#dev.off()
#dev.off()
pdf("HCA_cellPathway_sparsity_distribution.pdf",width =7,height =5)
patchwork::wrap_plots(plotlist = list(p1,p3,p6,p2,p4,p5), ncol = 3)
dev.off()

```

