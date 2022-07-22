 血细胞trait进行scPagwas的计算

[toc]

<html>
<!--在这里插入内容-->
</html>
GWAS summary statistics and fine-mapping analysis 
Blood cell traits. 
Summary statistics of 22 blood cell traits from Blood Cell Consortium 2 (BCX2) were processed 
as previously described

> (Vuckovic, D. et al. The Polygenic and Monogenic Basis of Blood Traits and Diseases. Cell 182, 1214–1231.e11 (2020)). 

### blood traits：
Data and Code Availability Summary statistics are available to download from: ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/UKBB_blood_cell_traits/ for UK Biobank and http://www.mhi-humangenetics.org/en/resources for the meta-analysis. The accession numbers for the UK Biobank summary statistics reported in this paper are: 

- basophil cell count
  Dataset: ieu-b-29
- eosinophil cell count
  Dataset: ieu-b-33
- lymphocyte cell count
  Dataset: ebi-a-GCST004627
- monocyte cell count
  Dataset: ieu-b-31
- neutrophil cell count
  Dataset: ieu-b-34
- white blood cell count
  Dataset: ieu-b-30
- Lymphocyte percentage of white cells
  Dataset: ebi-a-GCST004632
- Hemoglobin concentration
  Dataset: ebi-a-GCST90002311
- Mean corpuscular hemoglobin concentration
  Dataset: ebi-a-GCST90002329
- Mean corpuscular volume
  Dataset: ebi-a-GCST90002335



## PBMC and BMMC scRNA-seq data

### 1.prepare the singlecell result

#### PBMC single cell data

```R

load("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/phe_NM_Healthy_pbmc.Rdata")
NM_Healthy_pbmc
An object of class Seurat 
24737 features across 97039 samples within 1 assay 
Active assay: RNA (24737 features, 2000 variable features)
3 dimensional reductions calculated: protein_expression, pca, umap
table(NM_Healthy_pbmc$full_clustering)
                 ASDC           B_exhausted            B_immature 
                   16                   245                   650 
              B_naive B_non-switched_memory     B_switched_memory 
                 4887                   499                  1049 
         C1_CD16_mono             CD14_mono             CD16_mono 
                  198                  3774                  3007 
          HSC_CD38neg           HSC_CD38pos                CD4.CM 
                   11                    50                  7437 
               CD4.EM              CD4.IL22             CD4.Naive 
                  183                  6298                 12952 
           CD4.Prolif               CD4.Tfh               CD4.Th1 
                   29                   772                   101 
             CD4.Th17               CD4.Th2                CD8.EM 
                    1                    16                  5410 
            CD8.Naive            CD8.Prolif                CD8.TE 
                 7432                    67                  6145 
       CD83_CD14_mono                   DC1                   DC2 
                 6772                   104                  1052 
                  DC3             DC_prolif         HSC_erythroid 
                  975                     6                    67 
               ILC1_3                  ILC2                  MAIT 
                  213                    36                  3361 
          Mono_prolif           HSC_myeloid                   NKT 
                    4                     3                   823 
              NK_16hi               NK_56hi             NK_prolif 
                12375                  1959                   205 
      Plasma_cell_IgA       Plasma_cell_IgG       Plasma_cell_IgM 
                  121                   111                    54 
          Plasmablast             Platelets            HSC_prolif 
                   64                  1722                     6 
                  RBC                  Treg                   gdT 
                  291                    75                  4733 
                  pDC 
                  678 

table(NM_Healthy_pbmc$initial_clustering)
     B_cell          CD4          CD8         CD14         CD16          DCs 
        7568        28268        18796        10312         3466         2115 
         HSC Lymph_prolif         MAIT  Mono_prolif      NK_16hi      NK_56hi 
         144          282         2014            4        12640         2222 
 Plasmablast    Platelets          RBC         Treg          gdT          pDC 
         115         2145          327         2340         3592          689 

Idents(NM_Healthy_pbmc)<-NM_Healthy_pbmc$initial_clustering
saveRDS(NM_Healthy_pbmc,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_pbmc.rds")
```

#### BMMC single cell data

```R
library("scPagwas")
library("Seurat")
library("scRNAseq")
library("SingleCellExperiment")
library("stringr") 
scRNA_Healthy_Hema<-readRDS("E:/OneDrive/SingleCell/data/PBMCscATAC-seq/scRNA-Healthy-Hematopoiesis-191120.rds")
counts <- assay(scRNA_Healthy_Hema, "counts")
Seu_Healthy_Hema <- CreateSeuratObject(
               counts = counts, 
               meta.data=as.data.frame(colData(scRNA_Healthy_Hema)),
               min.cells = 3, 
               min.features = 200)

Idents(Seu_Healthy_Hema)<-scRNA_Healthy_Hema@colData$BioClassification
table(Idents(Seu_Healthy_Hema))

        01_HSC 02_Early.Eryth  03_Late.Eryth  04_Early.Baso    05_CMP.LMPP       06_CLP.1 
          1425           1653            446            111           2260            903 
        07_GMP    08_GMP.Neut         09_pDC         10_cDC 11_CD14.Mono.1 12_CD14.Mono.2 
          2097           1050            544            325           1800           4222 
  13_CD16.Mono         14_Unk       15_CLP.2       16_Pre.B           17_B      18_Plasma 
           292            520            377            710           1711             62 
      19_CD8.N      20_CD4.N1      21_CD4.N2       22_CD4.M      23_CD8.EM      24_CD8.CM 
          1521           2470           2364           3539            796           2080 
         25_NK         26_Unk 
          2143            161 

Seu_Healthy_Hema <- ScaleData(Seu_Healthy_Hema)
Seu_Healthy_Hema <- NormalizeData(Seu_Healthy_Hema, normalization.method = "LogNormalize", scale.factor = 10000)

```

#### Run the single_data_input

```R
####PBMC
setwd("/share/pub/dengcy/GWAS_Multiomics/compare")
 library(scPagwas)
 suppressMessages(library(Seurat))
Single_data=readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_pbmc.rds")
    Pagwas <- Single_data_input(Pagwas=NULL,
                                assay="RNA",
                                Single_data=Single_data,
                                Pathway_list=Genes_by_pathway_kegg)
  Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,
                                 Pathway_list=Genes_by_pathway_kegg
                                 )
save(Pagwas ,file="NM_Healthy_pbmc_prePagwas.RData")

##BMMC
library(scPagwas)
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")
  Pagwas <- Single_data_input(Pagwas=NULL,
                                assay="RNA",
                                Single_data=Single_data,
                                Pathway_list=Genes_by_pathway_kegg)
  Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,
                                 Pathway_list=Genes_by_pathway_kegg
                                 )
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Seu_Healthy_Hema_kegg_prePagwas.RData")
```

### 2.Run scPagwas

export OPENBLAS_NUM_THREADS=1

#### BMMC

```R
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
load("/share/pub/dengcy/GWAS_Multiomics/compare/Seu_Healthy_Hema_kegg_prePagwas.RData")

traits<-c("basophilcount","eosinophilcount" ,"Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")
for(i in traits){
Pagwas<-scPagwas_main(Pagwas =Pagwas,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"prune_gwas_data.txt"),
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds",
                     output.prefix=i,
                     Pathway_list=Genes_by_pathway_kegg,
                     output.dirs=paste0(i,"_scPagwas"),
                     ncores=2,
                      assay="RNA",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
 }
```

#### PBMC

```R
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))

suppressMessages(library(Seurat))
load("/share/pub/dengcy/GWAS_Multiomics/compare/NM_Healthy_pbmc_prePagwas.RData")
#################
#celltypes
#################
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/")
for(i in traits){
     Pagwas<-scPagwas_main(Pagwas =Pagwas,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_gwas_data.txt"),
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_pbmc.rds",
                     output.prefix=i,
                           singlecell=F,
                     Pathway_list=Genes_by_pathway_kegg,
                     output.dirs=paste0(i,"_pbmc_scPagwasv1.9"),
                     ncores=5,
                     assay="RNA",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/",i,"_pbmc_scPagwas.RData"))
    }

#################
#single cell
##################
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc")
suppressMessages(library(Seurat))
load("/share/pub/dengcy/GWAS_Multiomics/compare/NM_Healthy_pbmc_prePagwas.RData")

traits2<-c("basophilcount",,"Lymphocytecount3","monocytecount","Hemoglobinconcen","MeanCorpusVolume")

for(i in traits2){
     Pagwas<-scPagwas_main(Pagwas =Pagwas,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"),
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_pbmc.rds",
                     output.prefix=i,
                     celltype=F,
                     Pathway_list=Genes_by_pathway_kegg,
                     output.dirs=paste0(i,"_pbmc_scPagwasv1.9"),
                     ncores=1,
                     assay="RNA",
                     split_n=4,
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/",i,"_pbmc_scPagwas_singlecell.RData"))
}
```

### 4.Integrate bmmc result and visualize

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
#"RBCcount","Plateletcount",
traits<-c("basophilcount","eosinophilcount" ,"Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")
i<-traits[1]
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))

result_list<-lapply(traits,function(i){
   print(i)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
    a<-Pagwas@misc$bootstrap_results[-1,]
    return(a$bp_value)
})
names(result_list)<-traits
result_list<-as.data.frame(result_list)
#load("RBCcount_Hema_bmmc_scPagwas_v1.7.RData")
rownames(result_list)<-rownames(Pagwas@misc$bootstrap_results[-1,])
save(result_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/scPagwas_bmmc_list.RData")
write.csv(result_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/scPagwas_bmmc_list.csv")

setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc")

result_list<-lapply(traits,function(i){
   print(i)
    load(paste0(i,"_pbmc_scPagwas.RData"))
        a<-Pagwas$bootstrap_results[-1,]
        a<-a$bp_value
    return(a)
})

names(result_list)<-traits
result_list<-as.data.frame(result_list)
load("monocytecount_pbmc_scPagwas.RData")
rownames(result_list)<-rownames(Pagwas$bootstrap_results[-1,])

result_list<-result_list[,1:12]
save(result_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scPagwas_PBMC_celltypes_resultlist.RData")
write.csv(result_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scPagwas_pbmc_list.csv")
```

## Blood traits run MAGMA

### 1.Preprogress

```R
###
library(bigreadr)
traits<-c("basophilcount","eosinophilcount" ,"Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

for(i in traits){
gwas<-bigreadr::fread2(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"))
print(head(gwas))
magma_Input1<-gwas[,c("rsid","p","N")]
magma_Input2<-gwas[,c("rsid","chrom","pos","pos")]
write.table(magma_Input2,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_",i,"Input2.txt"),sep="\t",row.names=F,quote=F,col.names=F)
colnames(magma_Input1)<-c("SNP","P","N")
write.table(magma_Input1,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_",i,"Input1.txt"),sep="\t",row.names=F,quote=F)
}

```

### 2.SNP annotation

```sh
cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
for i in basophilcount  
#basophilcount eosinophilcount lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount LymphocytePercent Hemoglobinconcen MeanCorpuscularHemoglobin MeanCorpusVolume
do
./magma --annotate window=10,10 --snp-loc /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_${i}Input2.txt \
--gene-loc /share/pub/dengcy/Singlecell/COVID19/MAGMA/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded \
--out /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10_down
done

```

### 3. gene-level association

```sh
cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
#basophilcount 
for i in basophilcount 
do
./magma --bfile /share/pub/dengcy/Singlecell/COVID19/MAGMA/g1000_eur/g1000_eur \
--pval /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_${i}Input1.txt ncol=3 \
--gene-annot /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10_down.genes.annot \
--out /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10down
done

cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
for i in neutrophilcount WhiteBloodCellcount LymphocytePercent
do
./magma --bfile /share/pub/dengcy/Singlecell/COVID19/MAGMA/g1000_eur/g1000_eur \
--pval /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_${i}Input1.txt ncol=3 \
--gene-annot /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10_down.genes.annot \
--out /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10down
done

cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
for i in Hemoglobinconcen MeanCorpuscularHemoglobin MeanCorpusVolume
do
./magma --bfile /share/pub/dengcy/Singlecell/COVID19/MAGMA/g1000_eur/g1000_eur \
--pval /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_${i}Input1.txt ncol=3 \
--gene-annot /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10_down.genes.annot \
--out /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10down
done
```

### 4.Single data compute

```R
library(tidyverse)
library("rhdf5")
library("snow")
load("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data_scexpr.RData")
load("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_data_scexpr.RData")

cell_expr<-lapply(list(Hema=Seu_Hema_data_scexpr,pbmc=NM_Healthy_data_scexpr),function(expr){
#colnames(expr)<-annotation$V2
expr<-expr[apply(expr,1,sum)!=0,]
return(expr)
})

##########################
###1.gene cord
############################
#Filtered to remove extended MHC (chr6, 25Mb to 34Mb).
gene_coordinates <-read_tsv("/share/pub/dengcy/Singlecell/COVID19/MAGMA/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded",
col_names = FALSE,col_types = 'cciicc') %>%
  mutate(start=ifelse(X3-50000<0,0,X3-50000),end=X4+50000) %>%
  select(X2,start,end,6,1) %>%
  rename(chr="X2", Gene="X6",EntrezGene="X1") %>%
  mutate(chr=paste0("chr",chr))

#########top 10 genes
exp2<-lapply(exp,function(x){
  x[which(x==Inf)]<- max(x[which(x!=Inf)])  
    return(x)
})
exp2<-as.data.frame(exp2)

top10_function <-function(exp=mouse_brain_scexpr){
                exp$Gene <- rownames(exp)
                #Only keep genes with a unique name and tidy data.
                exp <- exp %>% add_count(Gene) %>%
                filter(n==1) %>%
                select(-n) %>%
                gather(key = column,value=Expr,-Gene) %>%
                as.tibble()
###2.Each cell type is scaled to the same total number of molecules.
                exp <- exp %>%
                group_by(column) %>%
                #mutate(Expr_sum_mean=Expr*1e6/sum(Expr))
                mutate(Expr_sum_mean=Expr/sum(Expr))
                exp<- exp %>%
                group_by(Gene) %>%
                mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean)) %>%
                ungroup()
                ###4.Get MAGMA genes
                exp2 <- inner_join(exp,gene_coordinates,by="Gene")
                n_genes <- length(unique(exp2$EntrezGene))
                n_genes_to_keep <- (n_genes * 0.1) %>% round()
                ###7.Write MAGMA/LDSC input files
                exp2 %>% filter(Expr_sum_mean>1) %>% magma_top10("column")
                #exp3 %>% filter(Expr_sum_mean>1) %>% ldsc_bedfile("column")
                print("sucess!")
                }

##########################
###6.Functions
#####################
##Get MAGMA input top10%
magma_top10 <- function(d,Cell_type,n_genes_to_keep=10){
  d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity)
  d_spe %>% do(write_group_magma(.,Cell_type))
}

write_group_magma  = function(df,Cell_type) {
  df <- select(df,column,EntrezGene)
  df_name <- make.names(unique(df[1]))
  colnames(df)[2] <- df_name  
  dir.create(paste0("MAGMA/"), showWarnings = FALSE)
  select(df,2) %>% t() %>% as.data.frame() %>% rownames_to_column("Cat") %>%
    write_tsv("/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/top10.txt",append=T)
  return(df)
}

#####top10_function
top10_function(exp=as.data.frame(cell_expr$Hema))
top10_function(as.data.frame(cell_expr$pbmc))
```

### 5.MAGMA main function

```R
#1.
cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
for i in basophilcount  
#basophilcount eosinophilcount Lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount LymphocytePercent Hemoglobinconcen MeanCorpuscularHemoglobin MeanCorpusVolume
do
int="/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10down.genes.raw"
cell_type="/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/pbmc_top10.txt"
./magma --gene-results  $int --set-annot  $cell_type --out /share/pub/dengcy/GWAS_Multiomics/compare/magma/${i}_magma_pbmc
done

cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
for i in lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount LymphocytePercent Hemoglobinconcen MeanCorpuscularHemoglobin MeanCorpusVolume
do
int="/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10down.genes.raw"
cell_type="/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/pbmc_top10.txt"
./magma --gene-results  $int --set-annot  $cell_type --out /share/pub/dengcy/GWAS_Multiomics/compare/magma/${i}_magma_pbmc
done

```

### 6.Integrate magma result

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/magma")
ds<-c("basophilcount","eosinophilcount" ,"Lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")
magma_result_list<-lapply(ds,function(i){
    message(i)
   result<-read.table(file=paste0(i,"_magma_Hema.gsa.out"),header = T)
    return(result)
  #cell_mild$VARIABLE<-annotation_cell$V2  
})
names(magma_result_list)<-ds
save(magma_result_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/magma_Hema_result_list.RData")

magma_result_list<-lapply(ds,function(i){
    message(i)
   result<-read.table(file=paste0(i,"_magma_pbmc.gsa.out"),header = T)
    return(result)
  #cell_mild$VARIABLE<-annotation_cell$V2  
})
names(magma_result_list)<-ds
save(magma_result_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/magma_pbmc_result_list.RData")

```

## Visualize the result for different methods

### 1.scPagwas

```R
library(gplots)
library(RColorBrewer)

setwd("E:/OneDrive/GWAS_Multiomics/Compare")
load("scPagwas_bmmc_list.RData")
result_list2<- -log2(result_list)
#magma_pbmc_p

coul <- colorRampPalette(brewer.pal(8, "Oranges"))(25)
p2star <- function(p){
  symnum(p,cutpoints = c(0,0.001,0.01,0.05,1),
         symbols = c('***','**','*',''),na = NA)
}
hM <- apply(result_list,2, function(x) as.character(p2star(x)))
rownames(hM)<-rownames(result_list2)

result_list2<-as.matrix(result_list2)
result_list2[which(result_list2>20)]<-20

pdf("E:/OneDrive/GWAS_Multiomics/Compare/Figure_scPagwas_bmmc_cellytpe_p.pdf",
    height =8,width =6)
heatmap.2(result_list2,
          trace="none",#
          col=coul,#
          density.info = "none",#
          key.xlab ='Correlation',
          key.title = "",
          cexRow = 1,cexCol = 1,#
          Rowv = F,Colv =T, #
          margins = c(6, 6),
          cellnote = hM,notecol='black'
)
dev.off()
```

### 2.magma

```R
##########################magma
load("magma_bmmc_list.RData")
result_list<-magma_bmmc_list[,3:12]
result_list2<- -log2(result_list)
#magma_pbmc_p

coul <- colorRampPalette(brewer.pal(8, "Oranges"))(25)#
#hM <- format(round(result_list2, 2))#
p2star <- function(p){
  symnum(p,cutpoints = c(0,0.001,0.01,0.05,1),
         symbols = c('***','**','*',''),na = NA)
}
hM <- apply(result_list,2, function(x) as.character(p2star(x)))
rownames(hM)<-rownames(result_list2)


result_list2<-as.matrix(result_list2)
result_list2[which(result_list2>20)]<-20

pdf("E:/OneDrive/GWAS_Multiomics/Compare/Figure_magma_bmmc_cellytpe_p.pdf",
    height =8,width =6)
heatmap.2(result_list2,
          trace="none",#
          col=coul,#
          density.info = "none",
          key.xlab ='Correlation',
          key.title = "",
          cexRow = 1,cexCol = 1,
          Rowv = F,Colv =T, #
          margins = c(6, 6),
          cellnote = hM,notecol='black'#
)
dev.off()
```



### 3.rolypoly

```R
setwd("/share2/pub/jiangdp/jiangdp/COVID/rolypoly")

traits<-c("basophilcount","eosinophilcount" ,"Lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

rolypoly_pbmc_list<-lapply(traits,function(i){
    message(i)
    load(paste0("roly_NM_",i,".RData"))
    return(rolypoly_result$bootstrap_results$bp_value[-1])
})
names(rolypoly_pbmc_list)<-traits
rolypoly_pbmc_list<-as.data.frame(rolypoly_pbmc_list)
load("roly_NM_RBCcount.RData")
rownames(rolypoly_pbmc_list)<-rownames(rolypoly_result$bootstrap_results[-1,])

save(rolypoly_pbmc_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/rolypoly_pbmc_list.RData")
write.csv(rolypoly_pbmc_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/rolypoly_pbmc_list.csv")

#############BMMC
setwd("/share2/pub/jiangdp/jiangdp/COVID/rolypoly")
rolypoly_bmmc_list<-lapply(traits,function(i){
    message(i)
    load(paste0("roly_Seu_",i,".RData"))
    return(rolypoly_result$bootstrap_results$bp_value[-1])
})
names(rolypoly_bmmc_list)<-traits
rolypoly_bmmc_list<-as.data.frame(rolypoly_bmmc_list)
load("roly_Seu_RBCcount.RData")
rownames(rolypoly_bmmc_list)<-rownames(rolypoly_result$bootstrap_results[-1,])

save(rolypoly_bmmc_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/rolypoly_BMMC_list.RData")
write.csv(rolypoly_bmmc_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/rolypoly_BMMC_list.csv")
#################rolypoly
load("rolypoly_bmmc_list.RData")
result_list<-rolypoly_bmmc_list[,3:12]
result_list2<- -log2(result_list)
#magma_pbmc_p

coul <- colorRampPalette(brewer.pal(8, "Oranges"))(25)#
#hM <- format(round(result_list2, 2))#
p2star <- function(p){
  symnum(p,cutpoints = c(0,0.001,0.01,0.05,1),
         symbols = c('***','**','*',''),na = NA)
}
hM <- apply(result_list,2, function(x) as.character(p2star(x)))
rownames(hM)<-rownames(result_list2)

result_list2<-as.matrix(result_list2)
result_list2[which(result_list2>20)]<-20

pdf("E:/OneDrive/GWAS_Multiomics/Compare/Figure_rolypoly_bmmc_cellytpe_p.pdf",
    height =8,width =6)
heatmap.2(result_list2,
          trace="none",#
          col=coul,#
          density.info = "none",#
          key.xlab ='Correlation',
          key.title = "",
          cexRow = 1,cexCol = 1,#
          Rowv = F,Colv =T, #
          margins = c(6, 6),
          cellnote = hM,notecol='black'#
)
dev.off()
```

### 4.ldsc

```R
setwd("/share2/pub/zhenggw/zhenggw/scPagwas_LDSC/NM_results")
traits<-c("basophilcount","eosinophilcount" ,"Lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

ldsc_pbmc_result<-lapply(traits,function(i){
    message(i)
   result<-read.table(file=paste0(i,".cell_type_results.txt"),header = T)
    return(result)
})
names(ldsc_pbmc_result)<-traits

save(ldsc_pbmc_result,file="/share/pub/dengcy/GWAS_Multiomics/compare/ldsc_brain_result.RData")
write.csv(ldsc_pbmc_result,file="/share/pub/dengcy/GWAS_Multiomics/compare/ldsc_pbmc_result.csv")

ldsc_pbmc_list<-as.data.frame(lapply(ldsc_pbmc_result,function(df){
    rownames(df)<-df$Name
    df<-df[ldsc_pbmc_result[[1]]$Name,]
  a<- df$Coefficient_P_value
  #a[which(a<4.32)]<-0
  return(a)
}))
rownames(ldsc_pbmc_list)<-ldsc_pbmc_result[[1]]$Name
write.csv(ldsc_pbmc_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/ldsc_pbmc_list.csv")
save(ldsc_pbmc_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/ldsc_pbmc_list.RData")
##############BMMC
setwd("/share2/pub/zhenggw/zhenggw/scPagwas_LDSC/Seu_results")
ldsc_bmmc_result<-lapply(traits,function(i){
    message(i)
   result<-read.table(file=paste0(i,".cell_type_results.txt"),header = T)
    return(result)
})
names(ldsc_bmmc_result)<-traits

ldsc_bmmc_list<-as.data.frame(lapply(ldsc_bmmc_result,function(df){
    rownames(df)<-df$Name
    df<-df[ldsc_bmmc_result[[1]]$Name,]
  a<- df$Coefficient_P_value
  #a[which(a<4.32)]<-0
  return(a)
}))
rownames(ldsc_bmmc_list)<-ldsc_bmmc_result[[1]]$Name
write.csv(ldsc_bmmc_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/ldsc_bmmc_list.csv")
save(ldsc_bmmc_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/ldsc_bmmc_list.RData")

#################rolypoly
##########################magma
load("ldsc_bmmc_list.RData")
result_list<-ldsc_bmmc_list[,3:12]
result_list2<- -log2(result_list)
#magma_pbmc_p

coul <- colorRampPalette(brewer.pal(8, "Oranges"))(25)#

p2star <- function(p){
  symnum(p,cutpoints = c(0,0.001,0.01,0.05,1),
         symbols = c('***','**','*',''),na = NA)
}
hM <- apply(result_list,2, function(x) as.character(p2star(x)))
rownames(hM)<-rownames(result_list2)


result_list2<-as.matrix(result_list2)
result_list2[which(result_list2>20)]<-20

pdf("E:/OneDrive/GWAS_Multiomics/Compare/Figure_ldsc_bmmc_cellytpe_p.pdf",
    height =8,width =6)
heatmap.2(result_list2,
          trace="none",#
          col=coul,#
          density.info = "none",#
          key.xlab ='Correlation',
          key.title = "",
          cexRow = 1,cexCol = 1,
          Rowv = F,Colv =T, #
          margins = c(6, 6),
          cellnote = hM,notecol='black'#
)
dev.off()
```



## Ranked visualized

### 1.Calculated the score

#### 1.1 meta data

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
 library(scPagwas)
 suppressMessages(library(Seurat))
library(org.Hs.eg.db)
library(dplyr)

 traits<-c("eosinophilcount","basophilcount","LymphocytePercent","Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")
g2s=toTable(org.Hs.egSYMBOL)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS")
for(i in traits){
    print(i)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
    Celltype_anno<-Pagwas@misc$Celltype_anno
magma_genes<-read.table(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/",i,"annotated_10kbup_10down.genes.out"),header=T)
#save(magma_genes,file=paste0(i,"_magma_genes.RData"))
magma_genes$gene_id<-magma_genes$GENE
magma_genes=merge(magma_genes,g2s,by="gene_id",all.x=T)
magma_genes<-na.omit(magma_genes)
magma_genes<-magma_genes[order(magma_genes$P,decreasing=F),]
save(magma_genes,file=paste0(i,"_magma_genes.RData"))
scPagwas_genes<-names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation,decreasing=T),])
    magmatop1000<-magma_genes$symbol
magmatop1000<-magmatop1000[1:1000]
scPagwastop1000<-scPagwas_genes[1:1000]
magmatop500<-magmatop1000[1:500]
scPagwastop500<-scPagwas_genes[1:500]
topgene<-list(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500)
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
save(topgene,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_",i,".RData")) 

magmatop1000<-paste(magmatop1000,collapse=",")
scPagwastop1000<-paste(scPagwastop1000,collapse=",")
magmatop500<-paste(magmatop500,collapse=",")
scPagwastop500<-paste(scPagwastop500,collapse=",")

a<-data.frame(genes=c(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500))
rownames(a)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
write.csv(a,file=paste0("scPagwas.magmatopgenes_scRDS",i,".csv"))

meta.data<-data.frame(scPagwas_score=Pagwas$scPagwas_score[colnames(Single_data)],
 CellsrankPvalue=Pagwas$CellsrankPvalue[colnames(Single_data),"pValueHigh"]
                     )
meta.data$positive<-rep("positive",ncol(Single_data))
meta.data$positive[Single_data$CellsrankPvalue>0.05]<-"negative"
save(meta.data,file=paste0(i,"_meta.data.RData"))
    
    }

###############################
for(i in traits){
    print(i)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
n<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")

a<-data.frame(topgene[[1]],topgene[[2]])
colnames(a)<-c("scPagwastop1000","magmatop1000")
write.csv(a,file=paste0("scPagwas.magmatopgenes_scRDS",i,".csv"),row.names = F)
}
#############################

#####load scRDS data
library(SeuratDisk)
library(Seurat)
DefaultAssay(Single_data) <- "RNA"
SaveH5Seurat(Single_data, "Seu_Hema_addata.h5seurat")
Convert("Seu_Hema_addata.h5seurat", dest="h5ad")
write.csv(Single_data@meta.data,file="Seu_Hema.metadata.csv")
```

#### 1.2 scDRS

```python
cd /share/pub/dengcy/GWAS_Multiomics/Pkg/scDRS
git checkout -b rv1 v1.0.0
pip install -e .
######################################
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
DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test"
H5AD_FILE = os.path.join(DATA_PATH, "Seu_Hema_addata.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)

traits = ['eosinophilcount','Lymphocytecount3','monocytecount','neutrophilcount','WhiteBloodCellcount','MeanCorpuscularHemoglobin','MeanCorpusVolume',"basophilcount","LymphocytePercent"]

scdrs.preprocess(adata) 
for i in traits:
 df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_scRDS"+i+".csv", index_col=0)
 df_gs = df_gs.loc[["scPagwastop1000","magmatop1000","scPagwastop500","magmatop500"],:]
 df_gs=  df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".geneset.gs", sep="\t", index=False)
 df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".geneset.gs")
 gene_list  = df_gs['scPagwastop1000'][0]
 gene_weight  = df_gs['scPagwastop1000'][1]
 df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=20)
 df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".1000scPagwas.df_res.csv", sep=",", index=False)
 gene_list  = df_gs['magmatop1000'][0]
 gene_weight  = df_gs['magmatop1000'][1]
 df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=20)
 df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".1000magma.df_res.csv", sep=",", index=False)
 gene_list  = df_gs['scPagwastop500'][0]
 gene_weight  = df_gs['scPagwastop500'][1]
 df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=20)
 df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".500scPagwasdf_res.csv", sep=",", index=False)
 gene_list  = df_gs['magmatop500'][0]
 gene_weight  = df_gs['magmatop500'][1]
 df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=20)
 df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".500magma.df_res.csv", sep=",", index=False)
 print (i)

```

#### 1.3 Other methods：

```R
library(testSctpa)
library(dplyr)
library(GSVA)
library(VISION)
library(AUCell)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/comparegroudtruth")
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")
i<-"monocytecount"
i<-"basophilcount"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_",i,".RData"))
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
Single_data <- Seurat::AddModuleScore(Single_data, assay = "RNA", topgene, name = c("scPagwastop1000.seruatscore","magmatop1000.seruatscore","scPagwastop500.seruatscore","magmatop500.seruatscore"))

auccell = cal_PAS(seurat_object = Single_data,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
gmt_file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_",i,".gmt"))
auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
colnames(auccell_df)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
save(auccell_df,file=paste0(i,"_auccell_df.RData"))

Visiondf = cal_PAS(seurat_object = Single_data,
                       tool = 'Vision',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
gmt_file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_",i,".gmt"))
Visiondf<- as.data.frame(t(GetAssayData(Visiondf, slot="data", assay="PAS")))
colnames(Visiondf)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
save(Visiondf,file=paste0(i,"_Vision_df.RData"))

gsva_df<-gsva(data.matrix(Single_data@assays$RNA@data), gset.idx.list=topgene,method="gsva")
save(df,file=paste0(i,"_gsva_df.RData"))

```

#### Compare with methods and AUC

```R
i<-"monocytecount"
library('ComplexHeatmap')
library(circlize)
library(dplyr)
require(pROC)
require(ggplot2)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_meta.data.RData"))
load("Vision_comparedf.RData")
load("auccell_comparedf.RData")
celltypes= c("12_CD14.Mono.2","13_CD16.Mono","24_CD8.CM","25_NK")

meta.data1<-Single_data@meta.data

df<-meta.data[meta.data$celltypes %in% celltypes,]
df1<-meta.data1[meta.data1$celltypes %in% celltypes,c("scPagwastop1000.seruatscore1","magmatop1000.seruatscore2","scPagwastop500.seruatscore3","magmatop500.seruatscore4")]

df1<-df1[rownames(df),]
auccell_df<-auccell_df[meta.data$celltypes %in% celltypes,]
Visiondf<-Visiondf[meta.data$celltypes %in% celltypes,]
df$celltype<-rep("Correct",nrow(df))
df$celltype[df$celltypes %in% c("24_CD8.CM","25_NK")]<-"non"
df<-cbind(df,df1)
#colors_celltypes=c("#125B50","#E04D01","#F8B400","#FAF5E4","#00AFC1"))

score_df1<-data.frame(scPagwas.scDRS=df$scPagwas1000_scdrs.raw_score,
                scPagwas.Vision=Visiondf$scPagwastop1000,
                scPagwas.auccell=auccell_df$scPagwastop1000,
                scPagwas.seruat=df$scPagwastop1000.seruatscore1,
                magma.scDRS=df$magma1000_scdrs.raw_score,
                magma.auccell=auccell_df$magmatop1000,
                magma.Vision=Visiondf$magmatop1000,
                   magma.seruat=df$magmatop1000.seruatscore2,  
                    index=df$celltype)
auc_list1<-lapply(score_df1[,1:8],function(x){
  roc(predictor=x,response=score_df1$index)  
})

names(auc_list1)<-colnames(score_df1)[1:8]
 auc_l1<- unlist(lapply(1:length(auc_list1),function(x) round(as.numeric(auc_list1[[x]]["auc"]),3)))
 
 names(auc_list1)<- paste0(names(auc_list1),"(AUC=",auc_l1,")")
 
 pdf("AUC_pagwastop1000_comparescpre.pdf",height = 6)
 ggroc(auc_list1, linetype = 2, size = 1,alpha=0.8)+
    ggtitle("ROC curve for top 1000 genes for monocytecount traits") + 
    theme_classic()+ggsci::scale_color_lancet()+
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
                 color="grey", linetype="dashed")
    #theme(panel.background=element_rect(fill="white",colour="blue"))
  dev.off()

                        
score_df2<-data.frame(scPagwas.scDRS=df$scPagwas500_scdrs.raw_score,
                scPagwas.Vision=Visiondf$scPagwastop500,
                scPagwas.auccell=auccell_df$scPagwastop500,
                scPagwas.seruat=df$scPagwastop500.seruatscore3,
                magma.scDRS=df$magma500_scdrs.raw_score,
                magma.auccell=auccell_df$magmatop500,
                magma.Vision=Visiondf$magmatop500,
                   magma.seruat=df$magmatop500.seruatscore4,  
                    index=df$celltype)
auc_list2<-lapply(score_df2[,1:8],function(x){
  roc(predictor=x,response=score_df2$index)  
})

names(auc_list2)<-colnames(score_df2)[1:8]
 auc_l2<- unlist(lapply(1:length(auc_list2),function(x) round(as.numeric(auc_list2[[x]]["auc"]),3)))
 
 names(auc_list2)<- paste0(names(auc_list2),"(AUC=",auc_l2,")")
 
 pdf("AUC_pagwastop500_comparescpre.pdf",height = 6)
 ggroc(auc_list2, linetype = 2, size = 1,alpha=0.8)+
    ggtitle("ROC curve for top 500 genes for monocytecount traits") + 
    theme_classic()+ggsci::scale_color_lancet()+
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
                 color="grey", linetype="dashed")
  dev.off()
save(df,file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/comparegroudtruth/monocytecounts_metadf.RData")
```

### 2.groudtruth

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/comparegroudtruth")
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")

traits<-c("eosinophilcount","basophilcount","LymphocytePercent","lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","MeanCorpuscularHemoglobin","MeanCorpusVolume")
for(i in traits){
 load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_meta.data.RData"))
pagwas1000_scDRS_re<-read.csv(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,".1000scPagwas.df_res.csv"))
magma1000_scDRS_re<-read.csv(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,".1000magma.df_res.csv"))
 pagwas500_scDRS_re<-read.csv(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,".500scPagwasdf_res.csv"))
magma500_scDRS_re<-read.csv(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,".500magma.df_res.csv"))
    
meta.data$scPagwas1000_scdrs.raw_score<-pagwas1000_scDRS_re$raw_score
meta.data$scPagwas1000_scdrs.zscore<-pagwas1000_scDRS_re$zscore
 meta.data$scPagwas1000_scdrs.pval<-pagwas1000_scDRS_re$pval
meta.data$magma1000_scdrs.raw_score<-magma1000_scDRS_re$raw_score
meta.data$magma1000_scdrs.zscore<-magma1000_scDRS_re$zscore
meta.data$magma1000_scdrs.pval<-magma1000_scDRS_re$pval

meta.data$scPagwas500_scdrs.raw_score<-pagwas500_scDRS_re$raw_score
meta.data$scPagwas500_scdrs.zscore<-pagwas500_scDRS_re$zscore
 meta.data$scPagwas500_scdrs.pval<-pagwas500_scDRS_re$pval
meta.data$magma500_scdrs.raw_score<-magma500_scDRS_re$raw_score
meta.data$magma500_scdrs.zscore<-magma500_scDRS_re$zscore
meta.data$magma500_scdrs.pval<-magma500_scDRS_re$pval
    
meta.data$celltypes<-Single_data$BioClassification
save(meta.data,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_meta.data.RData"))  
}

```

#### Compare with different methods

#### Ranked proportion：

```R
i<-"monocytecount"
library(ggpubr)
library(ggplot2)
load(paste0("E:/OneDrive/GWAS_Multiomics/Compare/5.16compareresult/",i,"_meta.data.RData"))
###########################
filecell="monocytecount"
celltypes= c("12_CD14.Mono.2","13_CD16.Mono","24_CD8.CM","25_NK")
colors_celltypes=c("#F38BA0","#FFBCBC","#EDF6E5","#B5EAEA")
df<-meta.data[meta.data$celltypes %in% celltypes,]
#df<-df[order(df$celltypes),]

###########################
filecell="MeanCorpusVolume"
celltypes=c("01_HSC","03_Late.Eryth","02_Early.Eryth","24_CD8.CM","22_CD4.M")
#colors_celltypes=c("#F38BA0","#FFBCBC","#EDF6E5","#B5EAEA")
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",filecell,"_meta.data.RData"))
df<-meta.data[meta.data$celltypes %in% celltypes,]

df1<-df[,c("scPagwas1000_scdrs.raw_score",
           "magma1000_scdrs.raw_score",
           "scPagwas500_scdrs.raw_score",
           "magma500_scdrs.raw_score","celltypes")]
  
df1<-df1[order(df1$scPagwas1000_scdrs.raw_score,decreasing=T),]
#a<-scPagwas_magma_genelist[[i]]
n<-nrow(df1)
b<-rep(1,n)
b[1:(0.25*n)]<-1
b[(0.25*n+1):(0.5*n)]<-2
b[(0.5*n+1):(0.75*n)]<-3
b[(0.75*n+1):n]<-4
df1$rg<-b
df1$celltype<-rep("Correct",nrow(df1))
df1$celltype[df1$celltypes %in% c("24_CD8.CM","22_CD4.M")]<-"non"

print(table(data.frame(df1$celltype,df1$rg)))
t1<-table(data.frame(df1$celltype,df1$rg))
percent1<-t1[,1]/sum(t1[,1])
#   Correct         non 
#0.998688811 0.001311189

df <- data.frame(
  group = names(percent1),
  value = percent1)

p1<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4") )

percent2<-t1[,2]/sum(t1[,2])

df <- data.frame(
  group = names(percent2),
  value = percent2)

p2<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4"))

percent3<-t1[,3]/sum(t1[,3])

df <- data.frame(
  group = names(percent3),
  value = percent3)

p3<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4"))

percent4<-t1[,4]/sum(t1[,4])

df <- data.frame(
  group = names(percent4),
  value = percent4)

p4<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4"))

setwd("E:/OneDrive/GWAS_Multiomics/Compare/5.16compareresult")

#####################500
df1<-df1[order(df1$scPagwas500_scdrs.raw_score,decreasing=T),]
#a<-scPagwas_magma_genelist[[i]]
n<-nrow(df1)
b<-rep(1,n)
b[1:(0.25*n)]<-1
b[(0.25*n+1):(0.5*n)]<-2
b[(0.5*n+1):(0.75*n)]<-3
b[(0.75*n+1):n]<-4
df1$rg<-b
df1$celltype<-rep("Correct",nrow(df1))
df1$celltype[df1$celltypes %in% c("24_CD8.CM","22_CD4.M")]<-"non"

print(table(data.frame(df1$celltype,df1$rg)))
t1<-table(data.frame(df1$celltype,df1$rg))
percent1<-t1[,1]/sum(t1[,1])

df <- data.frame(
  group = names(percent1),
  value = percent1)

p1.1<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4") )

percent2<-t1[,2]/sum(t1[,2])

df <- data.frame(
  group = names(percent2),
  value = percent2)

p1.2<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4"))

percent3<-t1[,3]/sum(t1[,3])

df <- data.frame(
  group = names(percent3),
  value = percent3)

p1.3<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4"))

percent4<-t1[,4]/sum(t1[,4])

df <- data.frame(
  group = names(percent4),
  value = percent4)

p1.4<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#DF7861", "#D4E2D4"))

#####################magma

df1<-df1[order(df1$magma1000_scdrs.raw_score,decreasing=T),]
#a<-scPagwas_magma_genelist[[i]]
n<-nrow(df1)
b<-rep(1,n)
b[1:(0.25*n)]<-1
b[(0.25*n+1):(0.5*n)]<-2
b[(0.5*n+1):(0.75*n)]<-3
b[(0.75*n+1):n]<-4
df1$rg<-b
df1$celltype<-rep("Correct",nrow(df1))
df1$celltype[df1$celltypes %in% c("24_CD8.CM","22_CD4.M")]<-"non"

print(table(data.frame(df1$celltype,df1$rg)))
t1<-table(data.frame(df1$celltype,df1$rg))
percent1<-t1[,1]/sum(t1[,1])
#Correct        non 
#0.91258741 0.08741259
df <- data.frame(
  group = names(percent1),
  value = percent1)

p2.1<-ggdonutchart(df, "value", label = "group",
                   fill = "group", color = "white",
                   palette = c("#DF7861", "#D4E2D4") )

percent2<-t1[,2]/sum(t1[,2])
#Correct       non 
#0.3041575 0.6958425
df <- data.frame(
  group = names(percent2),
  value = percent2)

p2.2<-ggdonutchart(df, "value", label = "group",
                   fill = "group", color = "white",
                   palette = c("#DF7861", "#D4E2D4"))

percent3<-t1[,3]/sum(t1[,3])
#Correct       non 
#0.1540481 0.8459519

df <- data.frame(
  group = names(percent3),
  value = percent3)

p2.3<-ggdonutchart(df, "value", label = "group",
                   fill = "group", color = "white",
                   palette = c("#DF7861", "#D4E2D4"))

percent4<-t1[,4]/sum(t1[,4])

df <- data.frame(
  group = names(percent4),
  value = percent4)

p2.4<-ggdonutchart(df, "value", label = "group",
                   fill = "group", color = "white",
                   palette = c("#DF7861", "#D4E2D4"))


######################500

df1<-df1[order(df1$magma500_scdrs.raw_score,decreasing=T),]
#a<-scPagwas_magma_genelist[[i]]
n<-nrow(df1)
b<-rep(1,n)
b[1:(0.25*n)]<-1
b[(0.25*n+1):(0.5*n)]<-2
b[(0.5*n+1):(0.75*n)]<-3
b[(0.75*n+1):n]<-4
df1$rg<-b
df1$celltype<-rep("Correct",nrow(df1))
df1$celltype[df1$celltypes %in% c("24_CD8.CM","22_CD4.M")]<-"non"

print(table(data.frame(df1$celltype,df1$rg)))
t1<-table(data.frame(df1$celltype,df1$rg))
percent1<-t1[,1]/sum(t1[,1])

df <- data.frame(
  group = names(percent1),
  value = percent1)

p3.1<-ggdonutchart(df, "value", label = "group",
                   fill = "group", color = "white",
                   palette = c("#DF7861", "#D4E2D4") )

percent2<-t1[,2]/sum(t1[,2])

df <- data.frame(
  group = names(percent2),
  value = percent2)

p3.2<-ggdonutchart(df, "value", label = "group",
                   fill = "group", color = "white",
                   palette = c("#DF7861", "#D4E2D4"))

percent3<-t1[,3]/sum(t1[,3])

df <- data.frame(
  group = names(percent3),
  value = percent3)

p3.3<-ggdonutchart(df, "value", label = "group",
                   fill = "group", color = "white",
                   palette = c("#DF7861", "#D4E2D4"))

percent4<-t1[,4]/sum(t1[,4])

df <- data.frame(
  group = names(percent4),
  value = percent4)

p3.4<-ggdonutchart(df, "value", label = "group",
                   fill = "group", color = "white",
                   palette = c("#DF7861", "#D4E2D4"))
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/comparegroudtruth")
pdf(paste0(filecell,"percent.1000_scdrs.pdf"),width = 10,height = 10)
ggpubr::ggarrange(p1,p2,p3,p4,p2.1,p2.2,p2.3,p2.4,nrow = 2,ncol = 4)
dev.off()

pdf(paste0(filecell,"percent.500_scdrs.pdf"),width = 10,height = 10)
ggpubr::ggarrange(p1.1,p1.2,p1.3,p1.4,p3.1,p3.2,p3.3,p3.4,nrow = 2,ncol = 4)
dev.off()
```



### 3.groudtruth : gene ranked plot

#### top gene enrichment analysis

```R

library(ggplot2)
library(reshape2)
library(ggpubr)
traits<-c("eosinophilcount","basophilcount","LymphocytePercent","Lymphocytecount3",
             "monocytecount","neutrophilcount","WhiteBloodCellcount","Hemoglobinconcen",
             "MeanCorpuscularHemoglobin","MeanCorpusVolume")

lapply(traits,function(i){
  print(i)
  if(i=="MeanCorpuscularHemoglobin"){
    magma_goresult<-NULL
    scPagwas_goresult<-read.delim(paste0("E:/OneDrive/GWAS_Multiomics/Compare/goanalysis/",i,"scPagwastop1000/enrichment_results_wg_result.txt"),header = T)
    a1<-0
    a2<-sum(scPagwas_goresult$FDR<0.01)
    a3<-0
    a4<-sum(scPagwas_goresult$FDR<0.001)
    a5<-0
    a6<-sum(scPagwas_goresult$FDR<0.0001)
    a7<-0
    a8<-sum(scPagwas_goresult$FDR<0.00001)
    
  }else{
    magma_goresult<-read.delim(paste0("E:/OneDrive/GWAS_Multiomics/Compare/goanalysis/",i,"magmatop1000/enrichment_results_wg_result.txt"),header = T)
    scPagwas_goresult<-read.delim(paste0("E:/OneDrive/GWAS_Multiomics/Compare/goanalysis/",i,"scPagwastop1000/enrichment_results_wg_result.txt"),header = T)
    a1<-sum(magma_goresult$FDR<0.01)
    a2<-sum(scPagwas_goresult$FDR<0.01)
    a3<-sum(magma_goresult$FDR<0.001)
    a4<-sum(scPagwas_goresult$FDR<0.001)
    a5<-sum(magma_goresult$FDR<0.0001)
    a6<-sum(scPagwas_goresult$FDR<0.0001)
    a7<-sum(magma_goresult$FDR<0.00001)
    a8<-sum(scPagwas_goresult$FDR<0.00001)
    
  }
  
  go_df <- data.frame(magma=c(a1,a3,a5,a7),scPagwas=c(a2,a4,a6,a8),FDR=c("1e-02","1e-03","1e-04","1e-05"))
  ##画图，棒棒图
  gg_go_df<-melt(go_df,id.vars = "FDR")
  gg_go_df$FDR<-factor(gg_go_df$FDR,levels = c("1e-05","1e-04","1e-03","1e-02"))
  
  setwd("E:/OneDrive/GWAS_Multiomics/Compare/goanalysis/")
  pdf(paste0("dorplot_",i,"_goranalysis_allgenes.pdf"),height =3,width=5)
  print(ggdotchart(gg_go_df, x="FDR", y="value", color = "variable",          
             palette = c("#1F4690","#FFA500"), # 配色
             sorting = "none",    # 排序  
             size =1,dot.size=5,
             rotate = T,label="value",
             add = "segments", #添加棒棒
             main=i, ylab="Numbers of significant GO terms",
             ggtheme = theme_pubr()) )
  
  dev.off()
  
})

lapply(traits,function(i){
 if(i=="MeanCorpuscularHemoglobin"){
    magma_goresult<-NULL
    scPagwas_goresult<-read.delim(paste0("E:/OneDrive/GWAS_Multiomics/Compare/goanalysis/",i,"scPagwastop1000/enrichment_results_wg_result.txt"),header = T)
    write.csv(scPagwas_goresult,file=paste0("E:/OneDrive/GWAS_Multiomics/Manuscripts/Supplementfiles/top1000goanalysis/",i,"scPagwas_top1000genes_goresult.csv"))

  }else{
    magma_goresult<-read.delim(paste0("E:/OneDrive/GWAS_Multiomics/Compare/goanalysis/",i,"magmatop1000/enrichment_results_wg_result.txt"),header = T)
    scPagwas_goresult<-read.delim(paste0("E:/OneDrive/GWAS_Multiomics/Compare/goanalysis/",i,"scPagwastop1000/enrichment_results_wg_result.txt"),header = T)
        write.csv(magma_goresult,file=paste0("E:/OneDrive/GWAS_Multiomics/Manuscripts/Supplementfiles/top1000goanalysis/",i,"magma_top1000genes_goresult.csv"))
        write.csv(scPagwas_goresult,file=paste0("E:/OneDrive/GWAS_Multiomics/Manuscripts/Supplementfiles/top1000goanalysis/",i,"scPagwas_top1000genes_goresult.csv"))

}
})

```

#### Enrichment the ranked gene

```R

library(scPagwas)
suppressMessages(library(Seurat))
traits<-c("eosinophilcount","basophilcount","LymphocytePercent","Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/")
for(i in traits){
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_",i,".RData"))
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
}

load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
geneList<-Pagwas@misc$gene_heritability_correlation[,1]
geneList=sort(geneList,decreasing = T) 
select_cells<-c("01_HSC","02_Early.Eryth","05_CMP.LMPP","09_pDC","12_CD14.Mono.2","19_CD8.N")
set.seed(1234)
a_list<-lapply(select_cells,function(x){
  a<-subset(Pagwas,idents =x) 
  a<-a[,sample(1:ncol(a),500)]
    return(a)
})
Pagwas_rank<-merge(x = a_list[[1]],y =  a_list[2:length(select_cells)])
save(Pagwas_rank,file="Pagwas_rank.RData")
data_mat<- GetAssayData(Pagwas_rank,assay ="RNA")
##################
a<-apply(data_mat,1,mean)
a<-sort(a,decreasing=T)                
rank_gene<-names(a)
save(rank_gene,file="expr_rank_gene.RData")
g1<-intersect(rank_gene[1:(1/2*length(a))],topgene$scPagwastop1000)
g2<-intersect(rank_gene[(1/2*length(a)):length(a)],topgene$scPagwastop1000)

m1<-intersect(rank_gene[1:(1/2*length(a))],topgene$magmatop1000)
m2<-intersect(rank_gene[(1/2*length(a)):length(a)],topgene$magmatop1000)

gene_list<-list(g1,g2,m1,m2)
names(gene_list)<-c("g1","g2","m1","m2")
lapply(names(gene_list), function(x){
  write.csv(gene_list[[x]],file=paste0(x,"split2_gene.csv"),row.names = F)
})
write.csv(degene,file="de_gene_monocyte.csv")
```

#### Ranked gene visualize

```R
setwd("E:/OneDrive/GWAS_Multiomics/Compare/5.16compareresult")
load("D:/OneDrive/GWAS_Multiomics/Compare/goanalysis/expr_rank_gene.RData")
load("D:/OneDrive/GWAS_Multiomics/Compare/goanalysis/Lymphocytecount3_genes_split2.RData")

percent1<-c(length(c(gene_list$g1))/length(unlist(gene_list[1:2])),
            length(c(gene_list$g2))/length(unlist(gene_list[1:2]))
)
percent2<-c(length(c(gene_list$m1))/length(unlist(gene_list[3:4])),
            length(c(gene_list$m2))/length(unlist(gene_list[3:4]))
)

b1<-rep(0,length(rank_gene))
b1[rank_gene %in% unlist(gene_list[1:2])]<-1
b2<-rep(0,length(rank_gene))
b2[rank_gene %in% unlist(gene_list[3:4])]<-1

a<-c(rep(1,length(rank_gene)*0.5),rep(2,length(rank_gene)*0.5))
rank_df<-data.frame(rank_gene,a,b1,b2)

library(ggpubr)
df <- data.frame(
  group = c("scPagwas","others"),
  value = c(percent1[1],1-percent1[1]))

p1<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#FFE194", "#4C4C6D") )
df <- data.frame(
  group = c("scPagwas","others"),
  value = c(percent1[2],1-percent1[2]))

p2<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#FFE194", "#4C4C6D")  )  
df <- data.frame(
  group = c("scPagwas","others"),
  value = c(percent2[1],1-percent2[1]))
p3<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#FFE194", "#4C4C6D")  )  

df <- data.frame(
  group = c("scPagwas","others"),
  value = c(percent2[2],1-percent2[2]))
p4<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#FFE194", "#4C4C6D")  )  
setwd("D:/OneDrive/GWAS_Multiomics/Compare")

pdf("percent.for.toprankgenes.lymphocyte.pdf",width = 10,height = 10)
ggpubr::ggarrange(p1,p2,p3,p4,nrow = 2,ncol = 4)
dev.off()

setwd("D:/OneDrive/GWAS_Multiomics/Compare")
library('ComplexHeatmap')
library(circlize)
col_fun = colorRamp2(c(0, 1), c("#F6F6F6", "#161D6F"))
col_fun(seq(0, 1))

pdf("lymphocyte.heatmap_rank_PAGWASgene.pdf",width = 2)
Heatmap(data.matrix(rank_df$b1),
        col = col_fun,
        #left_annotation = ha, 
        cluster_columns = F,
        cluster_rows = F,
        color_space="HLS",
        border=T,
        row_gap = unit(0.25, "mm"),
        show_parent_dend_line=T,
        #name = "-log2(p)",
        #row_order =order(rdf$types),
        show_row_names=T,
        show_column_names=T
        #row_split=rdf$phenotypes
)
dev.off()

pdf("lymphocyte.heatmap_rank_MAGMAgene.pdf",width = 2)
Heatmap(data.matrix(rank_df$b2),
        col = col_fun,
        #left_annotation = ha, 
        cluster_columns = F,
        cluster_rows = F,
        color_space="HLS",
        border=T,
        row_gap = unit(0.25, "mm"),
        show_parent_dend_line=T,
        #name = "-log2(p)",
        #row_order =order(rdf$types),
        show_row_names=T,
        show_column_names=T
        #row_split=rdf$phenotypes
)
dev.off()

```

### BMMC mean score for celltypes.plot

```R
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))

setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_pbmc.rds")

 traits<-c("eosinophilcount","basophilcount","LymphocytePercent","Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","MeanCorpuscularHemoglobin","MeanCorpusVolume")

 i<-traits[1]
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS")
 for(i in traits){
 load(paste0(i,"_meta.data.RData"))

 }
setwd("E:/OneDrive/GWAS_Multiomics/Compare")
load("scPagwas_Seu_Healthy_resultlist.RData")

scPagwas_Hema_p2<-as.data.frame(lapply(result_list,function(df){
  a<- -log2(df)
  a[which(a<4.32)]<-0
  return(a)
}))
rownames(scPagwas_Hema_p)<- rownames(result_list)

library('ComplexHeatmap')
library(circlize)
col_fun = colorRamp2(c(0, 10), c("#FBF8F1", "#FC4F4F"))
col_fun(seq(0, 10))
#magma_pbmc_p

pdf("E:/OneDrive/GWAS_Multiomics/Compare/Figure_scPagwas_hema_bloodtraits1.pdf",
    height =8,width =7)
Heatmap(as.matrix(scPagwas_Hema_p),
        col = col_fun,
        #left_annotation = ha, 
        cluster_columns = T,
        cluster_rows = T,
        color_space="HLS",
        border=T,
        row_gap = unit(0.25, "mm"),
        show_parent_dend_line=T,
        name = "-log2(p)",
        #row_order =order(rdf$types),
        show_row_names=T,
        show_column_names=T
        #row_split=rdf$phenotypes
)
dev.off()
```



