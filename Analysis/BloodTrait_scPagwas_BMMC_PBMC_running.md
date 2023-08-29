
# GWAS summary statistics and fine-mapping analysis 
Blood cell traits. 
Summary statistics of 22 blood cell traits from Blood Cell Consortium 2 (BCX2) were processed 
as previously described

## blood traits：
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

load("./phe_NM_Healthy_pbmc.Rdata")
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
saveRDS(NM_Healthy_pbmc,file="./NM_Healthy_pbmc.rds")
```

#### BMMC single cell data

```R
library("scPagwas")
library("Seurat")
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

#        01_HSC 02_Early.Eryth  03_Late.Eryth  04_Early.Baso    05_CMP.LMPP       06_CLP.1 
#          1425           1653            446            111           2260            903 
#        07_GMP    08_GMP.Neut         09_pDC         10_cDC 11_CD14.Mono.1 12_CD14.Mono.2 
#          2097           1050            544            325           1800           4222 
#  13_CD16.Mono         14_Unk       15_CLP.2       16_Pre.B           17_B      18_Plasma 
#           292            520            377            710           1711             62 
#      19_CD8.N      20_CD4.N1      21_CD4.N2       22_CD4.M      23_CD8.EM      24_CD8.CM 
#          1521           2470           2364           3539            796           2080 
#         25_NK         26_Unk 
#          2143            161 

Seu_Healthy_Hema <- ScaleData(Seu_Healthy_Hema)
Seu_Healthy_Hema <- NormalizeData(Seu_Healthy_Hema, normalization.method = "LogNormalize", scale.factor = 10000)

```

#### Run the single_data_input

```R
####PBMC
setwd("/share/pub/dengcy/GWAS_Multiomics/compare")
 library(scPagwas)
 suppressMessages(library(Seurat))
Single_data=readRDS("./NM_Healthy_pbmc.rds")
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
                      assay="RNA",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas.RData"))
 }
```

#### PBMC

```R
library(scPagwas)
library(ggplot2)
suppressMessages(library(Seurat))
suppressMessages(library("dplyr"))
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc")
load("/share/pub/dengcy/GWAS_Multiomics/compare/NM_Healthy_pbmc_prePagwas.RData")

traits<-c("basophilcount","eosinophilcount" ,"Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

for(i in traits){
     Pagwas<-scPagwas_main(Pagwas =Pagwas,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"),
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_pbmc.rds",
                     output.prefix=i,
                     Pathway_list=Genes_by_pathway_kegg,
                     output.dirs=paste0(i,"_pbmc_scPagwas"),
                     ncores=1,
                     assay="RNA",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/",i,"_pbmc_scPagwas_singlecell.RData"))
}
```

### 4.Integrate bmmc and PBMC celltype result

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
#"RBCcount","Plateletcount",
traits<-c("basophilcount","eosinophilcount" ,"Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")
result_list<-lapply(traits,function(i){
   print(i)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas.RData"))
    a<-Pagwas@misc$bootstrap_results[-1,]
    return(a$bp_value)
})
names(result_list)<-traits
result_list<-as.data.frame(result_list)
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

## Blood traits run MAGMA for trait-relevant celltypes pvalue result

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
for i in  basophilcount eosinophilcount lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount LymphocytePercent Hemoglobinconcen MeanCorpuscularHemoglobin MeanCorpusVolume
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
for i in basophilcount basophilcount eosinophilcount lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount LymphocytePercent Hemoglobinconcen MeanCorpuscularHemoglobin MeanCorpusVolume 
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

```shell
#1.
cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
for i in monocytecount neutrophilcount WhiteBloodCellcount LymphocytePercent Hemoglobinconcen MeanCorpuscularHemoglobin MeanCorpusVolume
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
```
## Ranked visualized data

### 1.Calculated the score

#### 1.1 meta data

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
library(scPagwas)
suppressMessages(library(Seurat))
library(org.Hs.eg.db)
library(dplyr)

traits<-c("eosinophilcount","basophilcount","LymphocytePercent","Lymphocytecount3","Lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")
g2s=toTable(org.Hs.egSYMBOL)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS")
for(i in traits){
  print(i)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas.RData"))
  Celltype_anno<-Pagwas@misc$Celltype_anno
  magma_genes<-read.table(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/magma/",i,"_10kbup_10down.genes.out"),header=T)
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

for(i in traits){
  print(i)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas.RData"))
Pagwas$scPagwas_score

###############################
for(i in traits){
    print(i)
    load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas.RData"))
    n<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")

    a<-data.frame(topgene[[1]],topgene[[2]])
    colnames(a)<-c("scPagwastop1000","magmatop1000")
    write.csv(a,file=paste0("scPagwas.magmatopgenes_scRDS",i,".csv"),row.names = F)
}

#####load scRDS data
library(SeuratDisk)
library(Seurat)
DefaultAssay(Single_data) <- "RNA"
SaveH5Seurat(Single_data, "Seu_Hema_addata.h5seurat")
Convert("Seu_Hema_addata.h5seurat", dest="h5ad")
write.csv(Single_data@meta.data,file="Seu_Hema.metadata.csv")

#########cov files
for(i in c("Lymphocytecount3","monocytecount","MeanCorpusVolume")){
    print(i)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
a<-data.frame(index=rownames(Pagwas@meta.data),const=rep(1,ncol(Pagwas)))
write.table(a,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scRDS_",i,".cov"),row.names = F,quote=F,sep="\t")
}

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

traits = ['Lymphocytecount3']
#['eosinophilcount','Lymphocytecount3','monocytecount','neutrophilcount','WhiteBloodCellcount','MeanCorpuscularHemoglobin','MeanCorpusVolume',"basophilcount","LymphocytePercent"]

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

```

cd /share/pub/dengcy/GWAS_Multiomics/Pkg/scDRS
python compute_score.py \
    --h5ad_file /share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Seu_Hema_addata.h5ad\
    --h5ad_species human\
    --cov_file /share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scRDS_monocytecount.cov\
    --gs_file /share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/monocytecount.geneset.gs\
    --gs_species human\
    --n_ctrl 20\
    --out_folder /share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/


python compute_downstream.py \
    --h5ad_file ${h5ad_file}.h5ad \
    --score_file @.full_score.gz \
    --cell_type celltypes \
    --cell_variable causal_variable,non_causal_variable,covariate\
    --flag_gene True\
    --flag_filter False\
    --flag_raw_count False\ # flag_raw_count is set to `False` because the toy data is already log-normalized, set to `True` if your data is not log-normalized
    --out_folder ${out_dir}
```



#### 1.3 Other methods：

```R
library(testSctpa)
library(dplyr)
library(GSVA)
library(VISION)
library(AUCell)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")
i<-"monocytecount"
i<-"Lymphocytecount3"
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

### 2.Scdrs top 1000genes score for magma and scPagwas

```R
#setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/comparegroudtruth")
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")

#traits<-c("eosinophilcount","basophilcount","LymphocytePercent","lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","MeanCorpuscularHemoglobin","MeanCorpusVolume")
traits<-"Lymphocytecount3"
load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas_v1.10.0.RData")
for(i in traits){
  #  i<-traits
# load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_meta.data.RData"))
    
pagwas1000_scDRS_re<-read.csv(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,".1000scPagwas.df_res.csv"))
magma1000_scDRS_re<-read.csv(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,".1000magma.df_res.csv"))
 pagwas500_scDRS_re<-read.csv(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,".500scPagwasdf_res.csv"))
magma500_scDRS_re<-read.csv(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,".500magma.df_res.csv"))
    
Pagwas$scPagwas1000_scdrs.raw_score<-pagwas1000_scDRS_re$raw_score
Pagwas$scPagwas1000_scdrs.zscore<-pagwas1000_scDRS_re$zscore
 Pagwas$scPagwas1000_scdrs.pval<-pagwas1000_scDRS_re$pval
Pagwas$magma1000_scdrs.raw_score<-magma1000_scDRS_re$raw_score
Pagwas$magma1000_scdrs.zscore<-magma1000_scDRS_re$zscore
Pagwas$magma1000_scdrs.pval<-magma1000_scDRS_re$pval

Pagwas$scPagwas500_scdrs.raw_score<-pagwas500_scDRS_re$raw_score
Pagwas$scPagwas500_scdrs.zscore<-pagwas500_scDRS_re$zscore
 Pagwas$scPagwas500_scdrs.pval<-pagwas500_scDRS_re$pval
Pagwas$magma500_scdrs.raw_score<-magma500_scDRS_re$raw_score
Pagwas$magma500_scdrs.zscore<-magma500_scDRS_re$zscore
Pagwas$magma500_scdrs.pval<-magma500_scDRS_re$pval
    
#Pagwas$celltypes<-Single_data$BioClassification
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas.RData")  
}

```
