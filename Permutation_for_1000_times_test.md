# permutation for 1000 times test

**Supplementary Figure S5**

## function for permutation

```R
Pagwas_chunyu_permutation <- function(data_su,seed){
	set.seed(seed)

if("OR" %in% colnames(data_su)) {sd_1 <-sd(log(data_su$OR))}
if("beta" %in% colnames(data_su)){sd_1 <-sd(data_su$beta)}
mean_1 <- 0
N = 10000000
y <- rnorm(N,mean = mean_1, sd = sd_1 )
data2 <- as.data.frame(y)
#Beta distribution for se
N <- 10000000

if("se" %in% colnames(data_su)){
m <-mean(data_su$se)
v <-var(data_su$se)}
if("SE" %in% colnames(data_su)){
m<-mean(data_su$SE)
v<-var(data_su$SE)
}
alpha <- (m*(m-2*m^2+m^3-v+m*v))/((1-m)*v)
beta <- (m-2*m^2+m^3-v+m*v)/v


alpha<- abs(alpha)
beta<-abs(beta)
dd<-as.data.frame(rbeta(N,alpha, beta))

N <- 10000000
maf_data <- as.data.frame(runif(N,min=0.1,max=1))
#Extract random se and Beta
num <- length(data_su$se)
beta_new <- sample(data2[,1],num,replace = FALSE)
SE_new <- sample(dd[,1],num,replace = FALSE)
MAF_new <- sample(maf_data[,1],num,replace = FALSE)

data_su_subset <- data_su[,1:3]

permuted_data_for_SCZ <- data_su_subset
permuted_data_for_SCZ$beta <- beta_new 
permuted_data_for_SCZ$se <- SE_new
permuted_data_for_SCZ$maf <-MAF_new 

return(permuted_data_for_SCZ)
}

```



## Run sample for different seeds

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/test/1000test")
library("dplyr")
library("foreach")
library("data.table")
library("Seurat")
library("Matrix")
library(stringr) 
library("irlba")
library("glmnet")
library(GenomicRanges)
library(utils)
library(ggplot2) 
library(ggthemes) 
library(ggpubr)

load("/share/pub/dengcy/GWAS_Multiomics/test/brain/mouse_brain_seu2.RData")
load("/share/pub/dengcy/GWAS_Multiomics/test/brain/gwas_files/gwas_PD.RData")
set.seed(1234)
gwas_PD_test<-gwas_PD[!duplicated(gwas_PD$rsid),]
gwas_PD_test<-gwas_PD_test[sample(1:nrow(gwas_PD_test),1000000),]
##use all gwas to run scpagwas
Pagwas <- Pagwas_main (Pagwas = NULL,
                       gwas_data = gwas_PD_test,
                       add_eqtls="FALSE",
                       block_annotation = block_annotation,
                       Single_data = mouse_brain_seu2,
                       Pathway_list= genes.by.pathway_kegg,
                       chrom_ld = chrom_ld)
save(Pagwas,file="mouse_FALSE_new.RData")

load("mouse_FALSE_new.RData")
source("/share/pub/dengcy/GWAS_Multiomics/test/1000test/Pagwas_chunyu_permutation.R")

timestart<-Sys.time()
perm_results1<-lapply(1101:1200,function(seed){
PD_test<-Pagwas_chunyu_permutation(data_su=gwas_PD_test,seed=seed)
Pagwas<-scPagwas_main(Pagwas =Pagwas,
                     gwas_data =PD_test,
                     Single_data =mouse_brain_seu2,
                     output.prefix=i,
                     celltype=T,
                     singlecell=F,
                     Pathway_list=Genes_by_pathway_kegg,
                     output.dirs=paste0(i,"_pbmc_scPagwasv1.9"),
                     ncores=1,
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  bootstrap_results<- Pagwas$bootstrap_results
  bootstrap_results<-bootstrap_results[-1,]
  bootstrap_results$bp_value<- p.adjust(bootstrap_results$bp_value)
  return(bootstrap_results)
  })
save(perm_results1,file="Mouse_perm_results1.RData")
 timeend<-Sys.time()
 runningtime<-timeend-timestart
 print(runningtime) 
```

