 血细胞trait进行scPagwas的计算

[toc]

<html>
<!--在这里插入内容-->
</html>


于福龙文章中：
GWAS summary statistics and fine-mapping analysis 
Blood cell traits. 
Summary statistics of 22 blood cell traits from Blood Cell Consortium 2 (BCX2) were processed 
as previously described
> (Vuckovic, D. et al. The Polygenic and Monogenic Basis of Blood Traits and Diseases. Cell 182, 1214–1231.e11 (2020)). 

Variants with fine-mapped posterior probability > 0.001 for a locus in 
one or more blood traits were retained and used for analysis.
### blood traits下载地址：
Data and Code Availability Summary statistics are available to download from: ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/UKBB_blood_cell_traits/ for UK Biobank and http://www.mhi-humangenetics.org/en/resources for the meta-analysis. The accession numbers for the UK Biobank summary statistics reported in this paper are GWAS Catalog: GCST90002379–GCST90002407.
The code generated during this study is publicly available at GitHub https://github.com/bloodcellgwas/manuscript_code/.

率先下载这三种血细胞的数据：
![image](BAB3ABCD7BAC4433989F16667ABFE84C)

各种缩写代表的意义：

英文 | 中文
---|---
red blood cell count (RBC count) | row 1 col 2
hemoglobin concentration (HGB) | row 2 col 2
hematocrit (HCT) | 
mean corpuscular hemoglobin (MCH) | 
mean corpuscular volume (MCV) | 
mean corpuscular hemoglobin concentration (MCHC) | 
RBC distribution width (RDW) | 
total white blood cell count (WBC count) | 
neutrophil count (Neutro) | 
lymphocyte count (Lympho) | 
monocyte count (Mono) | 
basophil count (Baso) | 
eosinophil count (Eosin)| 
platelet count (PLT count)| 
mean platelet volume (MPV) | 

- basophil cell count
  Dataset: ieu-b-29
- eosinophil cell count
  Dataset: ieu-b-33
- lymphocyte cell count
  Dataset: ieu-b-32
- monocyte cell count
  Dataset: ieu-b-31
- neutrophil cell count
  Dataset: ieu-b-34
- white blood cell count
  Dataset: ieu-b-30
- Red blood cell count
  Dataset: ieu-a-275
- Lymphocyte percentage of white cells
  Dataset: ebi-a-GCST004632
- Platelet count
  Dataset: ieu-a-1008
- Hemoglobin concentration
  Dataset: ebi-a-GCST90002311
- Mean corpuscular hemoglobin concentration
  Dataset: ebi-a-GCST90002329
- Mean corpuscular volume
  Dataset: ebi-a-GCST90002335

https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST004632/

Neutrophil percentage of white cells
Dataset: ebi-a-GCST004633

Eosinophil percentage of white cells

Dataset: ebi-a-GCST004600

**Lymphocyte percentage**
**Dataset: ukb-d-30180_irnt**



**Lymphocyte count**
**Dataset: ukb-d-30120_irnt**

**Lymphocyte count**
**Dataset: bbj-a-36**

**Lymphocyte counts**
**Dataset: ebi-a-GCST004627**





## 1.bloodtrait gwas处理

首先从 UCSC 下载纯文本格式的 dbSNP release 151 并解压，这里下载的是 hg19 版本：

注意，下载的文件的 `chromStart` 列是 0-based 的。0-based 指的是染色体坐标从 0 开始，第一个位置记为 0，而 1-based 则是从 1 开始算出。

```shell
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp151.txt.gz
gzip snp151.txt.gz -d
/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/snp151Common.txt.gz
gzip snp151Common.txt.gz -d

cd /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits
wget https://gwas.mrcieu.ac.uk/files/ieu-b-34/ieu-b-34.vcf.gz
wget https://gwas.mrcieu.ac.uk/files/ieu-b-30/ieu-b-30.vcf.gz
wget https://gwas.mrcieu.ac.uk/files/ieu-b-31/ieu-b-31.vcf.gz
wget https://gwas.mrcieu.ac.uk/files/ieu-b-29/ieu-b-29.vcf.gz
wget https://gwas.mrcieu.ac.uk/files/ieu-b-32/ieu-b-32.vcf.gz
wget https://gwas.mrcieu.ac.uk/files/ieu-b-33/ieu-b-33.vcf.gz
wget https://gwas.mrcieu.ac.uk/files/ieu-a-275/ieu-a-275.vcf.gz
https://gwas.mrcieu.ac.uk/files/ebi-a-GCST004601/ebi-a-GCST004601.vcf.gz

```

```sh

#sed -n ‘253,200p’ filename 
awk -F ' ' '{NF=5}1' /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/snp151.txt >  /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/snp2po.txt

#cat snp151Common.txt | grep '[[:blank:]]chr1[[:blank:]]10019[[:blank:]]' | awk {'print $2 $3 $4 $5 $8 $9'} > snp151.txt

#!/usr/bin/sh
#PBS -N read_ad
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
```



```R

library(bigreadr)
library(stringr)
snp2po <- bigreadr::fread2("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/snp151Common.txt")

snp2po<-snp2po[,c(2,3,5)]
snp2po<-unique(snp2po)
snp2po$rs_number<-paste0(snp2po[,1],snp2po[,2],collapse=":")
snp2po<-snp2po[,-c(1,2)]
save(snp2po,file="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/snp2po.RData")


#!/usr/bin/sh
#PBS -N bloodtraits
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe

#load("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/snp2po.RData")
library(bigreadr)
library(stringr)
for(j in c("MON","WBC","NEU","RDW","RBC","MCV","HCT","LYM","EOS","PLT","MPV","MCH","MCHC","BAS")){
 gwas_data <- bigreadr::fread2(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/BCX2_",j,"_EA_GWAMA.out.gz"))
gwas_data<-gwas_data[,c("rs_number","reference_allele","other_allele","eaf","beta","se","p-value")]
gwas_data<-gwas_data[gwas_data$eaf>0.1,]
gwas_data$chrom <- unlist(lapply(gwas_data$rs_number,function(x){
    str_split(x,pattern=":")[[1]][1]
}))
gwas_data$rs_number <- unlist(lapply(gwas_data$rs_number,function(x){
    str_split(x,pattern="_")[[1]][1]
}))

gwas_data<-gwas_data[,c("rs_number","chrom","eaf","beta","se","p-value")] 

write.table(gwas_data,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/BCX2_",j,"_EA_GWAMA.txt"),quote=F,row.names=F)    
}

```



## 2.基于covid Normal数据

### 2.1血小板PLT

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/")
#读取数据
BCX_PLT<-read.table("/share/pub/dengcy/GWAS_Multiomics/compare/BCX_exomechip_PLT_ALL_20160420.txt",header=T)

#单细胞数据
#"/share2/pub/jiangdp/jiangdp/COVID/data/Normal.rds"

gwasPre<-function(BCX_RBC=BCX_PLT,x="PLT"){
  BCX_RBC<-BCX_RBC[,c("CHR","POS","SNPNAME","Pvalue","se","beta")]
colnames(BCX_RBC)<-c("chrom","pos","rsid","p","se","beta")
##这里不对maf进行过滤，maf信息太少
BCX_RBC$maf<-0.1
#去掉NA
BCX_RBC<-na.omit(BCX_RBC)
write.table(BCX_RBC,file=paste0("BCX_",x,"_exomechip_Processed.txt"),quote=F,row.names=F)  
}
#处理gwas数据
gwasPre(BCX_RBC=BCX_PLT,x="PLT")

##计算scPagwas
 library(scPagwas)
 suppressMessages(library(Seurat))
 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)

 #1.start to run the wrapper functions for preprogress.
 Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share/pub/dengcy/GWAS_Multiomics/compare/BCX_PLT_exomechip_Processed.txt",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data ="/share2/pub/jiangdp/jiangdp/COVID/data/Normal.rds",
                     nfeatures =NULL,
                     Pathway_list=genes.by.reactome.pathway,
                     chrom_ld = chrom_ld)
 Pagwas<-Celltype_heritability_contributions(Pagwas,iters = 200)
 save(Pagwas,file="pbmc1_PLT_scPagwas.RData")
 
 Bootstrap_P_Barplot(Pagwas=Pagwas,
                    figurenames = "/share/pub/dengcy/GWAS_Multiomics/compare/pbmc1_PLT_scPagwas_barplot.pdf",
                    width = 5,
                    height = 7,
                    do_plot=F,
                    title = "PLT pbmc1 scPagwas")

```

![image](92B97492FC2A49818D8B5FE00D9A1421)
这个结果初步符合我们的预期，接下来进行单细胞的计算：

#### 单细胞计算：


```R
library(scPagwas)
suppressMessages(library(Seurat))
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/")
load("pbmc1_PLT_scPagwas.RData")

Pagwas<-Singlecell_heritability_contributions(Pagwas)
save(Pagwas,file="pbmc1_PLT_scPagwas.RData")
##血小板细胞实在太少了
cor_raw<-rep(0,length(Pagwas$Celltype_anno))
cor_raw[Pagwas$Celltype_anno$annotation=="Platelet"]<-1
auc_test(Pagwas$scPagwas_score,cor_raw)
#Area under the curve: 0.6673
```

#### scRDS计算auc

```R
library(bigreadr)
ds<-c("PLT","LYM","MON")
i<-"PLT"
#for(i in ds){
gwas<-bigreadr::fread2(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/BCX_exomechip_",i,"_ALL_20160420.txt"))
#文件中没有N
gwas$N<-100
print(head(gwas))
magma_Input1<-gwas[,c("SNPNAME","Pvalue","N")]
magma_Input2<-gwas[,c("SNPNAME","CHR","POS","POS")]
write.table(magma_Input2,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/magma_",i,"Input2.txt"),sep="\t",row.names=F,quote=F,col.names=F)
colnames(magma_Input1)<-c("SNP","P","N")
write.table(magma_Input1,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/magma_",i,"Input1.txt"),sep="\t",row.names=F,quote=F)
#}
##############
cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
./magma --annotate window=10,10 --snp-loc /share/pub/dengcy/GWAS_Multiomics/compare/magma_PLTInput2.txt \
--gene-loc /share/pub/dengcy/Singlecell/COVID19/MAGMA/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded \
--out /share/pub/dengcy/GWAS_Multiomics/compare/PLTannotated_10kbup_10_down


./magma --bfile /share/pub/dengcy/Singlecell/COVID19/MAGMA/g1000_eur/g1000_eur \
--pval /share/pub/dengcy/GWAS_Multiomics/compare/magma_PLTInput1.txt ncol=3 \
--gene-annot /share/pub/dengcy/GWAS_Multiomics/compare/PLTannotated_10kbup_10_down.genes.annot \
--out /share/pub/dengcy/GWAS_Multiomics/compare/PLTannotated_10kbup_10down

########################
library(testSctpa)
library(dplyr)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare")
PLT_magma<-read.table("PLTannotated_10kbup_10down.genes.out",header=T)
PLT_magma<-PLT_magma[order(PLT_magma$P,decreasing=F),]

Normal_pbmc =readRDS("/share2/pub/jiangdp/jiangdp/COVID/data/Normal.rds")

library(org.Hs.eg.db)
library(data.table)
library(Seurat)

a<-data.frame(gene_id=PLT_magma$GENE)
g2s=toTable(org.Hs.egSYMBOL)
rownames(PLT_magma)<-PLT_magma$GENE
rownames(g2s)<-g2s$gene_id
PLT_magma<-PLT_magma[intersect(PLT_magma$GENE,g2s$gene_id),]
g2s<-g2s[PLT_magma$GENE,]
PLT_magma$symbol<-g2s$symbol
PLT_magma<-na.omit(PLT_magma)
write.table(paste0(PLT_magma$symbol[1:500],collapse = ","),file="PLT_magma_top500genes.txt")

gset <- c("PLT_magma","NA",PLT_magma$symbol[1:500])
gset <- gset%>% 
     as.data.frame() %>% 
    t()
write.table(gset,file = "PLT_magma.gmt",sep = "\t",row.names = F,col.names = F,quote = F)

auccell_PLT = cal_PAS(seurat_object = Normal_pbmc,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
gmt_file="/share/pub/dengcy/GWAS_Multiomics/compare/PLT_magma.gmt")

ssGSEA_PLT = cal_PAS(seurat_object = Normal_pbmc,
                       tool = 'ssGSEA',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
gmt_file="/share/pub/dengcy/GWAS_Multiomics/compare/PLT_magma.gmt")

save(auccell_PLT,file="auccell_PLT.RData")

auccell_df<- as.data.frame(t(GetAssayData(auccell_PLT, slot="data", assay="PAS")))

auccell_auc<-lapply(auccell_df,function(y){
  roc(predictor=y,response=Simulat_5src2$raw_index)$auc
  })

###################
setwd("/share/pub/dengcy/GWAS_Multiomics/compare")
PLT_magma<-read.table("PLTannotated_10kbup_10down.genes.out",header=T)
PLT_magma<-PLT_magma[order(PLT_magma$P,decreasing=F),]
####把单细胞转为python格式数据
setwd("/share/pub/dengcy/GWAS_Multiomics/compare")
library(SeuratDisk)
library(Seurat)
##没有这一步后面会报错
Normal_pbmc =readRDS("/share2/pub/jiangdp/jiangdp/COVID/data/Normal.rds")
DefaultAssay(Normal_pbmc) <- "RNA"
SaveH5Seurat(Normal_pbmc, "Normal_pbmc.h5seurat")
Convert("Normal_pbmc.h5seurat", dest="h5ad")

write.csv(Normal_pbmc@meta.data,file="Normal_pbmc.metadata.csv")
counts<- as.data.frame(GetAssayData(Normal_pbmc,assay = "RNA", slot = "counts"))
counts$gene<-rownames(counts)
write.csv(counts,file="Normal_pbmc.counts.csv",row.names =F)

```

scRDS

```python
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
df_meta = pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Normal_pbmc.metadata.csv",  sep=",")
df_meta = df_meta.rename(
    index={"Unnamed: 0": "cell_id"}
)
df_meta = df_meta.set_index("Unnamed: 0")
df_expr = pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Normal_pbmc.counts.csv",  sep=",",index_col="gene")

df_expr = df_expr.T
df_expr =df_expr.astype(float)
raw_adata = AnnData(df_expr, obs=df_meta)
# assemble AnnData

sc.pp.normalize_total(raw_adata, target_sum=1e4)
sc.pp.log1p(raw_adata)
sc.pp.highly_variable_genes(raw_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

raw_adata = raw_adata[:, raw_adata.var.highly_variable]

sc.pp.scale(raw_adata, max_value=10)
sc.tl.pca(raw_adata, svd_solver="arpack")
sc.pp.neighbors(raw_adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(raw_adata)
sc.tl.leiden(raw_adata)

adata = raw_adata.raw.to_adata()
adata.obsp = raw_adata.obsp
adata.X = adata.X
# save to files
#adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
#2, deleting the backed up raw information
#del adata.raw
adata.write_h5ad("/share/pub/dengcy/GWAS_Multiomics/test/modeldata/Normal_pbmc.h5ad")

#adata = sc.read_h5ad("/share/pub/dengcy/GWAS_Multiomics/compare/Normal_pbmc.h5ad")

df_gs = pd.read_excel("/share/pub/dengcy/GWAS_Multiomics/compare/supp_tables.xlsx", sheet_name="MAGMA gene sets", index_col=0)
df_gs = df_gs.loc[["PLT_magma"], :]
#df_gs = scdrs.util.convert_gs_species(df_gs)
df_gs = df_gs.rename(
    index={"PLT_magma": "PLT"}
)
df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/geneset.gs", sep="\t", index=False)
# compile covariates
df_cov = pd.DataFrame(index=adata.obs.index)
#df_cov["index"] = adata.obs.index
df_cov["const"] = 1
df_cov.to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/cov.tsv", sep="\t")
#df_cov["pop"] = adata.obs["pop"]

python /share/pub/dengcy/GWAS_Multiomics/test/modeldata/scRDS/compute_score.py \
    --h5ad_file /share/pub/dengcy/GWAS_Multiomics/compare/Normal_pbmc.h5ad \
    --h5ad_species human \
    --gs_file /share/pub/dengcy/GWAS_Multiomics/compare/geneset.gs \
    --gs_species human \
    --cov_file /share/pub/dengcy/GWAS_Multiomics/compare/cov.tsv \
    --flag_filter True \
    --flag_raw_count True \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score True \
    --out_folder /share/pub/dengcy/GWAS_Multiomics/compare/scRDS/

###############################

```



### subfunctions:


```
auc_test<-function(scPagwas_score,cor_raw){
  if(Inf %in% scPagwas_score){
  scPagwas_score[which(scPagwas_score==Inf)]<- max(scPagwas_score[-which(scPagwas_score==Inf)])
  }
    if(-Inf %in% scPagwas_score){
  scPagwas_score[which(scPagwas_score==-Inf)]<- min(scPagwas_score[-which(scPagwas_score==-Inf)])
  }
  AUC <- roc(predictor=scPagwas_score,response=cor_raw)
  return(AUC)
}

gwasPre<-function(BCX_RBC=BCX_PLT,x="PLT"){
  BCX_RBC<-BCX_RBC[,c("CHROM","POS","SNPNAME","Pvalue","se","beta")]
colnames(BCX_RBC)<-c("chrom","pos","rsid","p","se","beta")
##这里不对maf进行过滤，maf信息太少
BCX_RBC$maf<-0.1
#去掉NA
BCX_RBC<-na.omit(BCX_RBC)
write.table(BCX_RBC,file=paste0("BCX_",x,"_exomechip_Processed.txt"),quote=F,row.names=F)  
}
```

### 2.2 monocyte count (Mono)

"BCX_exomechip_MON_ALL_20160420.txt"


```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/")
#读取数据
BCX_MON<-read.table("/share/pub/dengcy/GWAS_Multiomics/compare/BCX_exomechip_MON_ALL_20160420.txt",header=T)

#处理gwas数据
gwasPre(BCX_RBC=BCX_MON,x="MON")

 library(scPagwas)
 suppressMessages(library(Seurat))
 data(Genes_by_pathway_kegg)
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
 data(block_annotation)
 data(chrom_ld)

 #1.start to run the wrapper functions for preprogress.
 load("pbmc1_PLT_scPagwas.RData")

   suppressMessages(gwas_data <- bigreadr::fread2("/share/pub/dengcy/GWAS_Multiomics/compare/BCX_MON_exomechip_Processed.txt"))
   Pagwas <- GWAS_summary_input(Pagwas=Pagwas,
                                  gwas_data=gwas_data)
    snp_gene_df<-Snp2Gene(snp=Pagwas$gwas_data,refGene=block_annotation)
    snp_gene_df$slope <- rep(1,nrow(snp_gene_df))
    Pagwas$snp_gene_df<-snp_gene_df[snp_gene_df$Disstance=="0",]
   Pagwas <- Pathway_annotation_input(Pagwas=Pagwas,
                                       block_annotation=block_annotation)
    Pagwas <- Link_pathway_blocks_gwas(Pagwas=Pagwas,
                                       chrom_ld=chrom_ld)
 
 Pagwas<-Celltype_heritability_contributions(Pagwas=Pagwas,iters = 200)
 save(Pagwas,file="pbmc1_MON_scPagwas.RData")
 
 Bootstrap_P_Barplot(Pagwas=Pagwas,
                    figurenames = "/share/pub/dengcy/GWAS_Multiomics/compare/pbmc1_MON_scPagwas_barplot.pdf",
                    width = 5,
                    height = 7,
                    do_plot=F,
                    title = "MON pbmc1 scPagwas")
```
![image](E1F880D15FD6425FA2A2D5938A32B5C7)
竟然没有CD16

单细胞数据计算：
<font size=6>

```R
library(dplyr)
require(pROC)
require(ggplot2)
Pagwas<-Singlecell_heritability_contributions(Pagwas)
 save(Pagwas,file="pbmc1_MON_scPagwas.RData")
 load("pbmc1_MON_scPagwas.RData")
 
Celltype_anno<- Pagwas$Celltype_anno[Pagwas$Celltype_anno$annotation %in% c("CD14+monocyte","NK","Platelet"),]
scPagwas_score<-Pagwas$scPagwas_score[Celltype_anno$cellnames]
 
cor_raw<-rep(0,dim(Celltype_anno)[1])

cor_raw[Celltype_anno$annotation %in% c("CD14+monocyte")]<-1
auc_test(scPagwas_score,cor_raw)
```

</font>

```R

```



## pbmc4k数据进行测试


```R
 library(scPagwas)
 suppressMessages(library(Seurat))
 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)

 #1.start to run the wrapper functions for preprogress.
 Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share/pub/dengcy/GWAS_Multiomics/compare/BCX_MON_exomechip_Processed.txt",
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seurat_pbmc4k.rds",
                     nfeatures =NULL,
                     Pathway_list=genes.by.reactome.pathway,
                     chrom_ld = chrom_ld)
 Pagwas<-Celltype_heritability_contributions(Pagwas,iters = 200)
 save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/compare/pbmc4k_MON_scPagwas.RData")
 
 Bootstrap_P_Barplot(Pagwas=Pagwas,
                    figurenames = "/share/pub/dengcy/GWAS_Multiomics/compare/pbmc4k_MON_scPagwas_barplot.pdf",
                    width = 5,
                    height = 7,
                    do_plot=F,
                    title = "MON_ pbmc4k scPagwas")



```

## PBMC胡桓得到的数据以及骨髓数据

| HSC - Hematopoietic Stem Cells  造血干细胞                   |
| ------------------------------------------------------------ |
| Early Eryth. - Early Erythroid Cells 早期红细胞              |
| Late Eryth. - Late Erythroid Cells 晚期红细胞                |
| Early Basophil Cells 早期嗜碱细胞                            |
| CMP/LMPP - Common Myeloid Progenitor / lymphoid-primed multipotent progenitors：常见的髓系祖/淋巴组织诱导的多能祖细胞 |
| CLP 1 - Common Lymphoid Progenitor 常见的淋巴祖              |
| GMP - Granulocyte-Monocyte progenitors 粒细胞单核细胞祖细胞  |
| GMP/Neut - Granulocyte-Monocyte progenitors / Neutrophil  粒细胞单核细胞祖细胞/嗜中性粒细胞 |
| pDC - Plasmacytoid Dendritic Cell 浆细胞样树突状细胞         |
| cDC - Classical Dendritic CElls 树突状细胞                   |
| CD14 Mono 1 - CD14+ Monocyte Cells                           |
| CD14 Mono 2 - CD14+ Monocyte Cells                           |
| CD16 Mono - CD16+ Monocyte Cells                             |
| Unk - Unkown                                                 |
| CLP 2 - Common Lymphoid Progenitor  常见的淋巴祖             |
| Pre B - Pre B-Cell Progenitor B细胞前体祖细胞                |
| B - B Cells                                                  |
| Plasma - Plasma Cells                                        |
| CD8 N - CD8+ Naïve Cells                                     |
| CD4 N 1 - CD4+ Naïve Cells                                   |
| CD4 N 2 - CD4+ Naïve Cells                                   |
| CD4 M - CD4+ Memory Cells                                    |
| CD8 EM - CD8+ Effector Memory Cells                          |
| CD8 CM - CD8+ Central Memory Cells                           |
| NK - Natural Killer Cells                                    |
| Unk - Unkown                                                 |

### 1.前期准备数据

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

###############计算差异基因
parallelFindAllMarkers <- function(obj){

  all_markers <- future_lapply(levels(Idents(obj)), function(x){ # Five expression patterns
    FindMarkers(obj, ident.1 = x, ident.2 = NULL, test.use = "MAST")
  })

  return(value(all_markers))
}

library(future)
library(future.apply)
plan("multiprocess", workers = 3)
  Allmarkers <- parallelFindAllMarkers(NM_Healthy_pbmc)
 save(Allmarkers,file="All_cluster_markers.RData")
#Allmarkers[[0+15]]<-NULL
Allmarkers<-lapply(0:(length(Allmarkers)-1),function(x){
  df<-Allmarkers[[x+1]]
  df<-df[order(-df$avg_log2FC),]
  df$cluster<-rep(x,nrow(df))
  return(df)
  })
Allmarkers<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),Allmarkers)
write.csv(Allmarkers,file="17cluster_Allmarkers.csv")
```



```R
######################################################
pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")
#i<-pbmccells[1]
setwd("/share/pub/dengcy/GWAS_Multiomics/compare")
 library(scPagwas)
 suppressMessages(library(Seurat))
 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 #load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld) 
Single_data=readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_pbmc.rds")
#NM_Healthy_data_scexpr <- Seurat::AverageExpression(Single_data,assays="RNA")[["RNA"]]
#save(NM_Healthy_data_scexpr,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_data_scexpr.RData")


    Pagwas <- Single_data_input(Pagwas=NULL,
                                assay="RNA",
                                Single_data=Single_data,
                                Pathway_list=Genes_by_pathway_kegg)
  Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,
                                 Pathway_list=Genes_by_pathway_kegg
                                 )
save(Pagwas ,file="NM_Healthy_pbmc_prePagwas.RData")

######骨髓单细胞数据
setwd("/share/pub/dengcy/GWAS_Multiomics/compare")
 library(scPagwas)
 suppressMessages(library(Seurat))
 #Input pathway gene list, you can construct with youself.
 #data(Genes_by_pathway_kegg)
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld) 
#Seu_Hema_data=readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Healthy_Hema.rds")
#Seu_Hema_data<-NormalizeData(Seu_Hema_data)
#Seu_Hema_data<-ScaleData(Seu_Hema_data)
#saveRDS(Seu_Hema_data,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.RData")
#Seu_Hema_data_scexpr <- Seurat::AverageExpression(Seu_Hema_data,assays="RNA")[["RNA"]]
#save(Seu_Hema_data_scexpr,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data_scexpr.RData")


Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")
Pagwas <- Single_data_input(Pagwas=NULL,
                                assay="RNA",
                                Single_data=Single_data,
                                Pathway_list=genes.by.reactome.pathway)
Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,
                                 Pathway_list=genes.by.reactome.pathway
                                 )
save(Pagwas ,file="Seu_Healthy_Hema_prePagwas.RData")

###########骨髓kegg
setwd("/share/pub/dengcy/GWAS_Multiomics/compare")
library(scPagwas)
suppressMessages(library(Seurat))
 #Input pathway gene list, you can construct with youself.
data(Genes_by_pathway_kegg)
 #load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
 #gene annotation files.
data(block_annotation)
 #LD data
data(chrom_ld) 
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")
  Pagwas <- Single_data_input(Pagwas=NULL,
                                assay="RNA",
                                Single_data=Single_data,
                                Pathway_list=Genes_by_pathway_kegg)
  Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,
                                 Pathway_list=Genes_by_pathway_kegg
                                 )
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Seu_Healthy_Hema_kegg_prePagwas.RData")


单细胞数据平均表达数据地址："/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data_scexpr.RData"
GWAS数据地址："/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/"
traits:"RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume"
gwas数据例子：RBCcount_prune_gwas_data.txt
```



### 2.接下来计算scPagwas

```R
pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

setwd("/share/pub/dengcy/GWAS_Multiomics/compare")
load("NM_Healthy_pbmc_prePagwas.RData")
 library(scPagwas)
 suppressMessages(library(Seurat))
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
 data(block_annotation)
 data(chrom_ld) 

for(i in pbmccells[c(3,4,7)]){
    print(i)
 Pagwas2<-scPagwas_main(Pagwas =Pagwas,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"),
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
 Pagwas2<-Celltype_heritability_contributions(Pagwas2,iters = 200)
 save(Pagwas2,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/",i,"_NM_pbmc_scPagwas.RData"))
 
 Bootstrap_P_Barplot(Pagwas=Pagwas2,
                    figurenames = paste0("/share/pub/dengcy/GWAS_Multiomics/compare/NM_pbmc_",i,"_scPagwas_barplot.pdf"),
                    width = 5,
                    height = 7,
                    do_plot=F,
                    title = paste0(i,"_NM_Healthy pbmc scPagwas"))
}


pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")
setwd("/share/pub/dengcy/GWAS_Multiomics/compare")
load("/share/pub/dengcy/GWAS_Multiomics/compare/Seu_Healthy_Hema_prePagwas.RData")

 library(scPagwas)
 suppressMessages(library(Seurat))
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
 data(block_annotation)
 data(chrom_ld) 

for(i in pbmccells[1:4]){
#i<-pbmccells[1]
   print(i)
 Pagwas2<-scPagwas_main(Pagwas =Pagwas,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"),
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
 Pagwas2<-Celltype_heritability_contributions(Pagwas2,iters = 200)
 save(Pagwas2,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/",i,"_Seu_Healthy_scPagwas.RData"))
 
 Bootstrap_P_Barplot(Pagwas=Pagwas2,
                    figurenames = paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Seu_Healthy_",i,"_scPagwas_barplot.pdf"),
                    width = 5,
                    height = 7,
                    do_plot=F,
                    title = paste0(i,"_Seu_Healthy pbmc scPagwas"))
}
#########################################
```

#### 重新计算PBMC的结果 v1.7.1

这里先只计算RBCcount得原始gwas

```R

#"Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc")
library(scPagwas)
suppressMessages(library(Seurat))
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)

load("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/NM_Healthy_pbmc_prePagwas.RData")
#i<-"RBCcount"
i<-"Plateletcount"
     Pagwas<-scPagwas_main(Pagwas =Pagwas,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_gwas_data.txt"),
                           Single_data ="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_pbmc.rds",
                     output.prefix=i,
                     singlecell=F,
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/",i,"_NM_pbmc_scPagwas.RData"))
```

#### 重新计算BMMC的结果 v1.6.1

```R
pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
library(scPagwas)
suppressMessages(library(Seurat))
load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Seu_Healthy_Hema_prePagwas.RData")
i<-pbmccells[12]
     Pagwas<-scPagwas_main(Pagwas =Pagwas,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"),
                     output.prefix=i,
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.6.RData"))
```

#### 重新计算BMMC的结果1.7.1

```R
pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
library(scPagwas)
suppressMessages(library(Seurat))
load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Seu_Healthy_Hema_prePagwas.RData")
i<-pbmccells[9]
     Pagwas<-scPagwas_main(Pagwas =Pagwas,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"),
                     output.prefix=i,
                     ncores=5,
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.7.RData"))
  
source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/3.r


```

原始gwas(最后只计算得到)

```R
  
source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/12.r


#############################################reactome

Single_data = readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")
Single_data1<-subset(Single_data,)
Single_datacluster1 <- subset(x = Single_data, idents = c("05_CMP.LMPP","08_GMP.Neut"     ,"01_HSC","06_CLP.1","15_CLP.2","02_Early.Eryth","07_GMP","09_pDC","04_Early.Baso"  ,"03_Late.Eryth"))
Single_datacluster2 <- subset(x = Single_data, idents = c("17_B","12_CD14.Mono.2"     ,"16_Pre.B","10_cDC","11_CD14.Mono.1","25_NK"))

Single_datacluster3 <- subset(x = Single_data, idents = c("21_CD4.N2","22_CD4.M","23_CD8.EM","19_CD8.N","24_CD8.CM","26_Unk","20_CD4.N1","14_Unk","13_CD16.Mono","18_Plasma"))
saveRDS(Single_datacluster1,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data1.rds")
saveRDS(Single_datacluster2,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data2.rds")
saveRDS(Single_datacluster3,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data3.rds")

pbmccells<-c("lymphocytecount","LymphocytePercent","Lymphocytepercentage2")
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
library(scPagwas)
data(Genes_by_pathway_kegg)
 data(block_annotation)
 data(chrom_ld) 

suppressMessages(library(Seurat))
i<-pbmccells[3]

Pagwas<-scPagwas_main(Pagwas =NULL,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_gwas_data.txt"),
                     Single_data = "/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data1.rds",
                     output.prefix=i,
                       assay="RNA",
                      output.dirs=paste0(i,"1"),
                      Pathway_list=Genes_by_pathway_kegg,
                     ncores=10,
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/raw1",i,"_Hema_bmmc_scPagwas_v1.7.3.RData"))

Pagwas<-scPagwas_main(Pagwas =NULL,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_gwas_data.txt"),
                     Single_data = "/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data2.rds",
                     output.prefix=i,
                       assay="RNA",
                      output.dirs=paste0(i,"2"),
                      Pathway_list=Genes_by_pathway_kegg,
                     ncores=10,
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/raw2",i,"_Hema_bmmc_scPagwas_v1.7.3.RData"))

Pagwas<-scPagwas_main(Pagwas =NULL,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_gwas_data.txt"),
                     Single_data = "/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data3.rds",
                     output.prefix=i,
                       assay="RNA",
                      output.dirs=paste0(i,"3"),
                      Pathway_list=Genes_by_pathway_kegg,
                     ncores=10,
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/raw3",i,"_Hema_bmmc_scPagwas_v1.7.3.RData"))
```

export OPENBLAS_NUM_THREADS=1

Lymphocytecount2计算

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
library(scPagwas)
data(Genes_by_pathway_kegg)
 data(block_annotation)
 data(chrom_ld) 
suppressMessages(library(Seurat))
Pagwas<-scPagwas_main(Pagwas =NULL,
                     gwas_data ="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/Lymphocytecount2_gwas_data.txt",
                     Single_data = "/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds",
                     output.prefix="Lymphocytecount2",
                       assay="RNA",
                      output.dirs="Lymphocytecount2",
                      Pathway_list=Genes_by_pathway_kegg,
                     ncores=10,
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount2_Hema_bmmc_scPagwas_v1.9.RData")


```

Lymphocytecount3计算

```
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
library(scPagwas)
data(Genes_by_pathway_kegg)
 data(block_annotation)
 data(chrom_ld) 

suppressMessages(library(Seurat))

Pagwas<-scPagwas_main(Pagwas =NULL,
                     gwas_data ="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/Lymphocytecount3_gwas_data.txt",
                     Single_data = "/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds",
                     output.prefix="Lymphocytecount3",
                       assay="RNA",
                      output.dirs="Lymphocytecount3",
                      Pathway_list=Genes_by_pathway_kegg,
                     ncores=10,
                     add_eqtls="OnlyTSS",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas_v1.7.3.RData")
```

#### 重新计算BMMC的结果1.9.1

```R
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
suppressMessages(library(Seurat))
load("/share/pub/dengcy/GWAS_Multiomics/compare/Seu_Healthy_Hema_kegg_prePagwas.RData")
i<-"Lymphocytecount3"

pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

Pagwas<-scPagwas_main(Pagwas =Pagwas,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_gwas_data.txt"),
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds",
                     output.prefix=i,
                           seruat_return=F,
                     Pathway_list=Genes_by_pathway_kegg,
                     output.dirs=paste0(i,"_scPagwas"),
                     ncores=2,
                      assay="RNA",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
 
#!/usr/bin/sh
#PBS -N monocyte
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=5
#PBS -j oe
source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scPagwas1.8/1.r

```

```R
pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scPagwas1.8")
library(scPagwas)
suppressMessages(library(Seurat))
load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Seu_Healthy_Hema_prePagwas.RData")
 data(Genes_by_pathway_kegg)
 data(block_annotation)
 data(chrom_ld)
i<-"Lymphocytecount3"
     Pagwas<-scPagwas_main(Pagwas =NULL,
                     gwas_data ="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/Lymphocytecount3_gwas_data.txt",
                     Pathway_list=Genes_by_pathway_kegg,             Single_data="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds",
                     output.dirs=paste0(i,"_scPagwasv1.9"),
                     output.prefix=i,
                     assay="RNA",
                     celltype = T,
                     seruat_return=F,
                     ncores=10,
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas_v1.7.RData")
 

i<-"RBCcount"
     Pagwas<-scPagwas_main(Pagwas =NULL,
                     gwas_data ="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/RBCcount_gwas_data.txt",
                     Pathway_list=Genes_by_pathway_kegg,             Single_data="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds",
                     output.dirs=paste0(i,"_scPagwasv1.9"),
                     output.prefix="RBCcount",
                     assay="RNA",
                     ncores=2,
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scPagwas1.8/RBCcount_Hema_bmmc_scPagwas_v1.9.RData")

source activate R4
export OPENBLAS_NUM_THREADS=1
Rscript /share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/3.r
```

##### 新增"RBCcount2","Lymphocytepercentage2"

```R
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
pbmccells<-c("RBCcount2","Lymphocytepercentage2")

 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)

suppressMessages(library(Seurat))
load("/share/pub/dengcy/GWAS_Multiomics/compare/Seu_Healthy_Hema_kegg_prePagwas.RData")

pbmccells<-c("RBCcount2","Lymphocytepercentage2")
i<-pbmccells[1]
Pagwas<-scPagwas_main(Pagwas =Pagwas,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"),
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds",
                     output.prefix=i,
                           seruat_return=F,
                     Pathway_list=Genes_by_pathway_kegg,
                     output.dirs=paste0(i,"_scPagwas"),
                     ncores=2,
                      assay="RNA",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
   # }
load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytepercentage2_Hema_bmmc_scPagwas_v1.9.1.RData")
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/")
scPagwas_Visualization(Single_data=Pagwas,
                                   p_thre = 0.05,
                                   FigureType = "tsne",
                                   width = 7,
                                   height = 7,
                                   lowColor = "#000957", 
                        highColor = "#EBE645",
                        output.dirs=paste0(i,"_scPagwas"),
                                   size = 0.5,
                                   do_plot = F)

```



#### 重新计算PBMC的结果v1.9.1

/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc

```R
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)

suppressMessages(library(Seurat))
load("/share/pub/dengcy/GWAS_Multiomics/compare/NM_Healthy_pbmc_prePagwas.RData")
#################
#跑细胞类型
#################
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/")
#for(i in pbmccells[11:12]){
    i<-pbmccells[13]
i<-"Lymphocytecount3"
     Pagwas<-scPagwas_main(Pagwas =Pagwas,
                     gwas_data =paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_gwas_data.txt"),
                     #Single_data ="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_pbmc.rds",
                     output.prefix=i,
                           singlecell=F,
                           seruat_return=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     output.dirs=paste0(i,"_pbmc_scPagwasv1.9"),
                     ncores=5,
                     assay="RNA",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/",i,"_pbmc_scPagwas.RData"))
   # }

#################
#跑单细胞
##################
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume","lymphocytecount3","lymphocytecount2")

 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc")
suppressMessages(library(Seurat))
load("/share/pub/dengcy/GWAS_Multiomics/compare/NM_Healthy_pbmc_prePagwas.RData")

i<-"RBCcount"
i<-"monocytecount"
i<-"lymphocytecount3"
i<-"lymphocytecount2"
i<-"Hemoglobinconcen"
i<-"MeanCorpusVolume"

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
```

### 3.画森林图

```R
pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

 library(scPagwas)
 suppressMessages(library(Seurat))

for(i in pbmccells){
#i<-pbmccells[1]
   print(i)
    load(paste0(i,"_Seu_Healthy_scPagwas.RData"))
    if("Pagwas2" %in% ls()){
        pdf( paste0(i,"_Seu_Healthy_scPagwas.pdf"))
     Bootstrap_estimate_Plot(Pagwas=Pagwas2,
                        figurenames =NULL ,
                        width = 9,
                        height = 7,
                        do_plot=T)   
        dev.off()
    }else{
         pdf( paste0(i,"_Seu_Healthy_scPagwas.pdf"))
        Bootstrap_estimate_Plot(Pagwas=Pagwas,
                        figurenames = paste0(i,"_Seu_Healthy_scPagwas.pdf"),
                        width = 9,
                        height = 7,
                        do_plot=T)  
         dev.off()
    }

    }

```

### 4.骨髓数据热图

#### 1.输出bmmc细胞类型p值结果

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
#"RBCcount","Plateletcount",
pbmccells<-c("basophilcount","eosinophilcount" ,"Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")
i<-pbmccells[1]
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))

result_list<-lapply(pbmccells,function(i){
   print(i)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))

    a<-Pagwas@misc$bootstrap_results[-1,]
    return(a$bp_value)
})
names(result_list)<-pbmccells
result_list<-as.data.frame(result_list)
#load("RBCcount_Hema_bmmc_scPagwas_v1.7.RData")
rownames(result_list)<-rownames(Pagwas@misc$bootstrap_results[-1,])
save(result_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/scPagwas_bmmc_list.RData")
write.csv(result_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/scPagwas_bmmc_list.csv")

setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc")
pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume","Lymphocytecount3","Lymphocytecount2")

result_list<-lapply(pbmccells,function(i){
   print(i)
    load(paste0(i,"_pbmc_scPagwas.RData"))
        a<-Pagwas$bootstrap_results[-1,]
        a<-a$bp_value
    return(a)
})

names(result_list)<-pbmccells
result_list<-as.data.frame(result_list)
load("RBCcount_pbmc_scPagwas.RData")
rownames(result_list)<-rownames(Pagwas$bootstrap_results[-1,])

result_list<-result_list[,1:12]
save(result_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scPagwas_PBMC_celltypes_resultlist.RData")
write.csv(result_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/scPagwas_pbmc_list.csv")

```

## 血液细胞traits用MAGMA进行计算

### 1.预处理

```R
###批量获取magma数据gwas文件
library(bigreadr)
pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

pbmccells<-c("RBCcount","Plateletcount")
for(i in pbmccells){
    gwas<-bigreadr::fread2(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"))
#gwas_AD$QUAL<-100
print(head(gwas))
magma_Input1<-gwas[,c("rsid","p","N")]
magma_Input2<-gwas[,c("rsid","chrom","pos","pos")]
write.table(magma_Input2,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_",i,"Input2.txt"),sep="\t",row.names=F,quote=F,col.names=F)
colnames(magma_Input1)<-c("SNP","P","N")
write.table(magma_Input1,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_",i,"Input1.txt"),sep="\t",row.names=F,quote=F)

}
############"RBCcount"进行原始gwas得计算
pbmccells<-c("RBCcount")
for(i in pbmccells){
    gwas<-bigreadr::fread2(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_gwas_data.txt"))
#gwas_AD$QUAL<-100
print(head(gwas))
magma_Input1<-gwas[,c("rsid","p","N")]
magma_Input2<-gwas[,c("rsid","chrom","pos","pos")]
write.table(magma_Input2,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_",i,"Input2.txt"),sep="\t",row.names=F,quote=F,col.names=F)
colnames(magma_Input1)<-c("SNP","P","N")
write.table(magma_Input1,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_",i,"Input1.txt"),sep="\t",row.names=F,quote=F)

}

```

### 2.注释基因区域的SNP

```sh
#!/usr/bin/sh
#PBS -N magma
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe

cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
for i in RBCcount 
#Plateletcount 
#basophilcount eosinophilcount lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount LymphocytePercent Hemoglobinconcen MeanCorpuscularHemoglobin MeanCorpusVolume
do
./magma --annotate window=10,10 --snp-loc /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_${i}Input2.txt \
--gene-loc /share/pub/dengcy/Singlecell/COVID19/MAGMA/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded \
--out /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10_down
done

```

### 3.计算得到基因水平的关联结果

```sh
cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
#RBCcount Plateletcount basophilcount 
for i in RBCcount 
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

### 4.单细胞数据输入数据计算

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
###1.加载基因坐标，并且将上下游延长50kb
############################
#Filtered to remove extended MHC (chr6, 25Mb to 34Mb).
gene_coordinates <-read_tsv("/share/pub/dengcy/Singlecell/COVID19/MAGMA/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded",
           col_names = FALSE,col_types = 'cciicc') %>%
  mutate(start=ifelse(X3-50000<0,0,X3-50000),end=X4+50000) %>%
  select(X2,start,end,6,1) %>%
  rename(chr="X2", Gene="X6",EntrezGene="X1") %>%
  mutate(chr=paste0("chr",chr))

#########函数计算top10基因，内嵌小函数
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
                #write_tsv(exp,"single_cell_clusterExpression.txt")
                ###3.Specificity Calculation
                #The specifitiy is defined as the proportion of total expression performed by the cell type of interest (x/sum(x)).
                exp<- exp %>%
                group_by(Gene) %>%
                mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean)) %>%
                ungroup()
                ###4.Get MAGMA genes
                #Only keep genes that are tested in MAGMA
                exp2 <- inner_join(exp,gene_coordinates,by="Gene")
                ###5.Get number of genes

                #Get number of genes that represent 10% of the dataset
                n_genes <- length(unique(exp2$EntrezGene))
                n_genes_to_keep <- (n_genes * 0.1) %>% round()
                ###7.Write MAGMA/LDSC input files
                #Filter out genes with expression below 1 TPM.
                #exp3<-exp2
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

#####运行top10_function
top10_function(exp=as.data.frame(cell_expr$Hema))
top10_function(as.data.frame(cell_expr$pbmc))
```

### 5.进行MAGMA的运算

```R
#1.
cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
for i in RBCcount 
#Plateletcount 
#basophilcount eosinophilcount lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount LymphocytePercent Hemoglobinconcen MeanCorpuscularHemoglobin MeanCorpusVolume
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



###############################
cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
for i in RBCcount Plateletcount basophilcount
do
int="/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10down.genes.raw"
cell_type="/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/Hema_top10.txt"
./magma --gene-results  $int --set-annot  $cell_type --out /share/pub/dengcy/GWAS_Multiomics/compare/magma/${i}_magma_Hema
done

cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
for i in eosinophilcount lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount LymphocytePercent Hemoglobinconcen MeanCorpuscularHemoglobin MeanCorpusVolume
do
int="/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10down.genes.raw"
cell_type="/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/Hema_top10.txt"
./magma --gene-results  $int --set-annot  $cell_type --out /share/pub/dengcy/GWAS_Multiomics/compare/magma/${i}_magma_Hema
done


```

### 6.magma结果的整合

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/magma")
ds<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")
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

## 可视化结果：热图

### 1.scPagwas

```R
library(gplots)
library(RColorBrewer)

setwd("E:/OneDrive/GWAS_Multiomics/Compare")
load("scPagwas_bmmc_list.RData")
result_list2<- -log2(result_list)
#magma_pbmc_p

coul <- colorRampPalette(brewer.pal(8, "Oranges"))(25)#换个好看的颜色
#hM <- format(round(result_list2, 2))#对数据保留2位小数
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
          trace="none",#不显示trace
          col=coul,#修改热图颜色
          density.info = "none",#图例取消density
          key.xlab ='Correlation',
          key.title = "",
          cexRow = 1,cexCol = 1,#修改横纵坐标字体
          Rowv = F,Colv =T, #去除聚类
          margins = c(6, 6),
          cellnote = hM,notecol='black'#添加相关系数的值及修改字体颜色
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

coul <- colorRampPalette(brewer.pal(8, "Oranges"))(25)#换个好看的颜色
#hM <- format(round(result_list2, 2))#对数据保留2位小数
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
          trace="none",#不显示trace
          col=coul,#修改热图颜色
          density.info = "none",#图例取消density
          key.xlab ='Correlation',
          key.title = "",
          cexRow = 1,cexCol = 1,#修改横纵坐标字体
          Rowv = F,Colv =T, #去除聚类
          margins = c(6, 6),
          cellnote = hM,notecol='black'#添加相关系数的值及修改字体颜色
)
dev.off()
```



### 3.rolypoly

```R
setwd("/share2/pub/jiangdp/jiangdp/COVID/rolypoly")

pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

rolypoly_pbmc_list<-lapply(pbmccells,function(i){
    message(i)
    load(paste0("roly_NM_",i,".RData"))
    return(rolypoly_result$bootstrap_results$bp_value[-1])
})
names(rolypoly_pbmc_list)<-pbmccells
rolypoly_pbmc_list<-as.data.frame(rolypoly_pbmc_list)
load("roly_NM_RBCcount.RData")
rownames(rolypoly_pbmc_list)<-rownames(rolypoly_result$bootstrap_results[-1,])

save(rolypoly_pbmc_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/rolypoly_pbmc_list.RData")
write.csv(rolypoly_pbmc_list,file="/share/pub/dengcy/GWAS_Multiomics/compare/rolypoly_pbmc_list.csv")

#############BMMC
setwd("/share2/pub/jiangdp/jiangdp/COVID/rolypoly")
rolypoly_bmmc_list<-lapply(pbmccells,function(i){
    message(i)
    load(paste0("roly_Seu_",i,".RData"))
    return(rolypoly_result$bootstrap_results$bp_value[-1])
})
names(rolypoly_bmmc_list)<-pbmccells
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

coul <- colorRampPalette(brewer.pal(8, "Oranges"))(25)#换个好看的颜色
#hM <- format(round(result_list2, 2))#对数据保留2位小数
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
          trace="none",#不显示trace
          col=coul,#修改热图颜色
          density.info = "none",#图例取消density
          key.xlab ='Correlation',
          key.title = "",
          cexRow = 1,cexCol = 1,#修改横纵坐标字体
          Rowv = F,Colv =T, #去除聚类
          margins = c(6, 6),
          cellnote = hM,notecol='black'#添加相关系数的值及修改字体颜色
)
dev.off()
```

### 4.ldsc

```R
setwd("/share2/pub/zhenggw/zhenggw/scPagwas_LDSC/NM_results")
pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")


ldsc_pbmc_result<-lapply(pbmccells,function(i){
    message(i)
   result<-read.table(file=paste0(i,".cell_type_results.txt"),header = T)
    return(result)
})
names(ldsc_pbmc_result)<-pbmccells

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
ldsc_bmmc_result<-lapply(pbmccells,function(i){
    message(i)
   result<-read.table(file=paste0(i,".cell_type_results.txt"),header = T)
    return(result)
})
names(ldsc_bmmc_result)<-pbmccells

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

#########


#################rolypoly
##########################magma
load("ldsc_bmmc_list.RData")
result_list<-ldsc_bmmc_list[,3:12]
result_list2<- -log2(result_list)
#magma_pbmc_p

coul <- colorRampPalette(brewer.pal(8, "Oranges"))(25)#换个好看的颜色
#hM <- format(round(result_list2, 2))#对数据保留2位小数
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
          trace="none",#不显示trace
          col=coul,#修改热图颜色
          density.info = "none",#图例取消density
          key.xlab ='Correlation',
          key.title = "",
          cexRow = 1,cexCol = 1,#修改横纵坐标字体
          Rowv = F,Colv =T, #去除聚类
          margins = c(6, 6),
          cellnote = hM,notecol='black'#添加相关系数的值及修改字体颜色
)
dev.off()
```



## 排秩比较可视化

### 1.计算pbmc单细胞数据结果

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc")
 pbmccells<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")
 library(scPagwas)
i<-pbmccells[1]

load(paste0(i,"_NM_pbmc_scPagwas.RData"))
pdf(paste0("forest",i,".pdf"))
Bootstrap_estimate_Plot(Pagwas=Pagwas,
                        figurenames = NULL,
                        width = 9,
                        height = 7,
                        do_plot=T)
dev.off()
Pagwas<-Singlecell_heritability_contributions(Pagwas)
save(Pagwas,file=paste0(i,"_NM_pbmc_scPagwas.RData"))

```

### 2.计算不同情况的打分v1.7.1

#### 2.1 得到画图的meta data

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
 library(scPagwas)
 suppressMessages(library(Seurat))
library(org.Hs.eg.db)
library(dplyr)

 pbmccells<-c("RBCcount","Plateletcount","eosinophilcount","basophilcount","LymphocytePercent","Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")
#pbmccells<-c("basophilcount","LymphocytePercent")
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")
g2s=toTable(org.Hs.egSYMBOL)

setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS")

for(i in pbmccells){
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
##计算magma得分

#magmatop1000<-intersect(magma_genes$symbol,rownames(Single_data))
    magmatop1000<-magma_genes$symbol
magmatop1000<-magmatop1000[1:1000]
scPagwastop1000<-scPagwas_genes[1:1000]
# magmatop500<-intersect(magma_genes$symbol,rownames(Single_data))
magmatop500<-magmatop1000[1:500]
scPagwastop500<-scPagwas_genes[1:500]
topgene<-list(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500)
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
    
save(topgene,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_",i,".RData"))
    
    }

magmatop1000<-paste(magmatop1000,collapse=",")
scPagwastop1000<-paste(scPagwastop1000,collapse=",")
magmatop500<-paste(magmatop500,collapse=",")
scPagwastop500<-paste(scPagwastop500,collapse=",")

a<-data.frame(genes=c(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500))
rownames(a)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")

#a<-data.frame(genes=c(magmatop1000,magmatop500))
#rownames(a)<-c("magmatop1000","magmatop500")

write.csv(a,file=paste0("scPagwas.magmatopgenes_scRDS",i,".csv"))

meta.data<-data.frame(scPagwas_score=Pagwas$scPagwas_score[colnames(Single_data)],
 CellsrankPvalue=Pagwas$CellsrankPvalue[colnames(Single_data),"pValueHigh"]
                     )
meta.data$positive<-rep("positive",ncol(Single_data))
meta.data$positive[Single_data$CellsrankPvalue>0.05]<-"negative"
save(meta.data,file=paste0(i,"_meta.data.RData"))
    
    }

###############################

i<-"monocytecount"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
n<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")

a<-data.frame(topgene[[1]],topgene[[2]])
colnames(a)<-c("scPagwastop1000","magmatop1000")
write.csv(a,file=paste0("scPagwas.magmatopgenes_scRDS",i,".csv"),row.names = F)

#############################

#####导出scRDS需要得数据
library(SeuratDisk)
library(Seurat)
DefaultAssay(Single_data) <- "RNA"
SaveH5Seurat(Single_data, "Seu_Hema_addata.h5seurat")
Convert("Seu_Hema_addata.h5seurat", dest="h5ad")
write.csv(Single_data@meta.data,file="Seu_Hema.metadata.csv")
```

#### 计算scRDS，基于上面产生得数据进行计算

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

pbmccells = ['RBCcount','Plateletcount','eosinophilcount','lymphocytecount','monocytecount','neutrophilcount','WhiteBloodCellcount','MeanCorpuscularHemoglobin','MeanCorpusVolume']
pbmccells = ["basophilcount","LymphocytePercent"]

scdrs.preprocess(adata) 
for i in pbmccells:
i="Lymphocytecount3"
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

#### 计算其他富集分析方法：

计算monocyte和basophilcount:

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
#for(j in 1:length(topgene)){
#    genes=paste(topgene[[j]],collapse="\t")
#    aa=paste(names(topgene)[j],names(topgene)[j],genes,sep="\t")
   # write.table(aa,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_",i,".gmt"),row.names=F,col.names=F,quote=F,append=T)
#}
Single_data <- Seurat::AddModuleScore(Single_data, assay = "RNA", topgene, name = c("scPagwastop1000.seruatscore","magmatop1000.seruatscore","scPagwastop500.seruatscore","magmatop500.seruatscore"))

####每次跑都要加上例子数据
counts = load_counts()
se_oj = CreateSeuratObject(counts)
se_oj = cal_PAS(seurat_object = se_oj,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = 'log',
              species = 'mouse', 
              pathway='kegg')

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

#ssgsea_df<-gsva(data.matrix(Single_data@assays$RNA@data), gset.idx.list=topgene, annotation,method="ssgsea")

#zscore_df<-gsva(data.matrix(Single_data@assays$RNA@data), gset.idx.list=topgene, annotation,method="zscore")

#plage_df<-gsva(data.matrix(Single_data@assays$RNA@data), gset.idx.list=topgene, annotation,method="plage")

#df<-data.frame(gsva=unlist(gsva_df[1,]),ssgsea=unlist(ssgsea_df[1,]),zscore=unlist(zscore_df[1,]),plage=unlist(plage_df[1,]),aucell=unlist(auccell_df[,1]))

#save(df,file=paste0(i,"_enrichscore_df.RData"))
#Index<- simulated_cell3@meta.data[rownames(df),"Index_positive"]


```

#### 比较其他方法和scDRS的比较,AUC

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

### 2.对groudtruth得结果进行画图

预处理结果，输入scDRS结果：

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/comparegroudtruth")
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")

pbmccells<-c("RBCcount","Plateletcount","eosinophilcount","basophilcount","LymphocytePercent","lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","MeanCorpuscularHemoglobin","MeanCorpusVolume")
for(i in pbmccells){
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

#### 分别对不同trait进行排秩比较和画图

```R
##################分别计算不同trait得结果
i<-"monocytecount"
i<-"basophilcount"

library('ComplexHeatmap')
library(circlize)

load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_meta.data.RData"))
GroundtruthRankplot(meta.data=meta.data,filecell="monocytecount",celltypes= c("12_CD14.Mono.2","13_CD16.Mono","24_CD8.CM","17_B","25_NK"),
         colors_celltypes=c("#125B50","#E04D01","#F8B400","#FAF5E4","#00AFC1"))
i<-"basophilcount"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_meta.data.RData"))
GroundtruthRankplot(meta.data,filecell="basophilcount",
                    celltypes= c("04_Early.Baso","05_CMP.LMPP","23_CD8.EM","20_CD4.N1"),
         colors_celltypes=c("#EE8572","#ECB390","#E5F4E7","#E5F4E7"))

i<-"MeanCorpuscularHemoglobin"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_meta.data.RData"))
GroundtruthRankplot(meta.data,filecell="MeanCorpuscularHemoglobin",celltypes= c("01_HSC","03_Late.Eryth","02_Early.Eryth","24_CD8.CM","22_CD4.M"),
         colors_celltypes=c("#EE8572","#ECB390","#ECB390","#E5F4E7","#E5F4E7"))

i<-"MeanCorpusVolume"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_meta.data.RData"))
GroundtruthRankplot(meta.data,filecell=i,celltypes= c("01_HSC","03_Late.Eryth","02_Early.Eryth","24_CD8.CM","22_CD4.M"),
         colors_celltypes=c("#EE8572","#ECB390","#ECB390","#E5F4E7","#E5F4E7"))

GroundtruthRankplot(meta.data,filecell="neutrophilcount",celltypes= c("01_HSC","05_CMP.LMPP","21_CD4.N2","22_CD4.M"),
         colors_celltypes=c("#125B50","#F8B400","#FAF5E4","#00AFC1"))

GroundtruthRankplot(meta.data,filecell="eosinophilcount",celltypes= c("01_HSC","16_Pre.B","21_CD4.N2","22_CD4.M"),
         colors_celltypes=c("#125B50","#F8B400","#FAF5E4","#00AFC1"))

GroundtruthRankplot(meta.data,filecell="WhiteBloodCellcount",celltypes= c("01_HSC","19_CD8.N","21_CD4.N2","22_CD4.M"),
         colors_celltypes=c("#125B50","#F8B400","#FAF5E4","#00AFC1"))

GroundtruthRankplot(meta.data,filecell="LymphocytePercent",celltypes= c("01_HSC","20_CD4.N1","11_CD14.Mono.1","13_CD16.Mono"),
         colors_celltypes=c("#125B50","#F8B400","#FAF5E4","#00AFC1"))

GroundtruthRankplot(meta.data,filecell="Hemoglobinconcen",celltypes= c("06_CLP.1","15_CLP.2","21_CD4.N2","22_CD4.M"),
         colors_celltypes=c("#125B50","#F8B400","#FAF5E4","#00AFC1"))

##############函数
GroundtruthRankplot<-function(meta.data,
filecell="monocytecount",compare_type=c("scPagwas_score","scPagwas1000_scdrs.raw_score","magma500_scdrs.raw_score","magma1000_scdrs.raw_score" ,"scPagwas500_scdrs.raw_score"),
celltypes= c("12_CD14.Mono.2","13_CD16.Mono","24_CD8.CM","25_NK"), 
colors_celltypes=c("#EE8572","#ECB390","#E5F4E7","#E5F4E7")){
    
df<-meta.data[meta.data$celltypes %in% celltypes,]

    if("scPagwas_score" %in% compare_type){
        df<-df[order(df$scPagwas_score,decreasing=T),]
pdf(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/comparegroudtruth/Figure_1000scPagwas_score_",filecell,".rankplot.pdf"),
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
      if("scPagwas1000_scdrs.raw_score" %in% compare_type){
df<-df[order(df$scPagwas1000_scdrs.raw_score,decreasing=T),]
pdf(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/comparegroudtruth/Figure_1000scPagwas_scDRS_",filecell,".rankplot.pdf"),
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
####500
     if("scPagwas500_scdrs.raw_score" %in% compare_type){
df<-df[order(df$scPagwas500_scdrs.raw_score,decreasing=T),]
pdf(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/comparegroudtruth/Figure_500scPagwas_scDRS_",filecell,".rankplot.pdf"),
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
    
    if("magma1000_scdrs.raw_score" %in% compare_type){
df<-df[order(df$magma1000_scdrs.raw_score,decreasing=T),]
pdf(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/comparegroudtruth/Figure_1000magma_scDRS_",filecell,".rankplot.pdf"),
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
 if("magma500_scdrs.raw_score" %in% compare_type){    
df<-df[order(df$magma500_scdrs.raw_score,decreasing=T),]
pdf(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/comparegroudtruth/Figure_500magma_scDRS_",filecell,".rankplot.pdf"),
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
    
}

save.image("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/.RData")

###########05_CMP.LMPP

```

![image-20220420105252424](E:\OneDrive\GWAS_Multiomics\scriptmd\image-20220420105252424.png)

#### 计算排秩比例：

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



### 3.groudtruth基因结果表达排序比例图

输出top基因文件：

基于bmmc整体细胞的基因平均表达水平作为排秩的基础，去掉极低表达的基因（无效基因）

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
#library(scPagwas)
#suppressMessages(library(Seurat))
library(org.Hs.eg.db)
library(dplyr)

 pbmccells<-c("RBCcount","Plateletcount","eosinophilcount","basophilcount","LymphocytePercent","lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","MeanCorpuscularHemoglobin","MeanCorpusVolume")
#g2s=toTable(org.Hs.egSYMBOL)
load("Seu_Hema_data_scexpr.RData")
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS")

scPagwas_magma_genelist<-lapply(pbmccells,function(i){
    print(i)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.7.RData"))
Celltype_anno<-Pagwas$Celltype_anno
Celltype_anno$scPagwas_score<-Pagwas$scPagwas_score[rownames(Pagwas$Celltype_anno)]
load(paste0(i,"_magma_genes.RData"))
scPagwas_genes<-names(Pagwas$gene_heritability_correlation[order(Pagwas$gene_heritability_correlation,decreasing=T),])

magmatop1000<-intersect(magma_genes$symbol,rownames(Pagwas$data_mat))
magmatop1000<-magmatop1000[1:1000]
scPagwastop1000<-scPagwas_genes[1:1000]
#magmatop1000<-paste(magmatop1000,collapse=",")
#scPagwastop1000<-paste(scPagwastop1000,collapse=",")
scPagwas_magma_gene<-data.frame(scPagwastop1000,magmatop1000)
colnames(scPagwas_magma_gene)<-c("scPagwastop1000","magmatop1000")
    return(scPagwas_magma_gene)
#save(scPagwas_magma_gene,file=paste0(i,"_scPagwas_magma_gene.RData"))
    })
save(scPagwas_magma_genelist,file="scPagwas_magma_genelist.RData")

##########################
#计算并导出单细胞基因平均表达数据以及var数据
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS")
load("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data_scexpr.RData")
#mean_gene<-apply(Seu_Hema_data_scexpr,1,mean)
i<-"monocytecount"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.7.RData"))
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_meta.data.RData"))

data_mat<-Pagwas$data_mat

#data_mat<-data_mat[apply(data_mat,1,mean)>0.1,]
var_mean_df<-data.frame(
    gene=rownames(data_mat),
 var_x_ge=unlist(apply(data_mat,1,var)),
 mean_x_ge=unlist(apply(data_mat,1,mean)))
 save(var_mean_df,file="bmmc_genevar_mean_df.RData")

#data_mat<-Pagwas$data_mat[names(mean_gene),]
cellvar_mean_df<-data.frame(
    cellid=colnames(data_mat),
 var_x_ge=unlist(apply(data_mat,2,var)),
 mean_x_ge=unlist(apply(data_mat,2,mean))
)
meta.data<-meta.data[colnames(data_mat),]
cellvar_mean_df$CellsrankPvalue<-meta.data$CellsrankPvalue
cellvar_mean_df$positive<-meta.data$positive
cellvar_mean_df$magma_scdrs.pval<-meta.data$magma_scdrs.pval
save(cellvar_mean_df,file="bmmc_cellvar_mean_df.RData")


#!/usr/bin/sh
#PBS -N temp
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
####mild

source activate R4
Rscript /share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/1.r
```



#### top基因富集分析

```R
#######################
#1.所有基因进行简单得基因集合富集分析
#####################
#读入go富集分析结果
library(ggplot2)
library(reshape2)
library(ggpubr)
pbmccells<-c("eosinophilcount","basophilcount","LymphocytePercent","Lymphocytecount3",
             "monocytecount","neutrophilcount","WhiteBloodCellcount","Hemoglobinconcen",
             "MeanCorpuscularHemoglobin","MeanCorpusVolume")

lapply(pbmccells,function(i){
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
#magma_goresult<-read.delim("E:/OneDrive/GWAS_Multiomics/Compare/goanalysis/Lymphocytemagma1000genesGoanalysis/enrichment_results_wg_result1656581733.txt",header = T)
#scPagwas_goresult<-read.delim("E:/OneDrive/GWAS_Multiomics/Compare/goanalysis/LymphocytescPagwas1000genesGoanalysis/enrichment_results_wg_result1656581597.txt",header = T)

###############输出结果
lapply(pbmccells,function(i){
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

###
```

#### 重新进行基因表达排秩和富集分析

```R
###########################6.29
library(scPagwas)
suppressMessages(library(Seurat))
pbmccells<-c("RBCcount","Plateletcount","eosinophilcount","basophilcount","LymphocytePercent","Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/")
 i<-"Lymphocytecount3"
#for(i in pbmccells){
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_",i,".RData"))
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")

#}

load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
geneList<-Pagwas@misc$gene_heritability_correlation[,1]
geneList=sort(geneList,decreasing = T) #从高到低排序
#save(geneList,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"rankgeneList.RData"))

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
#data_mat <- data_mat[apply(data_mat,1,function(x) sum(x!=0)>ncol(Pagwas_rank)*0.01),]
#剩余10824个基因
##################获得
a<-apply(data_mat,1,mean)
a<-sort(a,decreasing=T)
                           
rank_gene<-names(a)
save(rank_gene,file="expr_rank_gene.RData")
g1<-intersect(rank_gene[1:(1/2*length(a))],topgene$scPagwastop1000)
#g2<-intersect(rank_gene[(1/4*length(a)):(2/4*length(a))],topgene$scPagwastop1000)
g2<-intersect(rank_gene[(1/2*length(a)):length(a)],topgene$scPagwastop1000)
#g4<-intersect(rank_gene[(3/4*length(a)):(4/4*length(a))],topgene$scPagwastop1000)

                           

m1<-intersect(rank_gene[1:(1/2*length(a))],topgene$magmatop1000)
m2<-intersect(rank_gene[(1/2*length(a)):length(a)],topgene$magmatop1000)
#m3<-intersect(rank_gene[(2/4*length(a)):(3/4*length(a))],topgene$magmatop1000)
#m4<-intersect(rank_gene[(3/4*length(a)):(4/4*length(a))],topgene$magmatop1000)

#> length(c(m1,m2,m3,m4))
#[1] 691
#> length(c(g1,g2,g3,g4))
#[1] 761
gene_list<-list(g1,g2,m1,m2)
names(gene_list)<-c("g1","g2","m1","m2")
lapply(names(gene_list), function(x){
  write.csv(gene_list[[x]],file=paste0(x,"split2_gene.csv"),row.names = F)
})
#################计算差异分析基因       
#degene<- FindMarkers(Pagwas_rank, ident.1 = "19_CD8.N",logfc.threshold =0,min.pct = 0)
                           
##degene$gene<-rownames(degene)
write.csv(degene,file="de_gene_monocyte.csv")

                          
#write.table(b,file=paste0(x,"_degene.txt"),sep="\t",quote=F,row.names=F)
                           
#lapply(names(gene_list),function(x){    
#    b<- degene[gene_list[[x]],c("gene","avg_log2FC")]
#    write.table(b,file=paste0(x,"_degene.txt"),sep="\t",quote=F,row.names=F)
#})  

#save(gene_list,file=paste0(i,"_genes_split2.RData"))
#}
```

#### 基因排秩比例画图

```R
setwd("E:/OneDrive/GWAS_Multiomics/Compare/5.16compareresult")
#load("bmmc_cellvar_mean_df.RData")
#load("bmmc_genevar_mean_df.RData")
load("D:/OneDrive/GWAS_Multiomics/Compare/goanalysis/expr_rank_gene.RData")
#load("E:/OneDrive/GWAS_Multiomics/Compare/goanalysis/monocytecount_genes_split4.RData")
load("D:/OneDrive/GWAS_Multiomics/Compare/goanalysis/Lymphocytecount3_genes_split2.RData")

percent1<-c(length(c(gene_list$g1))/length(unlist(gene_list[1:2])),
            length(c(gene_list$g2))/length(unlist(gene_list[1:2]))
)
percent2<-c(length(c(gene_list$m1))/length(unlist(gene_list[3:4])),
            #length(gene_list$m2)/length(unlist(gene_list[5:8])),
            length(c(gene_list$m2))/length(unlist(gene_list[3:4]))
            #length(gene_list$m4)/length(unlist(gene_list[5:8]))
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



### 4.对不同的trait的top基因打分在umap图上的映射

```R
suppressMessages(library(Seurat))

setwd("/share/pub/dengcy/GWAS_Multiomics/compare/umap_test")
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")
#suppressMessages(library(Seurat))
Single_data <- FindVariableFeatures(Single_data,nfeatures = 4000)
Single_data <- NormalizeData(Single_data, normalization.method = "LogNormalize", scale.factor = 10000)
Single_data <- ScaleData(Single_data)
Single_data <- RunPCA(object = Single_data, assay = "RNA", npcs = 50)
Single_data <- RunTSNE(object = Single_data,assay =  "RNA", reduction = "pca",dims = 1:50)
Single_data <- RunUMAP(object = Single_data, assay = "RNA", reduction = "pca",dims = 1:50)

fortify.Seurat <- function(x){
  xy <- as.data.frame(Embeddings(x, reduction = "umap"))
  colnames(xy) <- c("UMAP_1", "UMAP_2")
  xy$UMAP_1 <- as.numeric(xy$UMAP_1)
  xy$UMAP_2 <- as.numeric(xy$UMAP_2)

  xy2 <- as.data.frame(Embeddings(x, reduction = "tsne"))
  colnames(xy2) <- c("TSNE_1", "TSNE_2")
  xy2$TSNE_1 <- as.numeric(xy2$TSNE_1)
  xy2$TSNE_2 <- as.numeric(xy2$TSNE_2)
  xy<-cbind(xy,xy2)
  return(cbind(xy, as.data.frame(x@meta.data)))
}

 pbmccells<-c("RBCcount","Plateletcount","eosinophilcount","basophilcount","LymphocytePercent","lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","MeanCorpuscularHemoglobin","MeanCorpusVolume")

for(i in pbmccells){
    print(i)
for()
i<-
Single_meta_data<-fortify.Seurat(Single_data)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_meta.data.RData"))

Immune_anti_plot <- ggplot() +
        geom_point(data = T_fortify_can,
                   aes(x = TSNE_1, y = TSNE_2,color =Anti_inflammatory1), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#FFFCDC",high="#516BEB")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        #ggtitle("G2M score")

        pdf(file="Immune_Anti_inflammatory_plot.pdf",width = 5, height = 5)
        print(Immune_anti_plot)
        dev.off()

```

### BMMC单细胞数据平均得分热图展示

```R
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))

setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test")
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/NM_Healthy_pbmc.rds")

 pbmccells<-c("RBCcount","Plateletcount","eosinophilcount","basophilcount","LymphocytePercent","Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","MeanCorpuscularHemoglobin","MeanCorpusVolume")

 i<-pbmccells[1]
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS")
 for(i in pbmccells){
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



