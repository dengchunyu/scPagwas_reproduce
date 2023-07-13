
## sclinker

```R
library(readr)
SNP<- read_table("/share/pub/dengcy/GWAS_Multiomics/covid19/COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt")
SNP2<-SNP[,c(1:4,7:9,11:13)]
SNP2<-SNP2[,c("CHR","POS","rsid","REF","ALT","all_inv_var_meta_beta","all_inv_var_meta_sebeta","all_inv_var_meta_p","all_meta_AF","all_meta_sample_N")]
colnames(SNP2)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf","N")
colnames(SNP2)<-c("chrom","pos","SNP","A1","A2","BETA","se","pval","maf")
write.table(SNP2,file= "/share/pub/dengcy/GWAS_Multiomics/covid19/COVID19_HGI_gwas_data_ldsc.txt",sep="\t",row.names=F,quote=F)


SNP<- read_table("/share/pub/dengcy/GWAS_Multiomics/ad_test/ieu-b-2.vcf.gz",skip=107)
colnames(SNP)<-c("CHROM","POS","ID","ALT","REF","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
value_list<-lapply(1:nrow(SNP),function(i){
  #index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  #index<-unlist(index)
  value <- strsplit(SNP$IEU[i], split = ":",fixed=T)
  value<-unlist(value)[1:4]
  names(value)<-c("ES","SE","LP","ID")
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","ID")
 value_df<-value_df[,c("ES","SE","LP")]
 gwas_data<-SNP[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$p <-  10^(-gwas_data$LP)
 colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
gwas_data<-gwas_data[,c("CHROM", "POS","ID","ALT","REF","ES","SE","LP")]
colnames(gwas_data)<-c("chrom","pos","SNP","A1","A2","BETA","se","pval")
write.table(gwas_data,file= "/share/pub/dengcy/GWAS_Multiomics/ad_test/AD_gwas_data_ldsc.txt",sep="\t",row.names=F,quote=F)
```

```shell
source activate ldsc
/share/pub/dengcy/software/ldsc-master/munge_sumstats.py \
--sumstats /share/pub/dengcy/GWAS_Multiomics/covid19/COVID19_HGI_gwas_data_ldsc.txt \
--merge-alleles /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/w_hm3.snplist \
--signed-sumstat BETA,0 \
--N  905674 \
--chunksize 500000 \
--out /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI

/share/pub/dengcy/software/ldsc-master/munge_sumstats.py \
--sumstats /share/pub/dengcy/GWAS_Multiomics/ad_test/AD_gwas_data_ldsc.txt \
--merge-alleles /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/w_hm3.snplist \
--signed-sumstat BETA,0 \
--N 63926 \
--chunksize 500000 \
--out /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/AD
```


```R
library(tidyverse)
library("Seurat")

###########################################
library(tidyverse)
library(data.table)
library(dplyr)
gene_coordinates0 <-
  read_tsv("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/NCBI37.3.gene.loc.extendedMHCexcluded",
           col_names = FALSE,col_types = 'cciicc') %>%
  mutate(start=ifelse(X3-50000<0,0,X3-50000),end=X4+50000) %>%
  select(X2,start,end,6,1) %>%
  rename(chr="X2", Gene="X6",EntrezGene="X1") %>%
  mutate(chr=paste0("chr",chr))


#gene_coordinates1 = data.frame(fread("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/ABCpaper_NasserFulcoEngreitz2020_Blood_AvgHiC.txt.gz"))

 df_pre = data.frame(fread("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/AllPredictions.AvgHiC.ABC0.015.minus150.withcolnames.ForABCPaper.txt.gz"))
  df_pre = df_pre[which(df_pre$class == "intergenic" | df_pre$clas == "genic"), ]
  tissuename ="BLD"
  tissuenames2 = as.character(read.table("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/ABC.listbloodQC.txt", header=F)[,1])
 tissue_ids = as.numeric(unlist(sapply(tissuenames2, function(x) return(grep(x, df_pre$CellType)))))
 df = df_pre[tissue_ids, ]
  df2 = cbind.data.frame(df$chr, df$start, df$end, df$TargetGene)
  colnames(df2) = c("chr", "start", "end", "Gene")

  write.table(final_bed, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


Enhancer = read.table("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/Roadmap_Enhancers_Blood.txt",header=F)
df3<-Enhancer[,1:4]
colnames(df3)<-c("chr","start","end","Gene")
df3 = rbind(df3,df2)
gene_coordinates<-df3

gene_coordinates<-merge(gene_coordinates,gene_coordinates0[,c("Gene","EntrezGene")],by="Gene",all.x=TRUE)
gene_coordinates<-gene_coordinates[!is.na(gene_coordinates$EntrezGene),]
save(gene_coordinates,file="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/ABC_gene_coordinates.RData")

###########################################

library(tidyverse)
library(data.table)
library(dplyr)
library(Seurat)
adata1<-readRDS("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds")
Avg_exp_bed(scdata=adata1,celltype="celltype",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/modeldata/")
adata2<-readRDS("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds")
Avg_exp_bed(scdata=adata2,celltype="celltypes",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/realdata/")

adata1<-readRDS("/share/pub/dengcy/GWAS_Multiomics/ad_test/GSE160936.rds")
Avg_exp_bed(scdata=adata1,celltype="annotation",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/AD/")

#covid19
##moderate
load("/share/pub/dengcy/GWAS_Multiomics/covid19/scpagwas_moderate.v1.10.RData")
Avg_exp_bed(scdata=Pagwas,celltype="annotation",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/moderate/")
##Normal
load("/share/pub/dengcy/GWAS_Multiomics/covid19/scpagwas_Normal.v1.10.RData")
Avg_exp_bed(scdata=Pagwas,celltype="annotation",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/Normal/")
##severe
load("/share/pub/dengcy/GWAS_Multiomics/covid19/scpagwas_severe.v1.91.RData")
Avg_exp_bed(scdata=Pagwas,celltype="annotation",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/severe/")

##mild
load("/share/pub/dengcy/GWAS_Multiomics/covid19/scpagwas_mild.v1.10.RData")
Avg_exp_bed(scdata=Pagwas,celltype="annotation",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/mild/")


Avg_exp_bed<-function(scdata,celltype,bedfile,gene_coordinates=gene_coordinates){
average_expression<-AverageExpression(scdata, group.by = celltype)
average_expression<-average_expression$RNA
if(!file.exists(bedfile)){
dir.create(bedfile)
}
ave_func(average_expression,path=bedfile)
}


ave_func<-function(average_expression,path){
average_expression<-na.omit(average_expression)
average_expression<-average_expression[apply(average_expression,1,sum)!=0,]
average_expression<-as.data.frame(average_expression)
top10_function(average_expression,path=path)

}
###top10 genes
top10_function <-function(exp,path){
exp$Gene <- rownames(exp)
exp <- exp %>% add_count(Gene) %>%
  filter(n==1) %>%
  select(-n) %>%
  gather(key = column,value=Expr,-Gene) %>%
  as.tibble()

#############################
###2.Each cell type is scaled to the same total number of molecules.
############################
exp <- exp %>%
  group_by(column) %>%
  mutate(Expr_sum_mean=Expr*1e6/sum(Expr))
##############################
###3.Specificity Calculation
#####################################
exp<- exp %>%
  group_by(Gene) %>%
  mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean)) %>%
  ungroup()
###################
###4.Get MAGMA genes
######################
#Only keep genes that are tested in MAGMA
exp2 <- inner_join(exp,gene_coordinates,by="Gene")
##################
###5.Get number of genes
#################
#Get number of genes that represent 10% of the dataset
n_genes <- length(unique(exp2$EntrezGene))
n_genes_to_keep <- (n_genes * 0.1) %>% round()
##################
###7.Write MAGMA/LDSC input files
#########################
exp2 %>% filter(Expr_sum_mean>1) %>% ldsc_bedfile("column",n_genes_to_keep,path)
print("sucess!")
}

write_group  = function(df,Cell_type,path) {
  df <- select(df,column,chr,start,end,EntrezGene)
  dir.create(paste0(path), showWarnings = FALSE,recursive = TRUE)
  write_tsv(df[-1],paste0(path,make.names(unique(df[1])),".bed"),col_names = F)
  return(df)
}
ldsc_bedfile <- function(d,Cell_type,n_genes_to_keep,path){
  d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity)
  d_spe %>% do(write_group(.,Cell_type,path))
}
```

### 3.5 get annot files

```shell
cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/
source activate mypy
bunzip2 w_hm3.snplist.bz2
tail -n +2 w_hm3.snplist | cut -f 1 > hm_snp.txt

cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/1000G_EUR_Phase3_baseline/
#gunzip baseline.*.annot.gz
for file in baseline.*.annot
do
awk 'NR > 1{print "chr"$1"\t"$2"\t"$2"\t"$3}' $file >> tmp.bed
done
sortBed -i tmp.bed > 1000genomes_phase3_SNPs.bed2
rm tmp.bed

######shell
#!/bin/bash
#SBATCH -e Tcells.err
#SBATCH -o Tcells.out
#SBATCH -J Tcells
#SBATCH -w in010
#SBATCH --mem=50000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate ldsc
#B NK DC monocytes Tcells

test=$1
f=$2
path_name="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/"
all_snps="1000G_EUR_Phase3_baseline/1000genomes_phase3_SNPs.bed2"
all_annotations="1000G_EUR_Phase3_baseline/"
plink_file="1000G_EUR_Phase3_plink/"
hapmap_snps="hm_snp.txt"
weights="1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
frq="1000G_Phase3_frq/1000G.EUR.QC."
echo $test
mkdir "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/"$test"/results/"
mkdir "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/"$test"/results/"
cd "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/"$test
cd "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/"$test

#for f in *.bed
#do
echo $f
intersectBed -c -a $path_name$all_snps -b $f > $f".1000genomes.intersect"
awk '{if($5!=0) print $4}' $f".1000genomes.intersect" > $f".1000genomes.intersect.snp"
mkdir $f"_tissue_dir"
rm $f".1000genomes.intersect"
cd $f"_tissue_dir"
for j in $path_name$all_annotations/*.annot
do
echo $j
file_name=`basename $j`
perl $path_name/fast_match2_minimal.pl ../$f".1000genomes.intersect.snp" $f $j > $file_name
done
gzip *annot
for i in {1..22}
do
/share/pub/dengcy/software/ldsc-master/ldsc.py --l2 --bfile $path_name$plink_file/1000G.EUR.QC.$i --ld-wind-cm 1 --print-snps $path_name$hapmap_snps --annot baseline.$i.annot.gz --out ./baseline.$i
done
cd ..
rm $f".1000genomes.intersect.snp"
```


### 3.6 ldsc

```shell
#!/bin/bash
#SBATCH -e Tcells.err
#SBATCH -o Tcells.out
#SBATCH -J Tcells
#SBATCH -w in010
#SBATCH --mem=50000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate ldsc
sumstats="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/monocytecount_prune.sumstats.gz"
path_name="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/"
weights="1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
frq="1000G_Phase3_frq/1000G.EUR.QC."
all_annotations="1000G_EUR_Phase3_baseline"

cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/modeldata/

for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done
mkdir log_get_pvalues


cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/realdata/

for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done
mkdir log_get_pvalues

sumstats="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/AD.sumstats.gz"
cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/AD/

for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done

sumstats="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI.sumstats.gz"
cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/moderate/

for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done



cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/mild/

for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done

cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/severe/
for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done

cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/Normal/
for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done

```

```R
library("tidyverse")
library("stringr")
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/modeldata/")
files <- list.files(".",pattern=".tissue_dir.results",full.names = TRUE,recursive=T)
#files <- files[grepl("age",files)]

d <- data_frame(filename=files) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames=gsub(".bed_tissue_dir.results","",basename(filename)),
           makenames=gsub(".bed_continuous_tissue_dir.results","",basename(makenames))) %>% unnest()
  d <- d %>%  filter(Category=="L2_1") %>% mutate(P=1-pnorm(`Coefficient_z-score`)) %>%
  mutate(Trait=sub("_.*","",makenames),Cell_Type=gsub("^_","",str_extract(makenames, "_.*"))) %>%
  select(Trait,Cell_Type,Enrichment,Enrichment_std_error,Enrichment_p,P) %>% arrange(Enrichment_p)
    
write_tsv(d,path="modeldata_cell_types_pvalues.txt")


setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/realdata")
files <- list.files(".",pattern=".tissue_dir.results",full.names = TRUE,recursive=T)
d <- data_frame(filename=files) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames=gsub(".bed_tissue_dir.results","",basename(filename)),
           makenames=gsub(".bed_continuous_tissue_dir.results","",basename(makenames))) %>% unnest()
  d <- d %>%  filter(Category=="L2_1") %>% mutate(P=1-pnorm(`Coefficient_z-score`)) %>%
  mutate(Trait=sub("_.*","",makenames),Cell_Type=gsub("^_","",str_extract(makenames, "_.*"))) %>%
  select(Trait,Cell_Type,Enrichment,Enrichment_std_error,Enrichment_p,P) %>% arrange(Enrichment_p)
write_tsv(d,path="realdata_cell_types_pvalues.txt")


setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/AD")
files <- list.files(".",pattern=".tissue_dir.results",full.names = TRUE,recursive=T)
#files <- files[grepl("age",files)]

d <- data_frame(filename=files) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames=gsub(".bed_tissue_dir.results","",basename(filename)),
           makenames=gsub(".bed_continuous_tissue_dir.results","",basename(makenames))) %>% unnest()
  d <- d %>%  filter(Category=="L2_1") %>% mutate(P=1-pnorm(`Coefficient_z-score`)) %>%
  mutate(Trait=sub("_.*","",makenames),Cell_Type=gsub("^_","",str_extract(makenames, "_.*"))) %>%
  select(Trait,Cell_Type,Enrichment,Enrichment_std_error,Enrichment_p,P) %>% arrange(Enrichment_p)
write_tsv(d,path="AD_cell_types_pvalues.txt")

for(i in c("mild","moderate","severe","Normal")){
setwd(paste0("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/",i))
files <- list.files(".",pattern=".tissue_dir.results",full.names = TRUE,recursive=T)
d <- data_frame(filename=files) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames=gsub(".bed_tissue_dir.results","",basename(filename)),
           makenames=gsub(".bed_continuous_tissue_dir.results","",basename(makenames))) %>% unnest()
  d <- d %>%  filter(Category=="L2_1") %>% mutate(P=1-pnorm(`Coefficient_z-score`)) %>%
  mutate(Trait=sub("_.*","",makenames),Cell_Type=gsub("^_","",str_extract(makenames, "_.*"))) %>%
  select(Trait,Cell_Type,Enrichment,Enrichment_std_error,Enrichment_p,P) %>% arrange(Enrichment_p)
  path=paste0(i,"_cell_types_pvalues.txt")
write_tsv(d,path=path)
}
```

```shell
cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/
#for i in Normal severe moderate mild
#cd $i
#ehco $i
for f in *_tissue_dir
do
cd $f
rm -rf monocytecount_prune*
cd ..
done
```




```R
for(i in c("mild","moderate","severe","Normal")){
setwd(paste0("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/",i))
df<-read.table(paste0(i,"_cell_types_pvalues.txt"),header=T)
df$Cell_Type<-gsub("prune_","",df$Cell_Type)
write.table(df,paste0("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/sclinker_",i,"_pvalues.txt"),sep="\t",quote=F,row.names=F)
}

for(i in c("modeldata","realdata","AD")){
setwd(paste0("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/",i))
df<-read.table(paste0(i,"_cell_types_pvalues.txt"),header=T)
df$Cell_Type<-gsub("prune_","",df$Cell_Type)
write.table(df,paste0("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/sclinker_",i,"_pvalues.txt"),sep="\t",quote=F,row.names=F)
}

```

## EPIC

```R
source("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/package/R/EPIC.R")
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC")
gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt"
library(data.table)
gwas_data<-read.table(gwas_data,header=T)
colnames(gwas_data)<-c("chr", "pos", "rsid", "A1", "A2", "beta", "se", "P", "MAF", "N")
gwas_data$Zscore <- gwas_data$beta/gwas_data$se
gwas_data$MAC <- 2* gwas_data$N * gwas_data$MAF
snp_to_gene.Demo <- map_snp_to_gene(gwas = gwas_data,anno_path = "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/package/inst/extdata/annotation_dictionary.txt")
snp_division_obj.Demo <- divide_common_rare(gwas = gwas_data, snp_to_gene = snp_to_gene.Demo)
common_snp_to_gene.Demo <- snp_division_obj.Demo$common_snp_to_gene
rare_snp_to_gene.Demo <- snp_division_obj.Demo$rare_snp_to_gene
library(data.table)
prunein_gene_pval.Demo <- read_in(gwas = gwas_data, snp_to_gene = common_snp_to_gene.Demo)
rare_gene_pval.Demo <- read_in(gwas = gwas_data, snp_to_gene = rare_snp_to_gene.Demo)
save.image(file = "model_data.Demo.RData")
save.image(file = "real_data.Demo.RData")
#Single-cell RNA-seq

library(biomaRt)
library(Seurat)
#model data
Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds"
#Single_data <- readRDS(Single_data)
#scRNA.object <- process_scRNA(SeuratObject = Single_data, meta_ct_col = "celltype")
#realdata
load("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_groundtruth_bmmc_monocyte_v10.RData")
genes<-intersect(rownames(Single_data),rownames(Pagwas_groundtruth))
Pagwas_groundtruth<-Pagwas_groundtruth[genes,]
load("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/grch37.RData")
scRNA.object <- process_scRNA(SeuratObject = Pagwas_groundtruth, meta_ct_col = "celltypes",grch37=grch37)
rm(Single_data)

process_scRNA <- function(SeuratObject, meta_ct_col,grch37) {
  gene.loc = read.table(file = "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/package/inst/extdata/gene.noMHC.loc", sep = "\t",
                        header = FALSE, stringsAsFactors = FALSE)

  HGNC_to_ENSE=getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', "strand"),
                     filters = 'hgnc_symbol',
                     values = rownames(SeuratObject),
                     mart = grch37)

  dictionary = as.data.frame(HGNC_to_ENSE)
  dictionary <- dictionary[dictionary$ensembl_gene_id %in% gene.loc$V1, ]

  dictionary$chromosome_name <- as.numeric(dictionary$chromosome_name)
  o <- order(dictionary$chromosome_name, dictionary$start_position, dictionary$end_position)
  dictionary <- dictionary[o,]

  rawCount <- GetAssayData(object = SeuratObject, assay = "RNA", slot = "counts")
  keep.idx <- which(!is.na(match(rownames(rawCount), dictionary$hgnc_symbol)))
  rawCount <- rawCount[keep.idx, ]
  rawCount <- rawCount[match(dictionary$hgnc_symbol, rownames(rawCount)), ]

  # RPKM
  gene.length <- (dictionary$end_position - dictionary$start_position)/1000
  library.size <- apply(rawCount, 2, sum)
  rpkm <- (10^6)*t(t(as.matrix(rawCount))/library.size)/gene.length
  cts <- sort(unique(SeuratObject@meta.data[[meta_ct_col]]))

  log2_rpkm_ct <- matrix(NA, nrow = nrow(rpkm), ncol = length(cts))
  rownames(log2_rpkm_ct) <- rownames(rpkm)
  colnames(log2_rpkm_ct) <- cts
  for (ct in cts) {
    cell.idx = rownames(SeuratObject@meta.data)[which(SeuratObject@meta.data[[meta_ct_col]] == ct)]
    log2_rpkm_ct[,ct] <- log2(1 + apply(rpkm[, cell.idx, drop = FALSE], 1, mean))
  }

  SeuratObject <- subset(SeuratObject, features = rownames(rpkm))

  ngene = 8000
  var.features <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = ngene)

  common.genes.ngene <- var.features@assays$RNA@var.features
  dictionary.ngene <- dictionary[dictionary$hgnc_symbol %in% common.genes.ngene, ]
  log2_rpkm_ct.ngene <- log2_rpkm_ct[rownames(log2_rpkm_ct) %in% common.genes.ngene, ]
  colnames(log2_rpkm_ct.ngene) <- gsub(" ", "_", colnames(log2_rpkm_ct.ngene))

  log2_rpkm_ct.ngene = cbind(log2_rpkm_ct.ngene, apply(log2_rpkm_ct.ngene, 1, mean))
  colnames(log2_rpkm_ct.ngene)[ncol(log2_rpkm_ct.ngene)] = "Average"
  rownames(log2_rpkm_ct.ngene) <- dictionary$ensembl_gene_id[match(rownames(log2_rpkm_ct.ngene), dictionary$hgnc_symbol)]

  log2_rpkm_ct.ngene = as.data.frame(log2_rpkm_ct.ngene)

  dictionary.ngene <- dictionary.ngene[,-1]
  colnames(dictionary.ngene) <- c("V1","V2","V3","V4","V5")

  return(list(scRNA.rpkm = log2_rpkm_ct.ngene, scRNA.loc = dictionary.ngene))
}


scRNA.rpkm <- scRNA.object$scRNA.rpkm
scRNA.loc <- scRNA.object$scRNA.loc
head(scRNA.rpkm)
library(genio)
file_bed <- file.path("/share/pub/dengcy/1000G_EUR_Phase3_plink", "all_chr_1000G.bed")
file_bim <- file.path("/share/pub/dengcy/1000G_EUR_Phase3_plink", "all_chr_1000G.bim")
file_fam <- file.path("/share/pub/dengcy/1000G_EUR_Phase3_plink", "all_chr_1000G.fam")
bim <- read_bim(file_bim)
fam <- read_fam(file_fam)
genotype <- read_bed(file_bed, bim$id, fam$id)

X_ref <- genotype[match(bim$id, rownames(genotype)),]
rm(bim,fam )
save.image(file="real_data.Demo.RData")
save.image(file="model_data.Demo.RData")

inter_dir <-"/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/modeldata/inter"
inter_dir <-"/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/realdata/inter"
save(prunein_gene_pval.Demo,X_ref,scRNA.loc,file="real_poet_sw.RData")
for(chr in 1:22){

calculate_POET_sw(genotype = X_ref, gene_pval = prunein_gene_pval.Demo, 
                  gene.loc = scRNA.loc, # For scRNA-seq data, specify gene.loc = scRNA.loc, which is returned by function process_scRNA()
                  chr =chr, type = "POET", inter_dir = inter_dir)
}
load("real_data.Demo.RData")
pruneinObject.Demo <- get_gene_chisq(gene_pval = prunein_gene_pval.Demo, 
                                     gene.loc = scRNA.loc, # For scRNA-seq data, specify gene.loc = scRNA.loc, which is returned by function process_scRNA() 
                                     type = "POET", inter_dir = inter_dir)
save.image(file = "model_data.Demo.RData")
save.image(file = "real_data.Demo.RData")

prunein_POET_gene_pval.Demo <- pruneinObject.Demo$gene_POET_pval
prunein_gene_corr.mat.Demo <- pruneinObject.Demo$gene_corr.mat
#4.2 Burden test for rare variants
rare_gene_pval.Demo <- get_burden(genotype = X_ref, gene_pval = rare_gene_pval.Demo)

#4.3 Joint analysis for common and rare variants
combine_gene_pval.Demo <- get_combined(prunein_gene_pval = prunein_gene_pval.Demo, rare_gene_pval = rare_gene_pval.Demo)
X_super.Demo <- construct_X_super(genotype = X_ref, rare_gene_pval = rare_gene_pval.Demo)
save.image(file = "model_data.Demo.RData")
save(X_ref,combine_gene_pval.Demo,scRNA.loc,X_super.Demo,file = "model_data.poet2.RData")

save(X_ref,combine_gene_pval.Demo,scRNA.loc,X_super.Demo,file = "real_data.poet2.RData")

load( "model_data.poet2.RData")
for(chr in 1:22){
calculate_POET_sw(genotype = X_ref, gene_pval = combine_gene_pval.Demo, 
                  gene.loc = scRNA.loc, # For scRNA-seq data, specify gene.loc = scRNA.loc, which is returned by function process_scRNA()
                  chr =chr, type = "iPOET", X_super = X_super.Demo, inter_dir = inter_dir)
}
save(X_ref,combine_gene_pval.Demo,scRNA.loc,X_super.Demo,file = "model_data.poet2.RData")
save.image(file = "model_data.Demo.RData")
save.image(file = "real_data.Demo.RData")
combineObject.Demo <- get_gene_chisq(gene_pval = combine_gene_pval.Demo, 
                                     gene.loc = scRNA.loc, # For scRNA-seq data, specify gene.loc = scRNA.loc, which is returned by function process_scRNA() 
                                     type = "iPOET", inter_dir = inter_dir)
save.image(file = "model_data.Demo.RData")
combine_POET_gene_pval.Demo <- combineObject.Demo$gene_POET_pval
combine_gene_corr.mat.Demo <- combineObject.Demo$gene_corr.mat
#5. Prioritizing trait-relevant tissue(s) and cell type(s)
gtex.enrichment.Demo <- prioritize_relevance(gene_pval = prunein_POET_gene_pval.Demo, gene_corr.mat = prunein_gene_corr.mat.Demo, 
                                             gene_expr = scRNA.rpkm, gene.loc = scRNA.loc, chrs =1:22) 

gtex.enrichment.Demo <- prioritize_relevance(gene_pval = prunein_POET_gene_pval.Demo, gene_corr.mat = prunein_gene_corr.mat.Demo, 
                                             gene_expr = scRNA.rpkm, gene.loc = scRNA.loc, chrs =1:22) 

#gtex.enrichment.joint.Demo <- prioritize_relevance(gene_pval = #combine_POET_gene_pval.Demo, gene_corr.mat = combine_gene_corr.mat.Demo, 
                                                   #gene_expr = scRNA.rpkm, gene.loc = scRNA.loc, chrs = 1:22)
save(gtex.enrichment.joint.Demo,file = "epic_modeldata_result.RData")
save(gtex.enrichment.Demo,file = "epic_realdata_result.RData")
save.image(file = "model_data.Demo.RData")
save.image(file = "real_data.Demo.RData")

df<-as.data.frame(gtex.enrichment.joint.Demo)
colnames(df)<-"p"
df$adjp<-p.adjust(df$p,method = "BH")
write.csv(df,file = "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/EPIC_modeldata_cell_types_pvalues.txt")

df<-as.data.frame(gtex.enrichment.Demo)
colnames(df)<-"p"
df$adjp<-p.adjust(df$p,method = "BH")
write.csv(df,file = "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/EPIC_realdata_cell_types_pvalues.txt")
```

```shell
#!/bin/bash
#SBATCH -e epic2.err
#SBATCH -o epic2.out
#SBATCH -J epic2
#SBATCH -w in006
#SBATCH --mem=30000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate R4.2
Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/2.r
Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/1.r
```