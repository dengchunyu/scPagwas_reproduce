## TWAS analysis
### create sumstats file
```r
rm(list=ls())
library(data.table)
data<-fread("LymphocytePercent_gwas_data.txt",header=T,sep=" ")
data$Z<-data$beta/data$se
data<-data[,c(3:5,11)]
data1<-data[,c(1,3,2,4)]
colnames(data1)<-c("SNP","A1","A2","Z")
write.table(data1,file="LymphocytePercent_gwas_data.sumstats",sep=" ",quote=F)
``````
### run TWAS

```shell
#PBS -N fusion
#PBS -q workq
#PBS -l mem=150gb
#PBS -l ncpus=2
cd /share2/pub/chenchg/chenchg/TWAS/
source activate R4.2.2
for chr in $(seq 1 22);
do
Rscript /share2/pub/chenchg/chenchg/software/TWAS_fusion/fusion_twas-master/FUSION.assoc_test.R \
--sumstats /share2/pub/chenchg/chenchg/TWAS/WhiteBloodCellcount_gwas_data.sumstats \
--weights  /share2/pub/chenchg/chenchg/TWAS/WEIGHTS_Whole_Blood/GTExv8.ALL.Whole_Blood.pos \
--weights_dir /share2/pub/chenchg/chenchg/TWAS/WEIGHTS_Whole_Blood/ \
--ref_ld_chr /share2/pub/chenchg/chenchg/TWAS/LDREF/1000G.EUR. \
--chr $chr \
--out /share2/pub/chenchg/chenchg/TWAS/result/WhiteBloodCellcount_result/WhiteBloodCellcount_gwas_data.$chr.dat;
done

```