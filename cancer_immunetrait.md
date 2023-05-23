## cancer immunetrait

```R
library("scPagwas")
library("Seurat")
library("SingleCellExperiment")
library("stringr") 

setwd("/share/pub/dengcy/GWAS_Multiomics/gwasdata")
 load("GWAS_NK.activated_proportion_median.low.high.RData")
 result.maf.ls<-result.maf.ls[result.maf.ls$MAF>0.01, ]

 result.maf.ls<-result.maf.ls[,c("CHR","SNP","A1","A2","BP","MAF","STAT","P")]
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","maf")
 result.maf.ls
```

singlecell

```R
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/CRC")
crc_sc<-readRDS("CRC_merged_all.rds")
```

