# GWAS summary data progressed

## Deal with blood traits gwas

### basophil cell count

Dataset: ieu-b-29

```R
library(readr)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits")
############
basophilcount<-"ieu-b-29.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(basophilcount,skip=109)

colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
    
  value<-value[c("ES","SE","LP","AF","SS")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=5)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF","SS")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf","N")
# if(is.na(gwas_data$maf[1])) 
    # gwas_data$maf<-0.1
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/basophilcount_gwas_data.txt",row.names=F,quote=F)
```

### eosinophil cell count

Dataset: ieu-b-33

```R
library(readr)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits")
############
eosinophilcount<-"ieu-b-33.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(eosinophilcount,skip=109)

colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
    
  value<-value[c("ES","SE","LP","AF","SS")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=5)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF","SS")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf","N")
# if(is.na(gwas_data$maf[1])) 
    # gwas_data$maf<-0.1
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/eosinophilcount_gwas_data.txt",row.names=F,quote=F)
```

### monocyte cell count

Dataset: ieu-b-31

```R
library(readr)
setwd("/home/guofj/bloodtraits")
############
monocytecount<-"ieu-b-31.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table2(monocytecount,,comment = "#")

colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
    
  value<-value[c("ES","SE","LP","AF","SS")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=5)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF","SS")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf","N")
# if(is.na(gwas_data$maf[1])) 
    # gwas_data$maf<-0.1
write.table(gwas_data,file="monocytecount_gwas_data.txt",row.names=F,quote=F)
```

### neutrophil cell count

Dataset: ieu-b-34

```R
library(readr)
setwd("/home/guofj/bloodtraits")
############
neutrophilcount<-"ieu-b-34.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table2(neutrophilcount,comment = "#")

colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
    
  value<-value[c("ES","SE","LP","AF","SS")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=5)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF","SS")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf","N")
# if(is.na(gwas_data$maf[1])) 
    # gwas_data$maf<-0.1
write.table(gwas_data,file="neutrophilcount_gwas_data.txt",row.names=F,quote=F)
```

### white blood cell count

Dataset: ieu-b-30

```R
library(readr)
setwd("/home/guofj/bloodtraits")
############
WhiteBloodCellcount<-"ieu-b-30.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table2(WhiteBloodCellcount,comment = "#")

colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
    
  value<-value[c("ES","SE","LP","AF","SS")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=5)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF","SS")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf","N")
# if(is.na(gwas_data$maf[1])) 
    # gwas_data$maf<-0.1
write.table(gwas_data,file="WhiteBloodCellcount_gwas_data.txt",row.names=F,quote=F)
```

### 

### Hemoglobin concentration

Dataset: ebi-a-GCST90002311

```R
library(readr)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits")
############
Hemoglobinconcen<-"ebi-a-GCST90002311.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(Hemoglobinconcen,comment = "#")

colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
    
  value<-value[c("ES","SE","LP","AF")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf")
# if(is.na(gwas_data$maf[1])) 
    # gwas_data$maf<-0.1
 gwas_data$N<-1000
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/Hemoglobinconcen_gwas_data.txt",row.names=F,quote=F)
```



### Mean corpuscular hemoglobin concentration

Dataset: ebi-a-GCST90002329

```R
library(readr)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits")
############
MeanCorpuscularHemoglobin<-"ebi-a-GCST90002329.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(MeanCorpuscularHemoglobin,comment = "#")

colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
    
  value<-value[c("ES","SE","LP","AF")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf")
# if(is.na(gwas_data$maf[1])) 
     gwas_data$N<-1000
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/MeanCorpuscularHemoglobin_gwas_data.txt",row.names=F,quote=F)
```



### Mean corpuscular volume

Dataset: ebi-a-GCST90002335

```R
library(readr)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits")
############
MeanCorpusVolume<-"ebi-a-GCST90002335.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(MeanCorpusVolume,comment = "#")

colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
    
  value<-value[c("ES","SE","LP","AF")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf")
# if(is.na(gwas_data$maf[1])) 
    # gwas_data$maf<-0.1
 gwas_data$N<-1000
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/MeanCorpusVolume_gwas_data.txt",row.names=F,quote=F)
```

### **Lymphocyte count**

**Lymphocyte count**3
**Dataset: bbj-a-36**

```R
library(readr)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits")
############
Lymphocyte_count<-"bbj-a-36.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(Lymphocyte_count,comment = "#")

colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  #index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  #index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)[1:4]
  names(value)<-c("ES","SE","LP","AF")
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF")
 gwas_data<-GWAS_raw[,c(1:4)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$p <-  10^(-gwas_data$LP)
 colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
gwas_data<-gwas_data[,c("CHROM", "POS","ID","REF","ALT","ES","SE","p","AF")]
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf")

write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/Lymphocytecount3_gwas_data.txt",row.names=F,quote=F)
```

## blood traits snp pruning

This process is to reduce the number of snp to save time and memory

```shell
cd /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits
mkdir tempfile
#######################1.
#neutrophilcount WhiteBloodCellcount MeanCorpusVolume 
#Plateletcount basophilcount eosinophilcount lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount MeanCorpusVolume 
###
for i in  
do
awk  '{print $3 }' ${i}_gwas_data.txt  > /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/${i}_SNP_list.txt
done


######################2.
#for j in RBCcount Plateletcount basophilcount eosinophilcount lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount MeanCorpusVolume 

for j in RBCcount2 Lymphocytepercentage2 
do
echo $j
for i in $(seq 1 22)  
do 
echo $i
plink --bfile /share/pub/mayl/ABC_model/ldsc/data_ldsc/02_partitioned_LD_score_estimation/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i --extract /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/${j}_SNP_list.txt --noweb --make-bed --out /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/1000G.EUR.QC.${j}_${i}_filtered
done 
done

######################3.
#for j in eosinophilcount lymphocytecount monocytecount 
for j in Lymphocytepercentage2
do
echo $j
for i in $(seq 1 22)  
do 
echo $i
plink --bfile /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

######################4.
cd /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile
#for j in RBCcount Plateletcount basophilcount eosinophilcount lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount LymphocytePercent Hemoglobinconcen MeanCorpuscularHemoglobin MeanCorpusVolume
for j in RBCcount2 Lymphocytepercentage2
do
echo $j
cat [${j}]*.prune.in > ${j}_merge_plink_EUR_filtered_LD0.8.prune.in
done 
#####remove 
rm -f {*.prune.out}
##
for j in RBCcount Plateletcount basophilcount eosinophilcount lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount LymphocytePercent Hemoglobinconcen MeanCorpuscularHemoglobin MeanCorpusVolume
do
echo $j
wc -l ${j}_merge_plink_EUR_filtered_LD0.8.prune.in
wc -l /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/${j}_gwas_data.txt
done 
```

## Integrate blood traits snp files

```R
library(readr)
library(dplyr)
for (i in c("basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")) {
    print(i)
    gwas<-read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_gwas_data.txt"))
    SNP_prune<- read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/",i,"_merge_plink_EUR_filtered_LD0.8.prune.in"))
    SNP_prune<-SNP_prune[!duplicated(unlist(SNP_prune)),]
    colnames(SNP_prune)<-"rsid"
    #### Left Join using inner_join function 
   gwas= gwas %>% inner_join(SNP_prune,by="rsid")
    print(nrow(gwas))
   write.table(gwas,file= paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"),,row.names=F,quote=F)
        print(i)
#gc()
}
```

