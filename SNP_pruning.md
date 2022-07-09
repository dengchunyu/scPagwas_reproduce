# GWAS数据处理原始数据

处理成可应用与软件的数据（scPagwas,rolypoly,MAGMA,LDSC）

统计有42个GWAS数据，包含

- Neuropsychiatric disorders：MDD，ADHD，BIB，ASD；

- Neurodegenerative disorders：AD，PD，EP，Focal_EP，Generalized_EP，Juvenile_EP，MG，MS；

- Behavior-cognitive phenotype：SC，Smoking_status，CognitivePerformance，Openness，Conscientiousness，Subjective well-being，Happiness；

- Metabolic disease：BMI，WHR，WC，OB，CKD，TC，LDL ，HDL；

- Autoimmune disease：IBD，PBC，CoeliacDisease，Crohn'sDisease，T1D，ATD；

- Cardiovascular disease：CAD，DBP_European，DBP_EastAsian，PR_EastAsian，SBP，T2D；

等相关表型特征。

## Read the raw data for gwas summary statics

### 原始数据地址

ls查看文件夹下文件

Neuropsychiatric disorders：/share/pub/mayl/Sherlock/00_IEU_GWAS/01_neuropsychiatric_disorders

```
ADHD  AN  ASD  BIB  GWAS_SCZ_2020  GWAS_UK_Biobank_MDD
```

Neurodegenerative disorders：/share/pub/mayl/Sherlock/00_IEU_GWAS/02_neurodegenerative_disorders

```
AD  ALS  EP  Focal_EP  Generalized_EP  HD  insomnia  IntracranialVolume  IS  Juvenile_EP  Migraine  Multiple_sclerosis
```

Behavior-cognitive phenotype：/share/pub/mayl/Sherlock/00_IEU_GWAS/03_Behavior_cognitive_phenotype

```
AgeOfSmokingInitiation  Alcohol_per_weeks     Conscientiousness  DaytimeDozing   DepressiveSymptoms  EduYears      GettingUp  Intelligence  Neuroticism          SleepDuration     Snoring           TAG_smoking_GWAS
Agreeableness           CognitivePerformance  CPD                DaytimeNapping  EduAttainment       Extraversion  Happiness  Morningness   OpennessToExperence  smokingCessation  SubjectWellBeing
```

Metabolic disease：/share/pub/mayl/Sherlock/00_IEU_GWAS/04_Metabolic_disease

```
BMI  BodyFatPercentage  CKD  HDL  Height  HTY  LDL  Obesity  TC  TG  WC  WHR
```

Autoimmune disease：/share/pub/mayl/Sherlock/00_IEU_GWAS/05_Immunological_disease

```
#ATD  CoeliacDisease  #CrohnDisease  IBD  PBC  Psoriasis  RA  SLE  #T1D  T1D_Nature_2021  UlcerativeColitis
```

Cardiovascular disease：/share/pub/mayl/Sherlock/00_IEU_GWAS/06_Cardiovascular_disease

```
CAD  DBP_EastAsian  DBP_European  DBP_Hispanic  HBP  Hypertension  PR_EastAsian  SBP  T2D
```

### 分别对所有数据进行预处理

进入R对原始数据进行数据选择操作，并输出txt文件。

#### **1.Neuropsychiatric disorders**

/share/pub/mayl/Sherlock/00_IEU_GWAS/01_neuropsychiatric_disorders/GWAS_UK_Biobank_MDD/UK_biobank_data_for_depression/ukb-b-12064.vcf.gz

/share/pub/mayl/Sherlock/00_IEU_GWAS/01_neuropsychiatric_disorders/ADHD/ieu-a-1183.vcf.gz

/share/pub/mayl/Sherlock/00_IEU_GWAS/01_neuropsychiatric_disorders/ASD/ieu-a-1185.vcf.gz

/share/pub/mayl/Sherlock/00_IEU_GWAS/01_neuropsychiatric_disorders/BIB/ieu-b-41.vcf.gz

```R
library(readr)
setwd("/share/pub/mayl/Sherlock/00_IEU_GWAS/01_neuropsychiatric_disorders")
############
a<-c(ADHD="./ADHD/ieu-a-1183.vcf.gz",
ASD="./ASD/ieu-a-1185.vcf.gz"
MDD="./GWAS_UK_Biobank_MDD/UK_biobank_data_for_depression/ukb-b-12064.vcf.gz",
BIB="./BIB/ieu-b-41.vcf.gz")

##开始循环
lapply(names(a),function(x){
GWAS_raw <-read_table(a[x],skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
  #ES:SE:LP:AF:ID
  value<-value[c("ES","SE","LP","ID")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","ID")
 gwas_data<-GWAS_raw[,c(1:3,6)]
 gwas_data<- cbind(gwas_data,value_df)
  gwas_data$LP <-  10^(-gwas_data$LP)
 colnames(gwas_data)<-c("chrom","pos","rsid","N","beta","se","p","maf")
 if(is.na(gwas_data$maf[1])) gwas_data$maf<-0.1
 
write.table(gwas_data,file=paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/01_neuropsychiatric_disorders/",x,"_gwas_data.txt"),row.names=F,quote=F)
})

######################
b<-c(PD="/share/pub/dengcy/GWAS_Multiomics/test/brain/gwas_files/ieu-b-7.vcf.gz")

GWAS_raw <-read_table(b[1],skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","maf","FORMAT","IEU")
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
  #ES:SE:LP:AF:ID
  value<-value[c("ES","SE","LP","ID")]
  return(value)
  })

 value_df<- matrix(unlist(value_list),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","ID")
 
 if(GWAS_raw$maf[1] != "."){
 GWAS_raw$maf<-unlist(lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$maf[i], split = "=",fixed=T)
  index<-unlist(index)[2]
  return(index)
 }))
 }
 
 gwas_data<-GWAS_raw[,c(1:3,6,8)]
 gwas_data<- cbind(gwas_data,value_df[,1:3])
  gwas_data$LP <-  10^(-as.numeric(gwas_data$LP))
 colnames(gwas_data)<-c("chrom","pos","rsid","N","maf","beta","se","p")
# if(sum(is.na(gwas_data$maf))>2/3*nrow(gwas_data)) gwas_data$maf<-0.1
# gwas_data$maf<-0.1
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/PD_gwas_data.txt",row.names=F,quote=F)

rm(list=ls())
gc()
```

#### 2.Neurodegenerative disorders

/share/pub/mayl/Sherlock/00_IEU_GWAS/02_neurodegenerative_disorders/AD/ieu-b-2.vcf.gz

/share/pub/mayl/Sherlock/00_IEU_GWAS/02_neurodegenerative_disorders/EP/ieu-b-8.vcf.gz

/share/pub/mayl/Sherlock/00_IEU_GWAS/02_neurodegenerative_disorders/Focal_EP/ieu-b-10.vcf.gz

/share/pub/mayl/Sherlock/00_IEU_GWAS/02_neurodegenerative_disorders/Generalized_EP/ieu-b-9.vcf.gz

/share/pub/mayl/Sherlock/00_IEU_GWAS/02_neurodegenerative_disorders/Juvenile_EP/ieu-b-17.vcf.gz

/share/pub/mayl/Sherlock/00_IEU_GWAS/02_neurodegenerative_disorders/Migraine/ukb-b-16868.vcf.gz

/share/pub/mayl/Sherlock/00_IEU_GWAS/02_neurodegenerative_disorders/Multiple_sclerosis/ieu-b-18.vcf.gz

```R
library(readr)
setwd("/share/pub/mayl/Sherlock/00_IEU_GWAS/02_neurodegenerative_disorders")
############
a<-c(AD="AD/ieu-b-2.vcf.gz",
EP="EP/ieu-b-8.vcf.gz",
Focal_EP="Focal_EP/ieu-b-10.vcf.gz",
Generalized_EP="Generalized_EP/ieu-b-9.vcf.gz",
Juvenile_EP="Juvenile_EP/ieu-b-17.vcf.gz",
Migraine="Migraine/ukb-b-16868.vcf.gz",
Multiple_sclerosis="Multiple_sclerosis/ieu-b-18.vcf.gz")
##根据不同数据的格式进行整理

lapply(names(a),function(x){
GWAS_raw <-read_table(a[x],skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
  #ES:SE:LP:AF:ID
  value<-value[c("ES","SE","LP","AF")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF")
 gwas_data<-GWAS_raw[,c(1:3,6)]
 gwas_data<- cbind(gwas_data,value_df)
  gwas_data$LP <-  10^(-gwas_data$LP)
 colnames(gwas_data)<-c("chrom","pos","rsid","N","beta","se","p","maf")
 if(is.na(gwas_data$maf[1])) gwas_data$maf<-0.1
 
write.table(gwas_data,file=paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",x,"_gwas_data.txt"),row.names=F,quote=F)
})
rm(list=ls())
gc()
```

#### 3.Behavior-cognitive phenotype

```
library(readr)
setwd("/share/pub/mayl/Sherlock/00_IEU_GWAS/03_Behavior_cognitive_phenotype")
############
a<-c(smokingCessation="./smokingCessation/bbj-a-81.vcf.gz.vcf.gz",
OpennessToExperence="./OpennessToExperence/ieu-a-117.vcf.gz.vcf.gz",
Happiness="./Happiness/ukb-b-4062.vcf.gz.vcf.gz")

a<-c(CognitivePerformance="./CognitivePerformance/ebi-a-GCST006572.vcf.gz",
SubjectWellBeing="./SubjectWellBeing/ebi-a-GCST003766.vcf.gz")

##开始循环
lapply(names(a),function(x){
GWAS_raw <-read_table(a[x],skip=150)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
  #ES:SE:LP:AF:ID
  value<-value[c("ES","SE","LP","AF")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF")
 gwas_data<-GWAS_raw[,c(1:3,6)]
 gwas_data<- cbind(gwas_data,value_df)
  gwas_data$LP <-  10^(-gwas_data$LP)
 colnames(gwas_data)<-c("chrom","pos","rsid","N","beta","se","p","maf")
 if(is.na(gwas_data$maf[1])) gwas_data$maf<-0.1
 
write.table(gwas_data,file=paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",x,"_gwas_data.txt"),row.names=F,quote=F)
})
rm(list=ls())
gc()
############读取另一种类型的数据



```

#### 4.Metabolic_disease

/share/pub/mayl/Sherlock/00_IEU_GWAS/04_Metabolic_disease

```
library(readr)
setwd("/share/pub/mayl/Sherlock/00_IEU_GWAS/04_Metabolic_disease")
############
a<-c(BMI="./BMI/ukb-b-19953.vcf.gz",
WHR="./WHR/ieu-a-72.vcf.gz",
WC="./WC/ieu-a-60.vcf.gz",
Obesity="./Obesity/ieu-a-90.vcf.gz",
CKD="./CKD/ebi-a-GCST003374.vcf.gz",
TC="./TC/TC_2013/ieu-a-301.vcf.gz",
LDL="./LDL/LDL_2013/ieu-a-300.vcf.gz",
HDL="./HDL/HDL_2013/ieu-a-299.vcf.gz")

##开始循环
lapply(names(a),function(x){
GWAS_raw <-read_table(a[x],skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
  #ES:SE:LP:AF:ID
  value<-value[c("ES","SE","LP","AF")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF")
 gwas_data<-GWAS_raw[,c(1:3,6)]
 gwas_data<- cbind(gwas_data,value_df)
  gwas_data$LP <-  10^(-gwas_data$LP)
 colnames(gwas_data)<-c("chrom","pos","rsid","N","beta","se","p","maf")
 if(is.na(gwas_data$maf[1])) gwas_data$maf<-0.1
 
write.table(gwas_data,file=paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",x,"_gwas_data.txt"),row.names=F,quote=F)
})
rm(list=ls())
gc()
```

#### 5. Immunological_disease

```
library(readr)
setwd("/share/pub/mayl/Sherlock/00_IEU_GWAS/05_Immunological_disease")
############

a<-c(IBD="./IBD/ebi-a-GCST004131.vcf.gz",
PBC="./PBC/ebi-a-GCST003129.vcf",
UlcerativeColitis="./UlcerativeColitis/ebi-a-GCST004133.vcf.gz")

#lapply(names(a),function(x){
GWAS_raw <-read_table(a[3],skip=120)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","maf","FORMAT","IEU")
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
  #ES:SE:LP:AF:ID
  value<-value[c("ES","SE","LP","ID")]
  return(value)
  })

 value_df<- matrix(unlist(value_list),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","ID")
 if(GWAS_raw$maf[1] != "."){
 GWAS_raw$maf<-unlist(lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$maf[i], split = "=",fixed=T)
  index<-unlist(index)[2]
  return(index)
 }))}
 
 gwas_data<-GWAS_raw[,c(1:3,6,8)]
 gwas_data<- cbind(gwas_data,value_df[,1:3])
  gwas_data$LP <-  10^(-as.numeric(gwas_data$LP))
 colnames(gwas_data)<-c("chrom","pos","rsid","N","maf","beta","se","p")
# if(sum(is.na(gwas_data$maf))>2/3*nrow(gwas_data)) gwas_data$maf<-0.1
 gwas_data$maf<-0.1
write.table(gwas_data,file=paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",x,"_gwas_data.txt"),row.names=F,quote=F)
})
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/UlcerativeColitis_gwas_data.txt",row.names=F,quote=F)
###################################另一种类型的数据读取

a<-c(
CoeliacDisease="./CoeliacDisease/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz",
#CrohnDisease="./CrohnDisease/franke_2010_21102463_cro_efo0000384_1_gwas.sumstats.tsv.gz",#数据不合适
#ATD="./ATD/cooper_2012_22922229_atd_efo0006812_1_ichip.sumstats.tsv.gz"#数据不合适
T1D="./T1D_Nature_2021/GCST90014023_buildGRCh38.tsv"
)

####CoeliacDisease
GWAS_raw <-read_table(a[1])
###这个文件的格式和之前的完全不一样，更容易读取
#cols(
#  chrom = col_double(),
#  pos = col_double(),
#  rsid = col_character(),
#  other_allele = col_character(),
#  effect_allele = col_character(),
#  p = col_double(),
#  beta = col_double(),
#  se = col_double(),
#  OR = col_double(),
#  OR_lower = col_double(),
#  OR_upper = col_double()
#)

gwas_data<-GWAS_raw[,c(1,2,3,6,7,8)]
colnames(gwas_data)<-c("chrom","pos","rsid","p","beta","se")
  gwas_data$maf<-0.1
  gwas_data$N<-4533
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/CoeliacDisease_gwas_data.txt",row.names=F,quote=F)

####CrohnDisease
GWAS_raw <-read_table(a[2])
gwas_data<-GWAS_raw[,c(1,2,3,6,7,8)]
colnames(gwas_data)<-c("chrom","pos","rsid","p","beta","se")
  gwas_data$maf<-0.1
  gwas_data$N<-22027
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/CrohnDisease_gwas_data.txt",row.names=F,quote=F)

####T1D
GWAS_raw <-read_table(a[3])
#cols(
#  variant_id = col_character(),
#  p_value = col_double(),
#  chromosome = col_double(),
#  base_pair_location = col_double(),
#  effect_allele = col_character(),
#  other_allele = col_character(),
#  effect_allele_frequency = col_double(),
#  beta = col_double(),
#  standard_error = col_double(),
#  sample_size = col_double()
#)

gwas_data<-GWAS_raw[,c(1,2,3,4,7,8,9,10)]
colnames(gwas_data)<-c("rsid","p","chrom","pos","maf","beta","se","N")
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/gwasdata/T1D_gwas_data.txt",row.names=F,quote=F)

####ATD,beta为NA没法用的数据
GWAS_raw <-read_table(a[4])
gwas_data<-GWAS_raw[,c(1,2,3,6,7,8)]
colnames(gwas_data)<-c("chrom","pos","rsid","p","beta","se")
  gwas_data$maf<-0.1
  gwas_data$N<-"."
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/ATD_gwas_data.txt",row.names=F,quote=F)


rm(list=ls())
gc()
```

#### 6. Cardiovascular_disease

```

library(readr)
setwd("/share/pub/mayl/Sherlock/00_IEU_GWAS/06_Cardiovascular_disease")
############
#
a<-c(DBP_EastAsian="./DBP_EastAsian/bbj-a-17.vcf.gz",
DBP_European="./DBP_European/ieu-b-39.vcf.gz",
PR_EastAsian="./PR_EastAsian/bbj-a-46.vcf.gz",
SBP="./SBP/ieu-b-38.vcf.gz"
)

##开始循环
lapply(names(a),function(x){
GWAS_raw <-read_table(a[x],skip=119)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
  #ES:SE:LP:AF:ID
  value<-value[c("ES","SE","LP","AF")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF")
 gwas_data<-GWAS_raw[,c(1:3,6)]
 gwas_data<- cbind(gwas_data,value_df)
  gwas_data$LP <-  10^(-gwas_data$LP)
 colnames(gwas_data)<-c("chrom","pos","rsid","N","beta","se","p","maf")
 if(is.na(gwas_data$maf[1])) gwas_data$maf<-0.1
 
write.table(gwas_data,file=paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",x,"_gwas_data.txt"),row.names=F,quote=F)
})
#########另一种类型的数据读取
a<-c(CAD="./CAD/ebi-a-GCST005195.vcf.gz"

#T2D="./T2D/ebi-a-GCST006867.vcf.gz"#T2D数据很奇怪，很多问号不能用
)

#lapply(names(a),function(x){
GWAS_raw <-read_table(a[1],skip=120)
#  `1` = col_double(),
#  `754503` = col_double(),
#  rs3115859 = col_character(),
#  G = col_character(),
#  A = col_character(),
#  . = col_character(),
#  PASS = col_character(),
#  ._1 = col_character(),
#  `ES:SE:LP:ID` = col_character(),
#  `-0.000319708:0.00949475:0.0132283:rs3115859` = col_character()

colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","maf","FORMAT","IEU")
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
  #ES:SE:LP:AF:ID
  value<-value[c("ES","SE","LP","ID")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","ID")
 GWAS_raw$maf<-unlist(lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$maf[i], split = "=",fixed=T)
  index<-unlist(index)[2]
  return(index)
 }))
 
 gwas_data<-GWAS_raw[,c(1:3,6,8)]
 gwas_data<- cbind(gwas_data,value_df[,1:3])
  gwas_data$LP <-  10^(-gwas_data$LP)
 colnames(gwas_data)<-c("chrom","pos","rsid","N","maf","beta","se","p")
 if(is.na(gwas_data$maf[1])) gwas_data$maf<-0.1
 
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/CAD_gwas_data.txt",row.names=F,quote=F)
#})
rm(list=ls())
gc()

#"MDD" "ADHD" "ASD" "BIB" "EP" "Focal_EP" "Generalized_EP" "Juvenile_EP" "Migraine" "Multiple_sclerosis" "smokingCessation" "OpennessToExperence" "Happiness" "CognitivePerformance" "SubjectWellBeing" "BMI" "WHR" "WC" "Obesity" "CKD" "TC" "LDL" "HDL" "IBD" "PBC" "CoeliacDisease" "CrohnDisease" "T1D" "ATD" "DBP_EastAsian" "DBP_European" "PR_EastAsian" "SBP" "CAD"
```

#### 查看所有数据格式

```
cd /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata

for i in `ls`
do
echo $i
sed -n 1,4p $i
done

```

下面三个数据缺失了beta和se：

T1D_gwas_data.txt换成TIDnature的数据。
CrohnDisease_gwas_data.txt去掉
ATD_gwas_data.txt去掉

#很多数据没有maf这个数据，我们虚构了一列maf

# PLINK 软件进行SNP pruning

基于1000 Genome project欧洲人群样本生成plink 格式文件，用于后面的SNP pruning

#地址：/share/pub/mayl/Sherlock/plink
#修改PATH环境
#vim ~/.bashrc 
#//在最后一行添上：
#export PATH="/share/pub/mayl/Sherlock/plink":$PATH
#输入“source ~/.bashrc”命令，立即生效

## 第一部分
#### Refence LD based on the 1000Genome_European_Phase3
```
#地址：/share/pub/mayl/ABC_model/ldsc/data_ldsc/02_partitioned_LD_score_estimation/1000G_EUR_Phase3_plink
/share/pub/mayl/ABC_model/ldsc/data_ldsc/02_partitioned_LD_score_estimation/1000G_EUR_Phase3_plink
#Chr1-22
1000G.EUR.QC.9.bed
1000G.EUR.QC.9.bim
1000G.EUR.QC.9.fam

#The number of SNPs in each chromosome for this European reference
   510501 1000G.EUR.QC.10.bim
   493922 1000G.EUR.QC.11.bim
   480110 1000G.EUR.QC.12.bim
   366200 1000G.EUR.QC.13.bim
   324698 1000G.EUR.QC.14.bim
   287001 1000G.EUR.QC.15.bim
   316981 1000G.EUR.QC.16.bim
   269222 1000G.EUR.QC.17.bim
   285156 1000G.EUR.QC.18.bim
   232363 1000G.EUR.QC.19.bim
   779354 1000G.EUR.QC.1.bim
   221626 1000G.EUR.QC.20.bim
   138712 1000G.EUR.QC.21.bim
   141123 1000G.EUR.QC.22.bim
   839590 1000G.EUR.QC.2.bim
   706350 1000G.EUR.QC.3.bim
   729645 1000G.EUR.QC.4.bim
   633015 1000G.EUR.QC.5.bim
   664016 1000G.EUR.QC.6.bim
   589569 1000G.EUR.QC.7.bim
   549971 1000G.EUR.QC.8.bim
   438106 1000G.EUR.QC.9.bim
  9997231 total

#匹配22号染色体所在的行：
sed -n '/^22/p' /share/pub/jiangdp/Pagwas/GWAS/04_Metabolic_disease/WHR/WHR.txt
```



## 第二部分

得到所有表型的snp数据

```
#GWAS Summary statistics file
cd /share/pub/dengcy/GWAS_Multiomics/gwasdata
mkdir tempfile
cd ./predata

for i in MDD ADHD ASD PD BID EP Focal_EP Generalized_EP Juvenile_EP Migraine Multiple_sclerosis smokingCessation OpennessToExperence Happiness CognitivePerformance SubjectWellBeing BMI WHR WC Obesity TC LDL HDL IBD PBC CoeliacDisease T1D DBP_EastAsian DBP_European PR_EastAsian SBP CAD T2D UlcerativeColitis AD
do
awk  '{print $3 }' ${i}_gwas_data.txt  > /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${i}_SNP_list.txt
done
```



## 第三部分

Extract a subset of SNPs: "file-list" option
To extract only a subset of SNPs, it is possible to specify a list of required SNPs and make a new file, or perform an analysis on this subset, by using the command

```
plink --file data --extract mysnps.txt
```

#where the file is just a list of SNPs, one per line, e.g.
     snp005
     snp008
     snp101

#### 实际操作

在服务器上分别计算：

```
#!/usr/bin/env bash
#PBS -N plink1
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe

for j in MDD ADHD ASD PD
do
echo $j
for i in $(seq 1 22)  
do 
echo $i
plink --bfile /share/pub/mayl/ABC_model/ldsc/data_ldsc/02_partitioned_LD_score_estimation/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i --extract /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_SNP_list.txt --noweb --make-bed --out /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered
done 
done

#!/usr/bin/env bash
#PBS -N plink2
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
for j in BID EP Focal_EP Generalized_EP
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/mayl/ABC_model/ldsc/data_ldsc/02_partitioned_LD_score_estimation/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i --extract /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_SNP_list.txt --noweb --make-bed --out /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered
done 
done

#!/usr/bin/env bash
#PBS -N plink3
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe

for j in Juvenile_EP Migraine Multiple_sclerosis smokingCessation
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/mayl/ABC_model/ldsc/data_ldsc/02_partitioned_LD_score_estimation/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i --extract /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_SNP_list.txt --noweb --make-bed --out /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered
done 
done

#!/usr/bin/env bash
#PBS -N plink4
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
for j in OpennessToExperence Happiness CognitivePerformance SubjectWellBeing
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/mayl/ABC_model/ldsc/data_ldsc/02_partitioned_LD_score_estimation/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i --extract /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_SNP_list.txt --noweb --make-bed --out /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered
done 
done


#!/usr/bin/env bash
#PBS -N plink5
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
for j in BMI WHR WC Obesity
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/mayl/ABC_model/ldsc/data_ldsc/02_partitioned_LD_score_estimation/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i --extract /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_SNP_list.txt --noweb --make-bed --out /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered
done 
done

#!/usr/bin/env bash
#PBS -N plink6
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
for j in TC LDL HDL IBD
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/mayl/ABC_model/ldsc/data_ldsc/02_partitioned_LD_score_estimation/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i --extract /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_SNP_list.txt --noweb --make-bed --out /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered
done 
done

#!/usr/bin/env bash
#PBS -N plink7
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
for j in PBC CoeliacDisease T1D DBP_EastAsian
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/mayl/ABC_model/ldsc/data_ldsc/02_partitioned_LD_score_estimation/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i --extract /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_SNP_list.txt --noweb --make-bed --out /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered
done 
done

#!/usr/bin/env bash
#PBS -N plink8
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
for j in DBP_European PR_EastAsian SBP CAD
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/mayl/ABC_model/ldsc/data_ldsc/02_partitioned_LD_score_estimation/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i --extract /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_SNP_list.txt --noweb --make-bed --out /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered
done 
done

#!/usr/bin/env bash
#PBS -N plink9
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
for j in T2D UlcerativeColitis AD
do
  for i in $(seq 1 22)  
  do   
  plink --bfile   /share/pub/mayl/ABC_model/ldsc/data_ldsc/02_partitioned_LD_score_estimation/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i --extract /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_SNP_list.txt --noweb --make-bed --out /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered
  done 
done

```

## 第四部分

### 1 The VIF pruning routine is performed:

```
nohup plink --file data --indep 50,5,2
```

will create files

- plink.prune.in
- plink.prune.out

Each is a simlpe list of SNP IDs; 
both these files can subsequently be specified as the argument for a --extract or --exclude command.

1. The parameters for --indep are: window size in SNPs (e.g. 50), the number of SNPs to shift the window at each step (e.g. 5), the VIF threshold. 
2. The VIF is 1/(1-R^2) where R^2 is the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously. 
3. That is, this considers the correlations between SNPs but also between linear combinations of SNPs. 
4. A VIF of 10 is often taken to represent near collinearity problems in standard multiple regression analyses (i.e. implies R^2 of 0.9). 
5. A VIF of 1 would imply that the SNP is completely independent of all other SNPs. Practically, values between 1.5 and 2 should probably be used; particularly in small samples, if this threshold is too low and/or the window size is too large, too many SNPs may be removed.

### 2 The second procedure is performed:
```
plink --file data --indep-pairwise 50 5 0.5
```

1. This generates the same output files as the first version; the only difference is that a simple pairwise threshold is used. The first two parameters (50 and 5) are the same as above (window size and step); the third parameter represents the r^2 threshold. 
2. Note: this represents the pairwise SNP-SNP metric now, not the multiple correlation coefficient; also note, this is based on the genotypic correlation, i.e. it does not involve phasing.
3. To give a concrete example: the command above that specifies 50 5 0.5 would:
     a) consider a window of 50 SNPs, 
     b) calculate LD between each pair of SNPs in the window, remove one of a pair of SNPs if the LD is greater than 0.5, 
     c) shift the window 5 SNPs forward and repeat the procedure.

### 3 To make a new, pruned file, then use something like (in this example, we also convert the standard PED fileset to a binary one):
```
plink --file data --extract plink.prune.in --make-bed --out pruneddata
```



#### 实际操作

```
#!/usr/bin/env bash
#PBS -N plink1
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe

for j in MDD ADHD ASD
do
echo $j
for i in $(seq 1 22)  
do 
echo $i
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

#!/usr/bin/env bash
#PBS -N plink2
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
for j in BID EP Focal_EP
do
for i in $(seq 1 22)
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

#!/usr/bin/env bash
#PBS -N plink3
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe

for j in Juvenile_EP Migraine Multiple_sclerosis
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

#!/usr/bin/env bash
#PBS -N plink4
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
for j in OpennessToExperence Happiness CognitivePerformance
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done


#!/usr/bin/env bash
#PBS -N plink5
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
for j in BMI WHR WC
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

#!/usr/bin/env bash
#PBS -N plink6
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
for j in TC LDL HDL
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

#!/usr/bin/env bash
#PBS -N plink7
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
for j in PBC CoeliacDisease T1D
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

#!/usr/bin/env bash
#PBS -N plink8
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
for j in DBP_European PR_EastAsian SBP 
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

#!/usr/bin/env bash
#PBS -N plink9
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
for j in T2D UlcerativeColitis AD
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

#!/usr/bin/env bash
#PBS -N plink10
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
for j in CAD DBP_EastAsian IBD
do
for i in $(seq 1 22)  
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done


#!/usr/bin/env bash
#PBS -N plink11
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
for j in Obesity SubjectWellBeing smokingCessation
do
for i in $(seq 1 22)
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done


#!/usr/bin/env bash
#PBS -N plink12
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
for j in Generalized_EP PD
do
for i in $(seq 1 22)
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done


#!/usr/bin/env bash
#PBS -N plink1
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe

for i in $(seq 1 3)
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.PD_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/PD_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

#!/usr/bin/env bash
#PBS -N plink2
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe

for i in $(seq 4 6)
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.PD_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/PD_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

#!/usr/bin/env bash
#PBS -N plink3
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe

for i in $(seq 7 9)
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.PD_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/PD_${i}_plink_prune_EUR_filtered_LD0.8
done 


#!/usr/bin/env bash
#PBS -N plink4
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe

for i in $(seq 10 12)
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.PD_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/PD_${i}_plink_prune_EUR_filtered_LD0.8
done 




#!/usr/bin/env bash
#PBS -N plink5
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe

for i in $(seq 13 15)
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.PD_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/PD_${i}_plink_prune_EUR_filtered_LD0.8
done


#!/usr/bin/env bash
#PBS -N plink6
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe

for i in $(seq 16 18)
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.PD_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/PD_${i}_plink_prune_EUR_filtered_LD0.8
done


#!/usr/bin/env bash
#PBS -N plink7
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe

for i in $(seq 19 22)
do   
plink --bfile /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/1000G.EUR.QC.PD_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/PD_${i}_plink_prune_EUR_filtered_LD0.8
done 

```

## 第五部分

合并不同表型的SNP：

```
###合并SNPs
cd /share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/
for j in MDD ADHD ASD PD BID EP Focal_EP Generalized_EP Juvenile_EP Migraine Multiple_sclerosis smokingCessation OpennessToExperence Happiness CognitivePerformance SubjectWellBeing BMI WHR WC Obesity TC LDL HDL IBD PBC CoeliacDisease T1D DBP_EastAsian DBP_European PR_EastAsian SBP CAD T2D UlcerativeColitis AD
do
echo $j
cat [${j}]*.prune.in > ${j}_merge_plink_EUR_filtered_LD0.8.prune.in
done 
#####删除多余文件
rm -f {*.prune.in}
##文件剩余的行数
for j in MDD ADHD ASD PD BID EP Focal_EP Generalized_EP Juvenile_EP Migraine Multiple_sclerosis smokingCessation OpennessToExperence Happiness CognitivePerformance SubjectWellBeing BMI WHR WC Obesity TC LDL HDL IBD PBC CoeliacDisease T1D DBP_EastAsian DBP_European PR_EastAsian SBP CAD T2D UlcerativeColitis AD
do
wc -l ${j}_merge_plink_EUR_filtered_LD0.8.prune.out
wc -l /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/${j}_gwas_data.txt
done 

###########################文件剩余的行数（snp的数量）
23100115 MDD_merge_plink_EUR_filtered_LD0.8.prune.out
7338291 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/MDD_gwas_data.txt
32176215 ADHD_merge_plink_EUR_filtered_LD0.8.prune.out
7468462 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/ADHD_gwas_data.txt
63357810 ASD_merge_plink_EUR_filtered_LD0.8.prune.out
8370011 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/ASD_gwas_data.txt
19088919 PD_merge_plink_EUR_filtered_LD0.8.prune.out
17891935 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/PD_gwas_data.txt
25230627 BID_merge_plink_EUR_filtered_LD0.8.prune.out
12382416 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/BID_gwas_data.txt
32835373 EP_merge_plink_EUR_filtered_LD0.8.prune.out
4880491 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/EP_gwas_data.txt
69341006 Focal_EP_merge_plink_EUR_filtered_LD0.8.prune.out
4862781 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/Focal_EP_gwas_data.txt
69341267 Generalized_EP_merge_plink_EUR_filtered_LD0.8.prune.out
4867067 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/Generalized_EP_gwas_data.txt
69423834 Juvenile_EP_merge_plink_EUR_filtered_LD0.8.prune.out
4983224 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/Juvenile_EP_gwas_data.txt
37169388 Migraine_merge_plink_EUR_filtered_LD0.8.prune.out
6386694 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/Migraine_gwas_data.txt
78647570 Multiple_sclerosis_merge_plink_EUR_filtered_LD0.8.prune.out
6304357 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/Multiple_sclerosis_gwas_data.txt
15131768 smokingCessation_merge_plink_EUR_filtered_LD0.8.prune.out
5961481 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/smokingCessation_gwas_data.txt
65179717 OpennessToExperence_merge_plink_EUR_filtered_LD0.8.prune.out
2302978 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/OpennessToExperence_gwas_data.txt
26591761 Happiness_merge_plink_EUR_filtered_LD0.8.prune.out
9329578 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/Happiness_gwas_data.txt
39969970 CognitivePerformance_merge_plink_EUR_filtered_LD0.8.prune.out
10049917 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/CognitivePerformance_gwas_data.txt
45727421 SubjectWellBeing_merge_plink_EUR_filtered_LD0.8.prune.out
2053172 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/SubjectWellBeing_gwas_data.txt
194416758 BMI_merge_plink_EUR_filtered_LD0.8.prune.out
9851865 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/BMI_gwas_data.txt
37283386 WHR_merge_plink_EUR_filtered_LD0.8.prune.out
2529832 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/WHR_gwas_data.txt
91616756 WC_merge_plink_EUR_filtered_LD0.8.prune.out
2533038 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/WC_gwas_data.txt
87957838 Obesity_merge_plink_EUR_filtered_LD0.8.prune.out
2377583 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/Obesity_gwas_data.txt
56670790 TC_merge_plink_EUR_filtered_LD0.8.prune.out
2423571 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/TC_gwas_data.txt
10732104 LDL_merge_plink_EUR_filtered_LD0.8.prune.out
2414489 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/LDL_gwas_data.txt
55207168 HDL_merge_plink_EUR_filtered_LD0.8.prune.out
2424036 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/HDL_gwas_data.txt
244878012 IBD_merge_plink_EUR_filtered_LD0.8.prune.out
9603981 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/IBD_gwas_data.txt
310387114 PBC_merge_plink_EUR_filtered_LD0.8.prune.out
1121258 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/PBC_gwas_data.txt
79264348 CoeliacDisease_merge_plink_EUR_filtered_LD0.8.prune.out
523399 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/CoeliacDisease_gwas_data.txt
71579478 T1D_merge_plink_EUR_filtered_LD0.8.prune.out
62115238 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/T1D_gwas_data.txt
746504637 DBP_EastAsian_merge_plink_EUR_filtered_LD0.8.prune.out
6108954 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/DBP_EastAsian_gwas_data.txt
1362040513 DBP_European_merge_plink_EUR_filtered_LD0.8.prune.out
7160618 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/DBP_European_gwas_data.txt
507026621 PR_EastAsian_merge_plink_EUR_filtered_LD0.8.prune.out
6108954 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/PR_EastAsian_gwas_data.txt
1128891905 SBP_merge_plink_EUR_filtered_LD0.8.prune.out
7088082 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/SBP_gwas_data.txt
2359161483 CAD_merge_plink_EUR_filtered_LD0.8.prune.out
7918806 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/CAD_gwas_data.txt
2251704106 T2D_merge_plink_EUR_filtered_LD0.8.prune.out
5021051 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/T2D_gwas_data.txt
2514053491 UlcerativeColitis_merge_plink_EUR_filtered_LD0.8.prune.out
9459542 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/UlcerativeColitis_gwas_data.txt
2229104191 AD_merge_plink_EUR_filtered_LD0.8.prune.out
10528609 /share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/AD_gwas_data.txt

```

将过滤后的snp与之前得到的gwas数据进行合并：打开R 

```
library(readr)
library(dplyr)

for (i in c("MDD","ADHD","ASD","PD","BID","EP","Focal_EP","Generalized_EP","Juvenile_EP","Migraine","Multiple_sclerosis","smokingCessation","OpennessToExperence","Happiness","CognitivePerformance","DBP_EastAsian","SubjectWellBeing","BMI","WHR","WC","Obesity","TC","LDL","HDL","IBD","PBC","CoeliacDisease","T1D","DBP_European","PR_EastAsian","SBP","CAD","T2D","UlcerativeColitis","AD")) {
    gwas<-read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",i,"_gwas_data.txt"))
    SNP_prune<- read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/",i,"_merge_plink_EUR_filtered_LD0.8.prune.in"))
    SNP_prune<-SNP_prune[!duplicated(unlist(SNP_prune)),]
    colnames(SNP_prune)<-"rsid"
    #### Left Join using inner_join function 
   gwas= gwas %>% inner_join(SNP_prune,by="rsid")
        write.table(gwas,file= paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/prunedata/",i,"_prune_gwas_data.txt"),,row.names=F,quote=F)
        print(i)
gc()
}

#for (i in c("DBP_European","PR_EastAsian","SBP","CAD","AD","T2D","UlcerativeColitis","DBP_EastAsian")) {

i<-"DBP_European"
    gwas<-read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",i,"_gwas_data.txt"))
    SNP_prune<- read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/",i,"_merge_plink_EUR_filtered_LD0.8.prune.in"))
    SNP_prune<-SNP_prune[!duplicated(unlist(SNP_prune)),]
    colnames(SNP_prune)<-"rsid"
    #### Left Join using inner_join function 
   gwas= gwas %>% inner_join(SNP_prune,by="rsid")
        write.table(gwas,file= paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/prunedata/",i,"_prune_gwas_data.txt"),,row.names=F,quote=F)
        print(i)
gc()
#}
###########################################################

#!/usr/bin/env bash
#PBS -N 8
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
source activate R403
Rscript /share/pub/dengcy/GWAS_Multiomics/gwasdata/8.r

library(readr)
library(dplyr)
i<-"DBP_European"
    gwas<-read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",i,"_gwas_data.txt"))
    SNP_prune<- read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/",i,"_merge_plink_EUR_filtered_LD0.8.prune.in"))
    SNP_prune<-SNP_prune[!duplicated(unlist(SNP_prune)),]
    colnames(SNP_prune)<-"rsid"
    #### Left Join using inner_join function 
   gwas= gwas %>% inner_join(SNP_prune,by="rsid")
        write.table(gwas,file= paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/prunedata/",i,"_prune_gwas_data.txt"),,row.names=F,quote=F)
        print(i)
gc()

#!/usr/bin/env bash
#PBS -N 7
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
source activate R403
Rscript /share/pub/dengcy/GWAS_Multiomics/gwasdata/7.r

library(readr)
library(dplyr)
i<-"PR_EastAsian"
    gwas<-read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",i,"_gwas_data.txt"))
    SNP_prune<- read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/",i,"_merge_plink_EUR_filtered_LD0.8.prune.in"))
    SNP_prune<-SNP_prune[!duplicated(unlist(SNP_prune)),]
    colnames(SNP_prune)<-"rsid"
    #### Left Join using inner_join function 
   gwas= gwas %>% inner_join(SNP_prune,by="rsid")
        write.table(gwas,file= paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/prunedata/",i,"_prune_gwas_data.txt"),,row.names=F,quote=F)
        print(i)
gc()

#!/usr/bin/env bash
#PBS -N 6
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
source activate R403
Rscript /share/pub/dengcy/GWAS_Multiomics/gwasdata/6.r

library(readr)
library(dplyr)
i<-"SBP"
    gwas<-read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",i,"_gwas_data.txt"))
    SNP_prune<- read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/",i,"_merge_plink_EUR_filtered_LD0.8.prune.in"))
    SNP_prune<-SNP_prune[!duplicated(unlist(SNP_prune)),]
    colnames(SNP_prune)<-"rsid"
    #### Left Join using inner_join function 
   gwas= gwas %>% inner_join(SNP_prune,by="rsid")
        write.table(gwas,file= paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/prunedata/",i,"_prune_gwas_data.txt"),,row.names=F,quote=F)
        print(i)
gc()


#!/usr/bin/env bash
#PBS -N 5
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
source activate R403
Rscript /share/pub/dengcy/GWAS_Multiomics/gwasdata/5.r

library(readr)
library(dplyr)
i<-"CAD"
    gwas<-read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",i,"_gwas_data.txt"))
    SNP_prune<- read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/",i,"_merge_plink_EUR_filtered_LD0.8.prune.in"))
    SNP_prune<-SNP_prune[!duplicated(unlist(SNP_prune)),]
    colnames(SNP_prune)<-"rsid"
    #### Left Join using inner_join function 
   gwas= gwas %>% inner_join(SNP_prune,by="rsid")
        write.table(gwas,file= paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/prunedata/",i,"_prune_gwas_data.txt"),,row.names=F,quote=F)
        print(i)
gc()


#!/usr/bin/env bash
#PBS -N 4
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
source activate R403
Rscript /share/pub/dengcy/GWAS_Multiomics/gwasdata/4.r

library(readr)
library(dplyr)
i<-"AD"
    gwas<-read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",i,"_gwas_data.txt"))
    SNP_prune<- read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/",i,"_merge_plink_EUR_filtered_LD0.8.prune.in"))
    SNP_prune<-SNP_prune[!duplicated(unlist(SNP_prune)),]
    colnames(SNP_prune)<-"rsid"
    #### Left Join using inner_join function 
   gwas= gwas %>% inner_join(SNP_prune,by="rsid")
        write.table(gwas,file= paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/prunedata/",i,"_prune_gwas_data.txt"),,row.names=F,quote=F)
        print(i)
gc()

#!/usr/bin/env bash
#PBS -N 3
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
source activate R403
Rscript /share/pub/dengcy/GWAS_Multiomics/gwasdata/3.r

library(readr)
library(dplyr)
i<-"T2D"
    gwas<-read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",i,"_gwas_data.txt"))
    SNP_prune<- read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/",i,"_merge_plink_EUR_filtered_LD0.8.prune.in"))
    SNP_prune<-SNP_prune[!duplicated(unlist(SNP_prune)),]
    colnames(SNP_prune)<-"rsid"
    #### Left Join using inner_join function 
   gwas= gwas %>% inner_join(SNP_prune,by="rsid")
        write.table(gwas,file= paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/prunedata/",i,"_prune_gwas_data.txt"),,row.names=F,quote=F)
        print(i)
gc()


#!/usr/bin/env bash
#PBS -N 2
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
source activate R403
Rscript /share/pub/dengcy/GWAS_Multiomics/gwasdata/2.r

library(readr)
library(dplyr)
i<-"UlcerativeColitis"
    gwas<-read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",i,"_gwas_data.txt"))
    SNP_prune<- read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/",i,"_merge_plink_EUR_filtered_LD0.8.prune.in"))
    SNP_prune<-SNP_prune[!duplicated(unlist(SNP_prune)),]
    colnames(SNP_prune)<-"rsid"
    #### Left Join using inner_join function 
   gwas= gwas %>% inner_join(SNP_prune,by="rsid")
        write.table(gwas,file= paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/prunedata/",i,"_prune_gwas_data.txt"),,row.names=F,quote=F)
        print(i)
gc()


#!/usr/bin/env bash
#PBS -N 1
#PBS -q workq
#PBS -l nodes=node03
#PBS -l ncpus=1
#PBS -j oe
source activate R403
Rscript /share/pub/dengcy/GWAS_Multiomics/gwasdata/1.r


library(readr)
library(dplyr)
i<-"DBP_EastAsian"
    gwas<-read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/predata/",i,"_gwas_data.txt"))
    SNP_prune<- read_table(paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/tempfile/",i,"_merge_plink_EUR_filtered_LD0.8.prune.in"))
    SNP_prune<-SNP_prune[!duplicated(unlist(SNP_prune)),]
    colnames(SNP_prune)<-"rsid"
    #### Left Join using inner_join function 
   gwas= gwas %>% inner_join(SNP_prune,by="rsid")
        write.table(gwas,file= paste0("/share/pub/dengcy/GWAS_Multiomics/gwasdata/prunedata/",i,"_prune_gwas_data.txt"),,row.names=F,quote=F)
        print(i)
gc()

```

## blood traits gwas 数据的处理

### RBCcount

```R
library(readr)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits")
############
RBCcount<-"ieu-a-275.vcf.gz"

GWAS_raw <-read_table(RBCcount,skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
    
  value<-value[c("ES","SE","LP","SS")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","SS")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$LP <-  10^(-gwas_data$LP)
#CHROM    POS         ID REF ALT      ES     SE       LP    SS
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","N")
# if(is.na(gwas_data$maf[1])) 
     gwas_data$maf<-0.1
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/RBCcount_gwas_data.txt",row.names=F,quote=F)
gwas_data<-read.table("RBCcount_gwas_data.txt",header=T)
gwas_data$p[which(gwas_data$p>1)]<-1

######################ukb-d-30010_irnt
library(readr)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits")
############
RBCcount<-"ukb-d-30010_irnt.vcf.gz"
GWAS_raw <-read_table(RBCcount,skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$INFO[i], split = "=",fixed=T)[[1]][[2]]
  return(index)
  })
GWAS_raw$INFO<-unlist(value_list)

value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
    
  value<-value[c("ES","SE","LP","SS")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","SS")
 gwas_data<-GWAS_raw[,c(1:5,8)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$LP <-  10^(-gwas_data$LP)
#CHROM    POS         ID REF ALT      ES     SE       LP    SS
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","maf","beta","se","p","N")
# if(is.na(gwas_data$maf[1])) 
     #gwas_data$maf<-0.1
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/RBCcount2_gwas_data.txt",row.names=F,quote=F)
gwas_data<-read.table("RBCcount_gwas_data.txt",header=T)
gwas_data$p[which(gwas_data$p>1)]<-1

```

###  Plateletcount

```R

library(readr)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits")
############
Plateletcount<-"ieu-a-1008.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(Plateletcount,skip=109)
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
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/Plateletcount_gwas_data.txt",row.names=F,quote=F)
```

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

### lymphocyte cell count

Dataset: ieu-b-32

```R
library(readr)
setwd("/home/guofj/bloodtraits")
############
lymphocytecount<-"ieu-b-32.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table2(lymphocytecount,comment = "#")

colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:ID
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
write.table(gwas_data,file="lymphocytecount_gwas_data.txt",row.names=F,quote=F)
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
write.table(gwas_data,file="WhiteBloodCellcount_gwas_data.txt",row.names=F,quote=F)
```

### Lymphocyte percentage of white cells

Dataset: ebi-a-GCST004632

```R
library(readr)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits")
############
LymphocytePercent<-"ebi-a-GCST004632.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(LymphocytePercent,comment = "#")

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
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/LymphocytePercent_gwas_data.txt",row.names=F,quote=F)
```



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



### **Lymphocyte percentage**

Dataset: **Dataset: ukb-d-30180_irnt**

```
library(readr)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits")
############
MeanCorpusVolume<-"ukb-d-30180_irnt.vcf.gz"
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
write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/Lymphocytepercentage2_gwas_data.txt",row.names=F,quote=F)
```

**Lymphocyte count**
**Dataset: ukb-d-30120_irnt**

```
library(readr)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits")
############
Lymphocyte_count<-"ukb-d-30120_irnt.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(Lymphocyte_count,comment = "#")

colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  #index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  #index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)[1:5]
  names(value)<-c("ES","SE","LP","AF","SS")
    #ES:SE:LP:AF:SS:ID
  #value<-value[c("ES","SE","LP","AF","SS")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=5)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF","SS")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$p <-  10^(-gwas_data$LP)
 colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
gwas_data<-gwas_data[,c("CHROM", "POS","ID","REF","ALT","ES","SE","p","AF","SS")]
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf","N")

write.table(gwas_data,file="/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/Lymphocytecount2_gwas_data.txt",row.names=F,quote=F)
```

**Lymphocyte count**3
**Dataset: bbj-a-36**

```
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



## blood traits snp数据的pruning

```shell
cd /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits
mkdir tempfile
#######################1.
#neutrophilcount WhiteBloodCellcount MeanCorpusVolume 
#Plateletcount basophilcount eosinophilcount lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount MeanCorpusVolume 
###最新数据
for i in RBCcount2 Lymphocytepercentage2 
do
awk  '{print $3 }' ${i}_gwas_data.txt  > /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/${i}_SNP_list.txt
done

for i in LymphocytePercent  MeanCorpuscularHemoglobin Hemoglobinconcen
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
###
#LymphocytePercent  MeanCorpuscularHemoglobin
for j in Hemoglobinconcen
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
for j in RBCcount2 Lymphocytepercentage2
do
echo $j
for i in $(seq 1 22)  
do 
echo $i
plink --bfile /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

for j in LymphocytePercent  MeanCorpuscularHemoglobin Hemoglobinconcen
do
echo $j
for i in $(seq 1 22)  
do 
echo $i
plink --bfile /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

for j in RBCcount Plateletcount basophilcount 
do
echo $j
for i in $(seq 1 22)  
do 
echo $i
plink --bfile /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

for j in neutrophilcount WhiteBloodCellcount MeanCorpusVolume
do
echo $j
for i in $(seq 1 22)  
do 
echo $i
plink --bfile /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

##########上面的没跑完，需要再跑三个：
for j in MeanCorpusVolume
do
echo $j
for i in $(seq 1 22)  
do 
echo $i
plink --bfile /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

for j in monocytecount
do
echo $j
for i in $(seq 1 22)  
do 
echo $i
plink --bfile /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

for j in Hemoglobinconcen
do
echo $j
for i in $(seq 1 22)  
do 
echo $i
plink --bfile /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/1000G.EUR.QC.${j}_${i}_filtered --indep-pairwise 50 5 0.8 --out  /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
done

for j in MeanCorpusVolume
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
#####删除多余文件
rm -f {*.prune.out}
##文件剩余的行数
for j in RBCcount Plateletcount basophilcount eosinophilcount lymphocytecount monocytecount neutrophilcount WhiteBloodCellcount LymphocytePercent Hemoglobinconcen MeanCorpuscularHemoglobin MeanCorpusVolume
do
echo $j
wc -l ${j}_merge_plink_EUR_filtered_LD0.8.prune.in
wc -l /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/${j}_gwas_data.txt
done 
```

## blood traits snp数据的整合

```R
library(readr)
library(dplyr)
#for j in 
#"RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume"

for (i in c("RBCcount2","Lymphocytepercentage2")) {
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

###########################################################

#!/usr/bin/env bash
#PBS -N 8
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=1
#PBS -j oe
source activate R403
Rscript /share/pub/dengcy/GWAS_Multiomics/gwasdata/8.r


```
