########################################
#Get the random test gwas data
########################
setwd("C:\\Users\\MYL\\Desktop")
setwd(".\\03-Simulation-analysis\\GWAS_simulation_100SNP")
data_SCZ <- read.table("PGC3_SCZ_wave3_public_Pagwas",header = T)
data_SCZ_New <- data_SCZ[,c(1,3,2,9,10,6)]  #order the data
len <- length(data_SCZ_New[,1])
index <- sample(1:len,1000000,replace = F)
data_SCZ_New_100W <- data_SCZ_New[index,]
colnames(data_SCZ_New_100W) <- c("chrom","pos","rsid","beta","se","maf")
##chrom       pos       rsid        beta        se     maf

length(data_SCZ_New_100W$chrom)

##########################Çé¿öÒ»
set.seed(1234)
h=0.5
p=0.1
num = 1000000
sd<-sqrt(h^2/(num*p))
y <- rnorm(num,mean = 0, sd = sd)
plot(density(y))


#1 a = 2, p = 0.1
data_SCZ_New_100W_p0.1_a2 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.1_a2$beta<-y
a=2
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.1_a2[i,6]
  data_SCZ_New_100W_p0.1_a2[i,4]=data_SCZ_New_100W_p0.1_a2[i,4]*(2*maf*(1-maf))^a
}
data2_p0.1_a2 <-data_SCZ_New_100W_p0.1_a2$beta
plot(density(data2_p0.1_a2))

save(data_SCZ_New_100W_p0.1_a2,file="data_SCZ_New_100W_h0.5_p_0.1_a2.Rdata")


#2 a = 0.72, p = 0.1
data_SCZ_New_100W_p0.1_a0.72 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.1_a0.72$beta<-y
a=0.72
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.1_a0.72[i,6]
  data_SCZ_New_100W_p0.1_a0.72[i,4]=data_SCZ_New_100W_p0.1_a0.72[i,4]*(2*maf*(1-maf))^a
}
data2_p0.1_a0.72 <-data_SCZ_New_100W_p0.1_a0.72$beta
plot(density(data2_p0.1_a0.72))

save(data_SCZ_New_100W_p0.1_a0.72,file="data_SCZ_New_100W_h0.5_p_0.1_a0.72.Rdata")



#3 a = 1, p = 0.1
data_SCZ_New_100W_p0.1_a1 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.1_a1$beta<-y
a=1
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.1_a1[i,6]
  data_SCZ_New_100W_p0.1_a1[i,4]=data_SCZ_New_100W_p0.1_a1[i,4]*(2*maf*(1-maf))^a
}
data2_p0.1_a1 <-data_SCZ_New_100W_p0.1_a1$beta
plot(density(data2_p0.1_a1))

save(data_SCZ_New_100W_p0.1_a1,file="data_SCZ_New_100W_h0.5_p_0.1_a1.Rdata")


#4 a = -1, p = 0.1
data_SCZ_New_100W_p0.1_a_1 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.1_a_1$beta<-y
a=-1
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.1_a_1[i,6]
  data_SCZ_New_100W_p0.1_a_1[i,4]=data_SCZ_New_100W_p0.1_a_1[i,4]*(2*maf*(1-maf))^a
}
data2_p0.1_a_1 <-data_SCZ_New_100W_p0.1_a_1$beta
plot(density(data2_p0.1_a_1))

save(data_SCZ_New_100W_p0.1_a_1,file="data_SCZ_New_100W_h0.5_p_0.1_a_1.Rdata")



#5  a = -2, p = 0.1
data_SCZ_New_100W_p0.1_a_2 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.1_a_2$beta<-y
a=-2
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.1_a_2[i,6]
  data_SCZ_New_100W_p0.1_a_2[i,4]=data_SCZ_New_100W_p0.1_a_2[i,4]*(2*maf*(1-maf))^a
}
data2_p0.1_a_2 <-data_SCZ_New_100W_p0.1_a_2$beta
plot(density(data2_p0.1_a_2))

save(data_SCZ_New_100W_p0.1_a_2,file="data_SCZ_New_100W_h0.5_p_0.1_a_2.Rdata")


#6 a = 0, p = 0.1
data_SCZ_New_100W_p0.1_a_0 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.1_a_0$beta<-y

save(data_SCZ_New_100W_p0.1_a_0,file="data_SCZ_New_100W_h0.5_p_0.1_a_0.Rdata")


plot(density(data2_p0.1_a2))
plot(density(data2_p0.1_a1))
plot(density(data2_p0.1_a0.72))
plot(density(y))
plot(density(data2_p0.1_a_1))
plot(density(data2_p0.1_a_2))



plot(density(data2_p0.01_a2))
plot(density(data2_p0.01_a1))
plot(density(data2_p0.01_a0.72))
plot(density(y1))
plot(density(data2_p0.01_a_1))
plot(density(data2_p0.01_a_2))


#################################
##########################Çé¿ö¶þ
set.seed(1234)
h=0.5
p=0.01
num = 1000000
sd<-sqrt(h^2/(num*p))
y1 <- rnorm(num,mean = 0, sd = sd )
plot(density(y1))

#sd<-sqrt(h^2/(num*p))
y2 <- rnorm(num,mean = 0.5, sd = sd)

#1 a = 2, p = 0.01
data_SCZ_New_100W_p0.01_a2 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.01_a2$beta<-y1
data_SCZ_New_100W_p0.01_a2$p<-y2

se=sqrt(((data_SCZ_New_100W_p0.01_a2$beta)^2)/qchisq(data_SCZ_New_100W_p0.01_a2$p,1,lower.tail=F))
data_SCZ_New_100W_p0.01_a2$p<-se


a=2
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.01_a2[i,6]
  data_SCZ_New_100W_p0.01_a2[i,4]=data_SCZ_New_100W_p0.01_a2[i,4]*(2*maf*(1-maf))^a
}
#data2_p0.01_a2 <-data_SCZ_New_100W_p0.01_a2$beta
#plot(density(data2_p0.01_a2))

save(data_SCZ_New_100W_p0.01_a2,file="data_SCZ_New_100W_h0.5_p_0.01_a2.Rdata")

#2 a = 0.72, p = 0.01
data_SCZ_New_100W_p0.01_a0.72 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.01_a0.72$beta<-y1
a=0.72
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.01_a0.72[i,6]
  data_SCZ_New_100W_p0.01_a0.72[i,4]=data_SCZ_New_100W_p0.01_a0.72[i,4]*(2*maf*(1-maf))^a
}
data2_p0.01_a0.72 <-data_SCZ_New_100W_p0.01_a0.72$beta
plot(density(data2_p0.01_a0.72))

save(data_SCZ_New_100W_p0.01_a0.72,file="data_SCZ_New_100W_h0.5_p_0.01_a0.72.Rdata")

#3 a = 1, p = 0.01
data_SCZ_New_100W_p0.01_a1 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.01_a1$beta<-y1
a=1
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.01_a1[i,6]
  data_SCZ_New_100W_p0.01_a1[i,4]=data_SCZ_New_100W_p0.01_a1[i,4]*(2*maf*(1-maf))^a
}
data2_p0.01_a1 <-data_SCZ_New_100W_p0.01_a1$beta
plot(density(data2_p0.01_a1))

save(data_SCZ_New_100W_p0.01_a1,file="data_SCZ_New_100W_h0.5_p_0.01_a1.Rdata")

#4 a = -1, p = 0.01
data_SCZ_New_100W_p0.01_a_1 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.01_a_1$beta<-y1
a=-1
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.01_a_1[i,6]
  data_SCZ_New_100W_p0.01_a_1[i,4]=data_SCZ_New_100W_p0.01_a_1[i,4]*(2*maf*(1-maf))^a
}
data2_p0.01_a_1 <-data_SCZ_New_100W_p0.01_a_1$beta
plot(density(data2_p0.01_a_1))

save(data_SCZ_New_100W_p0.01_a_1,file="data_SCZ_New_100W_h0.5_p_0.01_a-1.Rdata")


#5 a = -2, p = 0.01
data_SCZ_New_100W_p0.01_a_2 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.01_a_2$beta<-y1
a=-2
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.01_a_2[i,6]
  data_SCZ_New_100W_p0.01_a_2[i,4]=data_SCZ_New_100W_p0.01_a_2[i,4]*(2*maf*(1-maf))^a
}
data2_p0.01_a_2 <-data_SCZ_New_100W_p0.01_a_2$beta
plot(density(data2_p0.01_a_2))

save(data_SCZ_New_100W_p0.01_a_2,file="data_SCZ_New_100W_h0.5_p_0.01_a-2.Rdata")


#6 a = 0, p = 0.01
data_SCZ_New_100W_p0.01_a_0 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.01_a_0$beta<-y1

save(data_SCZ_New_100W_p0.01_a_0,file="data_SCZ_New_100W_h0.5_p_0.01_a0.Rdata")



##########################-----------------------------------------------------

setwd("C:\\Users\\MYL\\Desktop")
setwd("\\03-Simulation-analysis\\GWAS_simulation_100SNP")
data_SCZ <- read.table("PGC3_SCZ_wave3_public_Pagwas",header = T)
data_SCZ_New <- data_SCZ[,c(1,3,2,9,10,6)]  #order the data
len <- length(data_SCZ_New[,1])
index <- sample(1:len,1000000,replace = F)
data_SCZ_New_100W <- data_SCZ_New[index,]
colnames(data_SCZ_New_100W) <- c("chrom","pos","rsid","beta","se","maf")
##chrom       pos       rsid        beta        se     maf

length(data_SCZ_New_100W$chrom)

##########################Çé¿ö3
set.seed(1234)
h=0.5
p=0.001
num = 1000000
sd<-sqrt(h^2/(num*p))
y <- rnorm(num,mean = 0, sd = sd)
plot(density(y))


#1 a = 2, p = 0.001
data_SCZ_New_100W_p0.001_a2 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.001_a2$beta<-y
a=2
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.001_a2[i,6]
  data_SCZ_New_100W_p0.001_a2[i,4]=data_SCZ_New_100W_p0.001_a2[i,4]*(2*maf*(1-maf))^a
}
data2_p0.001_a2 <-data_SCZ_New_100W_p0.001_a2$beta
plot(density(data2_p0.001_a2))

save(data_SCZ_New_100W_p0.001_a2,file="data_SCZ_New_100W_h0.5_p_0.001_a2.Rdata")


#2 a = 0.72, p = 0.001
data_SCZ_New_100W_p0.001_a0.72 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.001_a0.72$beta<-y
a=0.72
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.001_a0.72[i,6]
  data_SCZ_New_100W_p0.001_a0.72[i,4]=data_SCZ_New_100W_p0.001_a0.72[i,4]*(2*maf*(1-maf))^a
}
data2_p0.001_a0.72 <-data_SCZ_New_100W_p0.001_a0.72$beta
plot(density(data2_p0.001_a0.72))

save(data_SCZ_New_100W_p0.001_a0.72,file="data_SCZ_New_100W_h0.5_p_0.001_a0.72.Rdata")



#3 a = 1, p = 0.001
data_SCZ_New_100W_p0.001_a1 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.001_a1$beta<-y
a=1
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.001_a1[i,6]
  data_SCZ_New_100W_p0.001_a1[i,4]=data_SCZ_New_100W_p0.001_a1[i,4]*(2*maf*(1-maf))^a
}
data2_p0.001_a1 <-data_SCZ_New_100W_p0.001_a1$beta
plot(density(data2_p0.001_a1))

save(data_SCZ_New_100W_p0.001_a1,file="data_SCZ_New_100W_h0.5_p_0.001_a1.Rdata")



#4 a = -1, p = 0.001
data_SCZ_New_100W_p0.001_a_1 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.001_a_1$beta<-y
a=-1
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.001_a_1[i,6]
  data_SCZ_New_100W_p0.001_a_1[i,4]=data_SCZ_New_100W_p0.001_a_1[i,4]*(2*maf*(1-maf))^a
}
data2_p0.001_a_1 <-data_SCZ_New_100W_p0.001_a_1$beta
plot(density(data2_p0.001_a_1))

save(data_SCZ_New_100W_p0.001_a_1,file="data_SCZ_New_100W_h0.5_p_0.001_a_1.Rdata")



#5 a = -2, p = 0.001
data_SCZ_New_100W_p0.001_a_2 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.001_a_2$beta<-y
a=-2
for (i in 1:num){
  maf<- data_SCZ_New_100W_p0.001_a_2[i,6]
  data_SCZ_New_100W_p0.001_a_2[i,4]=data_SCZ_New_100W_p0.001_a_2[i,4]*(2*maf*(1-maf))^a
}
data2_p0.001_a_2 <-data_SCZ_New_100W_p0.001_a_2$beta
plot(density(data2_p0.001_a_2))

save(data_SCZ_New_100W_p0.001_a_2,file="data_SCZ_New_100W_h0.5_p_0.001_a_2.Rdata")


#6 a = 0, p = 0.001
data_SCZ_New_100W_p0.001_a_0 <- data_SCZ_New_100W
data_SCZ_New_100W_p0.001_a_0$beta<-y

save(data_SCZ_New_100W_p0.001_a_0,file="data_SCZ_New_100W_h0.5_p_0.001_a_0.Rdata")

####################
#compute the p and se

set.seed(1234)
h=0.5
p=0.01
num = 1000000
sd<-sqrt(h^2/(num*p))
y1 <- rnorm(num,mean = 0, sd = sd )
##plot(density(y1))
#sd<-sqrt(h^2/(num*p))
#1 a = 2, p = 0.01
load("data_SCZ_New_100W_p0.001_a0.72.RData")
load("data_SCZ_New_100W_p0.001_a_0.RData")
load("data_SCZ_New_100W_p0.01_a0.72.RData")
load("data_SCZ_New_100W_p0.01_a_0.RData")
load("data_SCZ_New_100W_p0.1_a0.72.RData")
load("data_SCZ_New_100W_p0.1_a_0.RData")

data_SCZ_New_100W_p0.001_a0.72$p<-unlist(lapply(rnorm(num,mean = 0, sd = 1 ),function(x) P_value(cdf=pnorm,x=x,paramet=c(0,1),side=0)))
data_SCZ_New_100W_p0.001_a_0$p<-unlist(lapply(rnorm(num,mean = 0, sd = 2 ),function(x) P_value(cdf=pnorm,x=x,paramet=c(0,2),side=0)))
data_SCZ_New_100W_p0.01_a0.72$p<-unlist(lapply(rnorm(num,mean = 0, sd = 3 ),function(x) P_value(cdf=pnorm,x=x,paramet=c(0,3),side=0)))
data_SCZ_New_100W_p0.01_a_0$p<-unlist(lapply(rnorm(num,mean = 0, sd = 4 ),function(x) P_value(cdf=pnorm,x=x,paramet=c(0,4),side=0)))
data_SCZ_New_100W_p0.1_a0.72$p<-unlist(lapply(rnorm(num,mean = 0, sd = 5 ),function(x) P_value(cdf=pnorm,x=x,paramet=c(0,5),side=0)))
data_SCZ_New_100W_p0.1_a_0$p<-unlist(lapply(rnorm(num,mean = 0, sd = 6 ),function(x) P_value(cdf=pnorm,x=x,paramet=c(0,6),side=0)))

data_SCZ_New_100W_p0.001_a0.72$se<-sqrt(((data_SCZ_New_100W_p0.001_a0.72$beta)^2)/qchisq(data_SCZ_New_100W_p0.001_a0.72$p,1,lower.tail=F))
data_SCZ_New_100W_p0.001_a_0$se<-sqrt(((data_SCZ_New_100W_p0.001_a_0$beta)^2)/qchisq(data_SCZ_New_100W_p0.001_a_0$p,1,lower.tail=F))
data_SCZ_New_100W_p0.01_a0.72$se<-sqrt(((data_SCZ_New_100W_p0.01_a0.72$beta)^2)/qchisq(data_SCZ_New_100W_p0.01_a0.72$p,1,lower.tail=F))
data_SCZ_New_100W_p0.01_a_0$se<-sqrt(((data_SCZ_New_100W_p0.01_a_0$beta)^2)/qchisq(data_SCZ_New_100W_p0.01_a_0$p,1,lower.tail=F))
data_SCZ_New_100W_p0.1_a0.72$se<-sqrt(((data_SCZ_New_100W_p0.1_a0.72$beta)^2)/qchisq(data_SCZ_New_100W_p0.1_a0.72$p,1,lower.tail=F))
data_SCZ_New_100W_p0.1_a_0$se<-sqrt(((data_SCZ_New_100W_p0.1_a_0$beta)^2)/qchisq(data_SCZ_New_100W_p0.1_a_0$p,1,lower.tail=F))

save(data_SCZ_New_100W_p0.001_a0.72,file="data_SCZ_New_100W_p0.001_a0.72.RData")
save(data_SCZ_New_100W_p0.001_a_0,file="data_SCZ_New_100W_p0.001_a_0.RData")
save(data_SCZ_New_100W_p0.01_a0.72,file="data_SCZ_New_100W_p0.01_a0.72.RData")
save(data_SCZ_New_100W_p0.01_a_0,file="data_SCZ_New_100W_p0.01_a_0.RData")
save(data_SCZ_New_100W_p0.1_a0.72,file="data_SCZ_New_100W_p0.1_a0.72.RData")
save(data_SCZ_New_100W_p0.1_a_0,file="data_SCZ_New_100W_p0.1_a_0.RData")

########################
#1000test for different data
setwd("/share/pub/dengcy/GWAS_Multiomics/test/1000test")
library(scPagwas)
library(ggplot2) 
library(ggthemes) 
library(ggpubr)
source("/share/pub/dengcy/GWAS_Multiomics/test/1000test/Pagwas_chunyu_permutation.R")

liver_data <- readRDS("/share/pub/qiuf/COVID/liver/Rpoly.rds")
liver_data@meta.data$annotation<-liver_data@meta.data$anno
Idents(liver_data)<-liver_data@meta.data$anno
liver_data <- NormalizeData(liver_data, normalization.method = "LogNormalize", scale.factor = 10000)
#select the 5000 cells randomly
set.seed(123)
liver_data <- liver_data[sample(1:ncol(liver_data),5000),]
Pagwas <- list()
Pagwas <- Single_data_input(Pagwas=Pagwas,
                            assay="RNA",
                            Single_data=liver_data,
                            Pathway_list=Genes_by_pathway_kegg)
Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,
                               Pathway_list=Genes_by_pathway_kegg)

for(i in c("p0.001_a0.72","p0.001_a_0","p0.01_a0.72","p0.01_a_0","p0.1_a0.72","p0.1_a_0")){
  load(paste0("data_SCZ_New_100W_",i,".RData"))
  assign(paste0("data_SCZ_New_100W_",i),get(paste0("data_SCZ_New_100W_",i)))
  perm_results<-lapply(1000:1999,function(seed){
    #set.seed(seed)
    gwas_test<-Pagwas_chunyu_permutation(data_su=get(paste0("data_SCZ_New_100W_",i)),seed=seed)
    
    Pagwas <- scPagwas_main (Pagwas = Pagwas,
                             gwas_data = gwas_test,
                             celltype=T,
                             single_cell=F,
                             Single_data = liver_data)
    bootstrap_results<- Pagwas$bootstrap_results
    bootstrap_results<-bootstrap_results[-1,]
    bootstrap_results$bp_value<- p.adjust(bootstrap_results$bp_value)
    return(bootstrap_results)
  })
  save(perm_results,file=paste0("liver_perm_results_",i,".RData"))
}

###########################
#pathway random test for different gene features
#############################

setwd("/share/pub/dengcy/GWAS_Multiomics/test/1000test")
library(scPagwas)
library(ggplot2)
library(ggthemes)
#compute the gene features 2000,4000,6000,8000,10000
for (gwas in c(c("p0.001_a0.72","p0.001_a_0","p0.01_a0.72","p0.01_a_0","p0.1_a0.72","p0.1_a_0")){
  load(paste0("data_SCZ_New_100W_",i,".RData"))
  assign(paste0("data_SCZ_New_100W_",i),get(paste0("data_SCZ_New_100W_",i)))
for(i in c(2000,4000,6000,8000,10000)){
  set.seed(123)
  liver_test<- Seurat::FindVariableFeatures(liver_test,
        selection.method = "vst",
        nfeatures = nfeatures
      )
  VariableFeatures <- Seurat::VariableFeatures(liver_test)
  liver_test2<-liver_test[,VariableFeatures]
  #select 50,100,150,200,250,300 pathways randomly
  for(j in c(50,100,150,200,250,300)){
    #20 times for each run
    perm_results<-lapply(sample(3000,20),function(m){
    set.seed(m)
    Pathway_list<-sample(Genes_by_pathway_kegg,j)
    Pagwas <- scPagwas_main(Pagwas=NULL,
                            gwas_data=get(paste0("data_SCZ_New_100W_",i)),
                            assay="RNA",
                            Single_data=liver_test2,                          celltype=T,
                             single_cell=F,
                            Pathway_list=Pathway_list)
    bootstrap_results<- Pagwas$bootstrap_results
    bootstrap_results<-bootstrap_results[-1,]
    bootstrap_results$bp_value<- p.adjust(bootstrap_results$bp_value)
    return(bootstrap_results)
    })
    save(perm_results,file=paste0("liver_perm_results_",i,"_",j,".RData"))
  }
}

##########raw result
for(i in c("p0.001_a0.72","p0.001_a_0","p0.01_a0.72","p0.01_a_0","p0.1_a0.72","p0.1_a_0")){
  load(paste0("data_SCZ_New_100W_",i,".RData"))
  assign(paste0("data_SCZ_New_100W_",i),get(paste0("data_SCZ_New_100W_",i)))
   Pagwas <- scPagwas_main(Pagwas=NULL,
                            gwas_data=get(paste0("data_SCZ_New_100W_",i)),
                            assay="RNA",
                            Single_data=liver_test,                          
                            celltype=T,
                             single_cell=F,
                            Pathway_list=Genes_by_pathway_kegg)
                            Pagwas$bootstrap_results$bp_value<- p.adjust(Pagwas$bootstrap_results$bp_value)
                            bootstrap_results<- Pagwas$bootstrap_results
save(bootstrap_results,file=paste0("liver_raw_results_",i,".RData"))
}