
#2022-11-28

#@' Permutation analysis for 1000 genes with correlation with RISmed_search across 1000 times of random selection

options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))

# Reference website
# https://amunategui.github.io/pubmed-query/

#gene<-read.table("E:/brain/01-data/healthy/QC/gene.txt",header=F)
gene<-read.table("gene.txt",header=F)
library(RISmed)
set.seed(123456)

temp<-gene
cor<-NULL
for(j in 1:10){
  temp[,2]<-gene[sample(1:10850,10850),2]
  traits<-temp[order(temp[,2],decreasing=T)[1:1000],]
  #row.names(traits)<-traits[,1]
  PMID<-NULL
  for(i in 1:length(traits[,1])){
    print(paste0("j=",j))
    print(paste0("i=",i))
    tryCatch({    
    #  Sys.sleep(1)
      res <- EUtilsSummary(paste0("(",traits[i,1],") AND (Alzheimer) "), type='esearch', db='pubmed')
      
      PMID<-rbind(PMID,c(traits[i,1],"AD",QueryCount(res)))
    },error=function(e){
      tryCatch({
        Sys.sleep(5)
        res <- EUtilsSummary(paste0("(",traits[i,1],") AND (Alzheimer) "), type='esearch', db='pubmed')
        PMID<-rbind(PMID,c(traits[i,1],"AD",QueryCount(res)))
      },error=function(e){
        Sys.sleep(20)
        res <- EUtilsSummary(paste0("(",traits[i,1],") AND (Alzheimer) "), type='esearch', db='pubmed')
        PMID<-rbind(PMID,c(traits[i,1],"AD",QueryCount(res)))
        
      })

      
    })
  }
  cor<-c(cor,cor(traits[na.omit(match(PMID[,1],traits[,1])),2],log2(as.numeric(PMID[,3])+1)))
}
write.table(cor,"cor_100.txt",quote = F,row.names = F)

#@ We need to do our best to get the points for improving all the softwares for human genetics



#########Permutation_plot for random selections

results_permut <- read.table("cor_all.txt",header = T)

hist(results_permut[,1], col="#D2691E",xlab="Counts of overlapped genes",xlim = c(-0.2,0.2),main=NULL)
abline(v=0.21,col="blue",lty="longdash")
P_value= length(results_permut$x[results_permut$x>0.21])/length(results_permut$x)















