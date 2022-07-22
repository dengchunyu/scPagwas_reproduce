# Multiple disease traits for scPagwas banchmark

immune-/metabolism-related diseases/traits

brain-related diseases

c("MDD","ADHD","ASD","PD","BID","EP","Focal_EP","Generalized_EP","Juvenile_EP","Migraine","Multiple_sclerosis","smokingCessation","OpennessToExperence","Happiness","CognitivePerformance","SubjectWellBeing","BMI","WHR","WC","Obesity","TC","LDL","HDL","IBD","PBC","CoeliacDisease","T1D","DBP_EastAsian","DBP_European","PR_EastAsian","SBP","CAD","T2D","UlcerativeColitis","AD")

## 1. The result for four methods

#### 1.1 scPagwas ：

/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/scPgwas/

bootstrap_results_df1.RData

bootstrap_results_df2.RData

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/scPgwas/")
ds<-c("MDD","ADHD","ASD","PD","BID","EP","Focal_EP","Generalized_EP","Juvenile_EP","Migraine","Multiple_sclerosis","smokingCessation","OpennessToExperence","Happiness","CognitivePerformance","SubjectWellBeing","BMI","WHR","WC","Obesity","TC","LDL","HDL","IBD","PBC","CoeliacDisease","T1D","DBP_EastAsian","DBP_European","PR_EastAsian","SBP","CAD","T2D","UlcerativeColitis","AD")
scPagwas_results_list<-lapply(ds,function(i){
    message(i)
    load(paste0("scpagwas_brain_",i,".RData"))
    return(Pagwas$bootstrap_results)
})
names(scPagwas_results_list)<-ds
save(scPagwas_results_list,file="/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/scPagwas_results_list.RData")

scPagwas_HCA_list<-lapply(ds,function(i){
    message(i)
    load(paste0("scpagwas_HCA_",i,".RData"))
    return(Pagwas$bootstrap_results)
})
names(scPagwas_HCA_list)<-ds
save(scPagwas_HCA_list,file="/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/scPagwas_HCA_list.RData")

```

#### 1.2.magma 

地址：/share/pub/dengcy/Singlecell/COVID19/MAGMA

PR_EastAsian_magma_mouse_brain.gsa.out

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/magma/outresult")
ds<-c("MDD","ADHD","ASD","PD","BID","EP","Focal_EP","Generalized_EP","Juvenile_EP","Migraine","Multiple_sclerosis","smokingCessation","OpennessToExperence","Happiness","CognitivePerformance","SubjectWellBeing","BMI","WHR","WC","Obesity","TC","LDL","HDL","IBD","PBC","CoeliacDisease","T1D","DBP_EastAsian","DBP_European","PR_EastAsian","SBP","CAD","T2D","UlcerativeColitis","AD")
magma_result_list<-lapply(ds,function(i){
    message(i)
   result<-read.table(file=paste0(i,"_magma_mouse_brain.gsa.out"),header = T)
    return(result)
  #cell_mild$VARIABLE<-annotation_cell$V2  
})
names(magma_result_list)<-ds
save(magma_result_list,file="/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/magma_result_list.RData")

magma_HCA_list<-lapply(ds,function(i){
    message(i)
   result<-read.table(file=paste0(i,"_magma_HCA.gsa.out"),header = T)
    return(result)
  #cell_mild$VARIABLE<-annotation_cell$V2  
})
names(magma_HCA_list)<-ds
save(magma_HCA_list,file="/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/magma_HCA_list.RData")
```

#### 1.3.rolypoly

/share2/pub/jiangdp/jiangdp/COVID/rolypoly/

```R
setwd("/share2/pub/jiangdp/jiangdp/COVID/rolypoly_mouse")
ds<-c("MDD","ADHD","ASD","PD","BID","EP","Focal_EP","Generalized_EP","Juvenile_EP","Migraine","Multiple_sclerosis","smokingCessation","OpennessToExperence","Happiness","CognitivePerformance","SubjectWellBeing","BMI","WHR","WC","Obesity","TC","LDL","HDL","IBD","PBC","CoeliacDisease","T1D","DBP_EastAsian","DBP_European","PR_EastAsian","SBP","CAD","T2D","UlcerativeColitis","AD")

rolypoly_brain_list<-lapply(ds,function(i){
    message(i)
    load(paste0("roly_",i,".RData"))
    return(rolypoly_result$bootstrap_results)
})
names(rolypoly_brain_list)<-ds
save(rolypoly_brain_list,file="/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/rolypoly_brain_list.RData")

setwd("/share2/pub/jiangdp/jiangdp/COVID/rolypoly")
ds<-c("MDD","ADHD","ASD","PD","BID","EP","Focal_EP","Generalized_EP","Juvenile_EP","Migraine","Multiple_sclerosis","smokingCessation","OpennessToExperence","Happiness","CognitivePerformance","SubjectWellBeing","BMI","WHR","WC","Obesity","TC","LDL","HDL","IBD","PBC","CoeliacDisease","T1D","DBP_EastAsian","DBP_European","PR_EastAsian","SBP","CAD","T2D","UlcerativeColitis","AD")
rolypoly_HCA_list<-lapply(ds,function(i){
    message(i)
    load(paste0("roly_",i,".RData"))
    return(rolypoly_result$bootstrap_results)
})
names(rolypoly_HCA_list)<-ds
save(rolypoly_HCA_list,file="/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/rolypoly_HCA_list.RData")
```

#### 1.4.ldsc

/share2/pub/zhenggw/zhenggw/anaconda3/envs/ldsc/

```R
setwd("/share2/pub/zhenggw/zhenggw/scPagwas_LDSC/brain_results")
ds<-c("MDD","ADHD","ASD","PD","BID","EP","Focal_EP","Generalized_EP","Juvenile_EP","Migraine","Multiple_sclerosis","smokingCessation","OpennessToExperence","Happiness","CognitivePerformance","SubjectWellBeing","BMI","WHR","WC","Obesity","TC","LDL","HDL","IBD","PBC","CoeliacDisease","DBP_EastAsian","DBP_European","PR_EastAsian","SBP","CAD","T2D","UlcerativeColitis","AD")

ldsc_brain_result<-lapply(ds,function(i){
    message(i)
   result<-read.table(file=paste0(i,".cell_type_results.txt"),header = T)
    return(result)
})
names(ldsc_brain_result)<-ds
save(ldsc_brain_result,file="/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/ldsc_brain_result.RData")

ldsc_brain_p<-as.data.frame(lapply(ldsc_brain_result,function(df){
    rownames(df)<-df$Name
    df<-df[ldsc_brain_result[[1]]$Name,]
  a<- -log2(df$Coefficient_P_value)
  #a[which(a<4.32)]<-0
  return(a)
}))
rownames(ldsc_brain_p)<-ldsc_brain_result[[1]]$Name
write.csv(ldsc_brain_p,file="/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/ldsc_brain_p.csv")

################
setwd("/share2/pub/zhenggw/zhenggw/scPagwas_LDSC/HCA_results")
ds<-c("MDD","ADHD","ASD","PD","BID","EP","Focal_EP","Generalized_EP","Juvenile_EP","Migraine","Multiple_sclerosis","smokingCessation","OpennessToExperence","Happiness","CognitivePerformance","SubjectWellBeing","BMI","WHR","WC","Obesity","TC","LDL","HDL","IBD","PBC","CoeliacDisease","DBP_EastAsian","DBP_European","PR_EastAsian","SBP","CAD","T2D","UlcerativeColitis","AD")
ldsc_hca_result<-lapply(ds,function(i){
    message(i)
   result<-read.table(file=paste0(i,".cell_type_results.txt"),header = T)
    return(result)
  #cell_mild$VARIABLE<-annotation_cell$V2  
})
names(ldsc_hca_result)<-ds
save(ldsc_hca_result,file="/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/ldsc_hca_result.RData")

ldsc_hca_p<-as.data.frame(lapply(ldsc_hca_result,function(df){
    rownames(df)<-df$Name
    df<-df[ldsc_hca_result[[1]]$Name,]
  a<- -log2(df$Coefficient_P_value)
   
  #a[which(a<4.32)]<-0
  return(a)
}))
rownames(ldsc_hca_p)<-ldsc_hca_result[[1]]$Name
write.csv(ldsc_hca_p,file="/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/ldsc_hca_p.csv")

```

## 2. Visualization for integrate results

```R
setwd("E:/OneDrive/GWAS_Multiomics/banchmarkresult")
load("scPagwas_brain_list.RData")
load("magma_brain_list.RData")
load("rolypoly_brain_list.RData")
load("ldsc_brain_result.RData")

scPagwas_brain_p<-as.data.frame(lapply(scPagwas_results_list,function(df){
  a<- -log2(df$bp_value)
  a[which(a<4.32)]<-0
  return(a)
}))
#
rolypoly_brain_p<-as.data.frame(lapply(rolypoly_brain_list,function(df){
  a<- -log2(df$bp_value)
  a[which(a<4.32)]<-0
  return(a)
}))

rownames(scPagwas_brain_p)<-rownames(scPagwas_results_list[[1]])
rownames(rolypoly_brain_p)<-rownames(rolypoly_brain_list[[1]])
scPagwas_brain_p<-scPagwas_brain_p[-1,]
rolypoly_brain_p<-rolypoly_brain_p[-1,]
scPagwas_brain_p<-as.data.frame(t(scPagwas_brain_p))
rolypoly_brain_p<-as.data.frame(t(rolypoly_brain_p))

magma_brain_p<-as.data.frame(lapply(magma_result_list,function(df){
  a<- -log2(df$P)
  a[which(a<4.32)]<-0
  return(a)
}))
rownames(magma_brain_p)<-rownames(scPagwas_results_list[[1]])[-1]
magma_brain_p<-as.data.frame(t(magma_brain_p))

ldsc_brain_p<-as.data.frame(lapply(ldsc_brain_result,function(df){
  a<- -log2(df$Coefficient_P_value)
  a[which(a<4.32)]<-0
  return(a)
}))

rownames(ldsc_brain_p)<-rownames(scPagwas_results_list[[1]])[-1]
ldsc_brain_p<-as.data.frame(t(ldsc_brain_p))

rolypoly_brain_p<-rolypoly_brain_p[rownames(ldsc_brain_p),]
scPagwas_brain_p<-scPagwas_brain_p[rownames(ldsc_brain_p),]
magma_brain_p<-magma_brain_p[rownames(ldsc_brain_p),]

###########################################
types<-c("Neuropsychiatric disorders",
         "Neuropsychiatric disorders",
         "Neuropsychiatric disorders",
         "Neurodegenerative disorders",
         "Neuropsychiatric disorders",
         "Neurodegenerative disorders",
         "Neurodegenerative disorders",
         "Neurodegenerative disorders",
         "Neurodegenerative disorders",
         "Neurodegenerative disorders",
         "Neurodegenerative disorders",
         "Behavior-cognitive phenotype",
         "Behavior-cognitive phenotype",
         "Behavior-cognitive phenotype",
         "Behavior-cognitive phenotype",
         "Behavior-cognitive phenotype",
         "Metabolic disease",
         "Metabolic disease",
         "Metabolic disease",
         "Metabolic disease",
         "Metabolic disease",
         "Metabolic disease",
         "Metabolic disease",
         "Autoimmune disease",
         "Autoimmune disease",
         "Autoimmune disease",
         "Cardiovascular disease",
         "Cardiovascular disease",
         "Cardiovascular disease",
         "Cardiovascular disease",
         "Cardiovascular disease",
         "Cardiovascular disease",
         "Autoimmune disease",
         "Neurodegenerative disorders")
scPagwas_brain_p$method<-"scPagwas"
scPagwas_brain_p$phenotypes<-rownames(scPagwas_brain_p)
scPagwas_brain_p$types<-types


magma_brain_p$method<-"Magma"
magma_brain_p$phenotypes<-rownames(scPagwas_brain_p)
magma_brain_p$types<-types

rolypoly_brain_p$method<-"Rolypoly"
rolypoly_brain_p$phenotypes<-rownames(scPagwas_brain_p)
rolypoly_brain_p$types<-types

ldsc_brain_p$method<-"LDSC"
ldsc_brain_p$phenotypes<-rownames(scPagwas_brain_p)
ldsc_brain_p$types<-types


colnames(magma_brain_p)<-colnames(scPagwas_brain_p)
colnames(rolypoly_brain_p)<-colnames(scPagwas_brain_p)
colnames(ldsc_brain_p)<-colnames(scPagwas_brain_p)

rdf<-rbind(scPagwas_brain_p,
           magma_brain_p,
           rolypoly_brain_p,
           ldsc_brain_p)
#
rdf$phenotypes[which(rdf$phenotypes=="Multiple_sclerosis")]<-"MS"
rdf$phenotypes[which(rdf$phenotypes=="smokingCessation")]<-"SC"
rdf$phenotypes[which(rdf$phenotypes=="OpennessToExperence")]<-"Openness"
rdf$phenotypes[which(rdf$phenotypes=="CognitivePerformance")]<-"CP"
rdf$phenotypes[which(rdf$phenotypes=="SubjectWellBeing")]<-"SWB"
rdf$phenotypes[which(rdf$phenotypes=="UlcerativeColitis")]<-"UC"
rdf$phenotypes[which(rdf$phenotypes=="CoeliacDisease")]<-"CD"

library('ComplexHeatmap')
library(circlize)
col_fun = colorRamp2(c(0, 20), c("#05595B", "#FFD32D"))
col_fun(seq(0, 20))

color2<-colorRamp2(c(1, length(unique(rdf$phenotypes))), c("#085E7D", "#FFD32D"))(seq(1,length(unique(rdf$phenotypes))))
names(color2)<-unique(rdf$phenotypes)

ha<-rowAnnotation(df=data.frame(method=rdf$method,
                                phenotypes=rdf$phenotypes,
                                type=rdf$types),
                  col=list(#phenotypes=as.vector(color2),
                           type=c(`Autoimmune disease`="#3E8E7E",
                                  `Behavior-cognitive phenotype`="#FABB51",
                                  `Cardiovascular disease`="#FAEDC6",
                                  `Metabolic disease`="#97BFB4",
                                  `Neurodegenerative disorders`="#FFBC97",
                                  `Neuropsychiatric disorders`="#009DAE"),
                           method=c(scPagwas="#FF5959",
                                    Magma="#676FA3",
                                    Rolypoly="#CDDEFF",
                                    LDSC="#E6BA95"
                           ))
                  )
unique(rdf$phenotypes)
pt<-unique(rdf[,c("phenotypes","types")])
pt<-pt[order(pt$types),]
rdf$phenotypes <- factor(rdf$phenotypes,levels = pt$phenotypes)

pdf("Figure_banchmark.pdf",height =20,width =13)

Heatmap(as.matrix(rdf[,1:39]),
              col = col_fun,
              left_annotation = ha, 
              cluster_columns = T,
              cluster_rows = FALSE,
              color_space="HLS",
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "-log2(p)",
              row_order =order(rdf$types),
              show_row_names=FALSE,
              show_column_names=T,
              
              row_split=rdf$phenotypes)

dev.off()


tapply(1:nrow(rdf),factor(rdf$types),function(x){
  a<-unique(rdf$types[x])
  ad<-rdf[x,]
  ha<-rowAnnotation(df=data.frame(method=ad$method
                                  #phenotypes=ad$phenotypes
                                  ),
                    col=list(#phenotypes=as.vector(color2),
                      method=c(scPagwas="#D29D2B",
                               Magma="#FFF1CE",
                               Rolypoly="#86C6F4",
                               LDSC="#139487"
                      )
                    )
  )
  pdf(paste0("Banchmark_brain_",a,".pdf"),height =length(rdf$types[x])/2,width =13)
  
  print(Heatmap(as.matrix(ad[,1:39]),
          col = col_fun,
          left_annotation = ha, 
          cluster_columns = T,
          cluster_rows = FALSE,
          color_space="HLS",
          border=T,row_gap = unit(0.25, "mm"),
          show_parent_dend_line=T,
          name = "-log2(p)",
          show_row_names=FALSE,
          show_column_names=T,
          row_split=ad$phenotypes))
  dev.off()
})

```


