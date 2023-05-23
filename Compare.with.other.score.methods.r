library(Seurat)
library(SeuratObject)
library("dplyr")
library("foreach")
library("ggplot2")
library(reshape2)  # melt
library(gridExtra) # arrangeGrob
library(ggpubr)    # as_ggplot
library(cowplot)   # draw_plot_label
library(ggtext)
library(testSctpa)
library(VISION)
library(AUCell)


##输出mgt文件
gwass<-c('eosinophilcount', 'basophilcount','LymphocytePercent','monocytecount' ,'neutrophilcount', 'WhiteBloodCellcount', 'Hemoglobinconcen', 'MeanCorpuscularHemoglobin', 'MeanCorpusVolume')
gwass<-c('Lymphocytecount3')
methods<-c('smultixcan', 'spredixcan' ,'TWAS')

for (gwas in gwass) {

	for (method in methods) {
		files<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.csv"))

genes<-read.csv(files)
rn= glue::glue("{method}_top")
b<-paste0(rn,1000)

a<-c(b,"NA",genes$gene_name)
a<-matrix(a,nrow=1)
out<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
write.table(a,file=out,append = T,col.names = F,row.names = F,quote = F,sep = "\t")
	}
}

##计算得分
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")
load("Pagwas_monocytecount_nk_v1.9.1.RData")
counts = load_counts()
se_oj = CreateSeuratObject(counts)
se_oj = cal_PAS(seurat_object = se_oj,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = 'log',
              species = 'mouse', 
              pathway='kegg')

DefaultAssay(object = Pagwas) <- "RNA"

gwas='monocytecount'
m_l<-lapply(methods,function(method){
  gmt<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
  DefaultAssay(object = Pagwas) <- "RNA"

  auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)

  rn= glue::glue("{method}-top")
  b<-paste0(rn,1000)
  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
  return(auccell_df[,b])
	})
auc_df<-as.data.frame(m_l)
colnames(auc_df)<-paste0(methods,"_monocytecount_AUCell")
Pagwas@meta.data<-cbind(Pagwas@meta.data,auc_df)


m_l<-lapply(methods,function(method){
  gmt<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
  DefaultAssay(object = Pagwas) <- "RNA"

  auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'Vision',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)

  rn= glue::glue("{method}-top")
  b<-paste0(rn,1000)
  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
  return(auccell_df[,b])
	})
df<-as.data.frame(m_l)
colnames(df)<-paste0(methods,"_monocytecount_Vision")
Pagwas@meta.data<-cbind(Pagwas@meta.data,df)


gwas='Lymphocytecount3'
m_l<-lapply(methods,function(method){
  gmt<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
  DefaultAssay(object = Pagwas) <- "RNA"

  auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)

  rn= glue::glue("{method}-top")
  b<-paste0(rn,1000)
  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
  return(auccell_df[,b])
	})
auc_df<-as.data.frame(m_l)
colnames(auc_df)<-paste0(methods,"_Lymphocytecount3_AUCell")
Pagwas@meta.data<-cbind(Pagwas@meta.data,auc_df)


m_l<-lapply(methods,function(method){
  gmt<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
  DefaultAssay(object = Pagwas) <- "RNA"

  auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'Vision',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)

  rn= glue::glue("{method}-top")
  b<-paste0(rn,1000)
  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
  return(auccell_df[,b])
	})
df<-as.data.frame(m_l)
colnames(df)<-paste0(methods,"_Lymphocytecount3_Vision")
Pagwas@meta.data<-cbind(Pagwas@meta.data,df)

save(Pagwas,file="Pagwas_monocytecount_nk_v1.9.1.RData")


############################
#Lymphocytecount3
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Pagwas_Lymphocytecount_modeldata_v1.10.0.RData")
gwas='monocytecount'
m_l<-lapply(methods,function(method){
  gmt<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
  DefaultAssay(object = Pagwas) <- "RNA"

  auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)

  rn= glue::glue("{method}-top")
  b<-paste0(rn,1000)
  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
  return(auccell_df[,b])
	})
auc_df<-as.data.frame(m_l)
colnames(auc_df)<-paste0(methods,"_monocytecount_AUCell")
Pagwas@meta.data<-cbind(Pagwas@meta.data,auc_df)


m_l<-lapply(methods,function(method){
  gmt<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
  DefaultAssay(object = Pagwas) <- "RNA"

  auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'Vision',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)

  rn= glue::glue("{method}-top")
  b<-paste0(rn,1000)
  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
  return(auccell_df[,b])
	})
df<-as.data.frame(m_l)
colnames(df)<-paste0(methods,"_monocytecount_Vision")
Pagwas@meta.data<-cbind(Pagwas@meta.data,df)


gwas='Lymphocytecount3'
m_l<-lapply(methods,function(method){
  gmt<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
  DefaultAssay(object = Pagwas) <- "RNA"

  auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)

  rn= glue::glue("{method}-top")
  b<-paste0(rn,1000)
  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
  return(auccell_df[,b])
	})
auc_df<-as.data.frame(m_l)
colnames(auc_df)<-paste0(methods,"_Lymphocytecount3_AUCell")
Pagwas@meta.data<-cbind(Pagwas@meta.data,auc_df)


m_l<-lapply(methods,function(method){
  gmt<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
  DefaultAssay(object = Pagwas) <- "RNA"

  auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'Vision',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)

  rn= glue::glue("{method}-top")
  b<-paste0(rn,1000)
  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
  return(auccell_df[,b])
	})
df<-as.data.frame(m_l)
colnames(df)<-paste0(methods,"_Lymphocytecount3_Vision")
Pagwas@meta.data<-cbind(Pagwas@meta.data,df)

save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Pagwas_Lymphocytecount_modeldata_v1.10.0.RData")


library(circlize)
library(dplyr)
require(pROC)
require(ggplot2)
##不同数量的top基因的打分

load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_nk_v1.9.1.RData")
all_fortify_can <- fortify.Seurat.umap(Pagwas)

auc_list1<-lapply(all_fortify_can[,c(29:32,12,33:37)],function(x){
  roc(predictor=x,response=all_fortify_can$type)  
})

names(auc_list1)<-colnames(all_fortify_can)[c(29:32,12,33:37)]
 auc_l1<- unlist(lapply(1:length(auc_list1),function(x) round(as.numeric(auc_list1[[x]]["auc"]),3)))
 names(auc_list1)<- paste0(names(auc_list1),"(AUC=",auc_l1,")")
 pdf("AUC_modeldata.pagwas.100.1000.pdf",height = 6)
 ggroc(auc_list1, linetype = 2, size = 1,alpha=0.8)+
    ggtitle("ROC curve for different number of gene") + 
    theme_classic()+ggsci::scale_color_lancet()+
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
                 color="grey", linetype="dashed")
    #theme(panel.background=element_rect(fill="white",colour="blue"))
  dev.off()

