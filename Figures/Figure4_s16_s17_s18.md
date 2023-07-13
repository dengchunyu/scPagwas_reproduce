# BMMC for Blood traits Integrate visualization

## visualization for BMMC

Figure 4A

### 1.load the result and plot

```R
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(stringr)
library(ggtext)
library(scPagwas)
suppressMessages(library(Seurat))
library(org.Hs.eg.db)
library(dplyr)

setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/IntegrateVisulize")

traits<-c("eosinophilcount","basophilcount","LymphocytePercent","Lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","MeanCorpuscularHemoglobin","MeanCorpusVolume")

Seu_Hema_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")

pdf(file="Figure4a.pdf",height=10,width=10)
DimPlot(Seu_Hema_data,group.by="celltypes",reduction="tsne",pt.size=0.5,label = TRUE, repel=TRUE,label.size = 4)+ 
umap_theme()+ labs(x="TSNE",y="")+
ggtitle("BMMC")+
        scale_colour_manual(name = "celltypes", values = color_seruat_viridis28) +
        theme(aspect.ratio=1)
dev.off()

```
### 2.Pheatmap: BMMC mean score for celltypes

Figure4B

```R
library(gplots)
library(RColorBrewer)
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
library('ComplexHeatmap')
library(circlize)
traits<-c("basophilcount","eosinophilcount" ,"Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

result_TRS<-lapply(traits,function(i){
   print(i)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas.RData"))
    a<-data.frame(celltypes=Pagwas$celltypes,trs=Pagwas$scPagwas.topgenes.Score1)
     a$celltypes<-factor(a$celltypes,levels=unique(a$celltypes))
    #a have tow column, celltypes and trs, celltypes are characters, compute the mean trs for each celltype in a dataframe
    df<-aggregate(trs ~ celltypes, a, mean)
    return(df$trs)
})
result_TRS<-as.data.frame(result_TRS)
colnames(result_TRS)<-traits
rownames(result_TRS)<- unique(Pagwas$celltypes)
save(result_TRS,file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/result_TRS.RData")
load("D:/OneDrive/GWAS_Multiomics/Compare/Hema_test2/result_TRS.RData")
load("D:/OneDrive/GWAS_Multiomics/Compare/scPagwas_bmmc_list.RData")
result_list<-result_list[rownames(result_TRS),]
result_list_adj<- apply(result_list,2, function(x) p.adjust(x,method = "BH"))
list<-lapply(1:10,function(i){
  a<-data.frame(result_TRS[,i],result_list_adj[,i])
  colnames(a)<-paste0(colnames(result_TRS)[i],c("_mean_TRS","_FDR"))
  return(a)
})
result<-do.call(cbind,list)
write.csv(result,file="D:/OneDrive/GWAS_Multiomics/Manuscripts/Revise_comments/CelltypeP/Figure4b_bloodtrait_celltype_p.csv")

#coul <- colorRampPalette(brewer.pal(8, "Oranges"))(25)
coul <- colorRampPalette(c("#F5E8C7", "#AC7088"))(6)
p2star <- function(p){
  symnum(p,cutpoints = c(0,0.001,0.01,0.05,1),
         symbols = c('***','**','*',''),na = NA)
}
hM <- apply(result_list_adj,2, function(x) as.character(p2star(x)))
rownames(hM)<-rownames(result_list)

hM[which(result_TRS < 0.2)]<-""

result_TRS<-as.matrix(result_TRS)
result_TRS<-result_TRS[!(rownames(result_TRS) %in% c("14_Unk","26_Unk")),]
result_TRS<-result_TRS[c("10_cDC","09_pDC","16_Pre.B","18_Plasma","17_B","13_CD16.Mono","12_CD14.Mono.2","11_CD14.Mono.1","25_NK","24_CD8.CM","22_CD4.M","23_CD8.EM","19_CD8.N","20_CD4.N1","21_CD4.N2","06_CLP.1","15_CLP.2","05_CMP.LMPP","07_GMP","01_HSC","04_Early.Baso",
"02_Early.Eryth","03_Late.Eryth","08_GMP.Neut"),c("Lymphocytecount3","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume","WhiteBloodCellcount","neutrophilcount","eosinophilcount","basophilcount","monocytecount")]
hM<-hM[c("10_cDC","09_pDC","16_Pre.B","18_Plasma","17_B","13_CD16.Mono","12_CD14.Mono.2","11_CD14.Mono.1","25_NK","24_CD8.CM","22_CD4.M","23_CD8.EM","19_CD8.N","20_CD4.N1","21_CD4.N2","06_CLP.1","15_CLP.2","05_CMP.LMPP","07_GMP","01_HSC","04_Early.Baso",
"02_Early.Eryth","03_Late.Eryth","08_GMP.Neut"),c("Lymphocytecount3","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume","WhiteBloodCellcount","neutrophilcount","eosinophilcount","basophilcount","monocytecount")]
colnames(result_TRS)<-c("Lymphocytecount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume","WhiteBloodCellcount","Neutrophilcount","Eosinophilcount","Basophilcount","Monocytecount")
pdf("D:/OneDrive/GWAS_Multiomics/Manuscripts/Revise_comments/CelltypeP/Figure_scPagwas_bmmc_cellytpe_p.pdf",
    height =8,width =6)
heatmap.2(result_TRS,
          trace="none",#
          col=coul,#
          density.info = "none",
            breaks = seq(0, 0.6, length.out = 7),
          key.xlab ='Correlation',
          key.title = "",
          cexRow = 1,cexCol = 1,#
          Rowv = F,Colv =F, #
          margins = c(10,10),
          cellnote = hM,notecol='black'
)
dev.off()
```

### 2.tsne plot for different traits
the tsne plot for Figure4CDE
Supplementary Figure S16

```R

for(i in c("Lymphocytecount3","monocytecount","MeanCorpusVolume")){
 load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas.RData"))

    all_fortify_can <- fortify.Seurat.tsne(Pagwas)
        p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = TSNE_1, y = TSNE_2,color =scPagwas.topgenes.Score1), size = 0.2, alpha = 1) +
        umap_theme() +
 scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle(i)

        pdf(file=paste0("TSNE.",i,".scPagwas.TRS.Score1.pdf"),width = 8, height = 8)
        print(p1)
        dev.off()

        p2<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = TSNE_1, y = TSNE_2,color =magma1000_scdrs.zscore), size = 0.2, alpha = 1) +
        umap_theme() +
 scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle(i)

        pdf(file=paste0("TSNE.",i,".magma1000_scdrs.zscore.pdf"),width = 8, height = 8)
        print(p2)
        dev.off()
}

```

### 3.Barplot and dotplot for scPagwas TRS score

the boxplot for Figure4CDE

```R

######################
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)
setwd("E:/OneDrive/GWAS_Multiomics/Compare/IntegrateVisulize")
######################
traits2<-c("Lymphocytecount3","MeanCorpuscularHemoglobin","MeanCorpusVolume","monocytecount","Plateletcount")

celltyperank<-c("01_HSC", "05_CMP.LMPP","06_CLP.1","15_CLP.2","07_GMP",
                "02_Early.Eryth" , "03_Late.Eryth" , "04_Early.Baso","08_GMP.Neut", 
                "09_pDC","10_cDC" ,
                "11_CD14.Mono.1","12_CD14.Mono.2","13_CD16.Mono","14_Unk",
                "16_Pre.B","17_B","18_Plasma",
                "19_CD8.N","20_CD4.N1","21_CD4.N2","22_CD4.M","23_CD8.EM","24_CD8.CM",
                "25_NK",
                "26_Unk"
)
celltypecolor<-c("#CB181D","#D7403E","#E3695F","#EF9280",
                 "#CD9B1D","#EEB422","#F1C148","#FFE557","#FFF3B0",
                 "#2171B5","#73A5D2",
                 "#7F2704","#A95426","#BE6A37","#D38148",
                 "#3F007D","#653698","#8C6DB3",
                 "#0E4F27","#2A653F","#39704C","#558664","#80A789","#BAD3BB",
                 "#763857",
                 "#555555"
)

for(i in traits2){
  load(paste0(i,"_meta.data.RData")) 
  head(meta.data)
  meta.data$celltypes<-factor(meta.data$celltypes,levels =celltyperank )
  p <- ggplot(meta.data, aes(x = celltypes, 
                             y =scPagwas1000_scdrs.raw_score,
                             fill=celltypes,
                             color=celltypes))+
    geom_boxplot(outlier.size = 0.2,alpha=1)+
    scale_fill_manual(values = celltypecolor)+
    scale_color_manual(values = celltypecolor)+
    theme_light()+
    theme(axis.text.x = element_text(angle =90,vjust = 1,hjust = 1),
          legend.position="none")+
    labs(x = "",y="",title=i)
  pdf(paste0(i,"boxplot_scPagwas.pdf"),width =8,height = 5)
  print(p)
  dev.off()
 }
```

## visualization for pbmc

### 1. PBMC Lymphocytecount

**Supplementary Figure S17**
Figure S17b
```R
library(scPagwas)
suppressMessages(library(Seurat))
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)
library(ggtext)
#Figure S17b
load("/share/pub/dengcy/GWAS_Multiomics/compare/Lymphocytecount3_pbmc_scPagwas_singlecell.RData")
all_fortify_can <- fortify.Seurat.umap(Pagwas)
    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas.topgenes.Score1), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D",
                             midpoint = median(all_fortify_can$scPagwas.topgenes.Score1))+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("Lymphocytecount pbmc")

        pdf(file="/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/pbmc.Umap.Lymphocytecount3.scPagwas1000.pdf",width = 8, height = 8)
        print(p1)
        dev.off()
```

### 2. PBMC Monocyte
Figure S17a
Supplementary Figure S18.
```R
library(scPagwas)
 suppressMessages(library(Seurat))
 library(ggplot2)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)
library(ggtext)

load("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/monocytecount_pbmc_scPagwas_singlecell.RData")

all_fortify_can <- fortify.Seurat.umap(Pagwas)
    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas.topgenes.Score1), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D",
                             midpoint = median(all_fortify_can$scPagwas.topgenes.Score1))+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("Monocytecount pbmc")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/pbmc.Umap.monocytecount.scPagwas1000.pdf",width = 8, height = 8)
print(p1)
dev.off()

#Supplementary Figure S18
all_fortify_can <- fortify.Seurat.umap(Pagwas)
    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas1000_scdrs.zscore), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D",
                             midpoint = median(all_fortify_can$scPagwas1000_scdrs.zscore))+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("scPagwas + scDRS")

    p2<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =magma1000_scdrs.zscore), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D",
                             midpoint = median(all_fortify_can$magma1000_scdrs.zscore))+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("magma + scDRS")
        
load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas.RData")
all_fortify_can <- fortify.Seurat.umap(Pagwas)
    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas1000_scdrs.zscore), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D",
                             midpoint = median(all_fortify_can$scPagwas1000_scdrs.zscore))+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("scPagwas + scDRS")

    p2<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =magma1000_scdrs.zscore), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D",
                             midpoint = median(all_fortify_can$magma1000_scdrs.zscore))+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("magma + scDRS")

```

## sub-function

```R
umap_theme <- function() {
  theme_grey() %+replace%
    theme(
      panel.background = element_rect(fill = "white", colour = "black", size = 2),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_blank()
    )
}


fortify.Seurat.umap <- function(x) {
  xy1 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "umap"))
  colnames(xy1) <- c("UMAP_1", "UMAP_2")
  xy1$UMAP_1 <- as.numeric(xy1$UMAP_1)
  xy1$UMAP_2 <- as.numeric(xy1$UMAP_2)

  return(cbind(xy1, as.data.frame(x@meta.data)))
}


#' fortify.Seurat.tsne
#' @description set data frame to ggplot
#' @param x seruat
#' @export
#' @return

fortify.Seurat.tsne <- function(x) {
  xy2 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "tsne"))
  colnames(xy2) <- c("TSNE_1", "TSNE_2")
  xy2$TSNE_1 <- as.numeric(xy2$TSNE_1)
  xy2$TSNE_2 <- as.numeric(xy2$TSNE_2)

  return(cbind(xy2, as.data.frame(x@meta.data)))
}

color_seruat_patient <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")

color_seruat_viridis28 <- c("#D9DD6B","#ECEFA4","#D54C4C","#8D2828","#FDD2BF","#E98580","#DF5E5E","#492F10","#334257","#476072","#548CA8",
"#00A19D","#ECD662","#5D8233","#284E78","#3E215D","#835151","#F08FC0","#C6B4CE","#BB8760","#FFDADA","#3C5186",
"#558776","#E99497","#FFBD9B","#0A1D37","#01937C","#464660","#368B85")
color_seruat_13 <- c("#FDD2BF","#DF5E5E","#492F10","#753422","#628395","#262A53","#ECD662","#5D8233","#284E78","#3E215D","#DF711B","#FFB740","#64C9CF")

color_seruat_type3 <-c("#D9DD6B","#D54C4C","#548CA8")
color_seruat_8 <- c("#FDD2BF","#753422","#297F87","#DF5E5E","#628395","#262A53","#ECD662","#5D8233")
```

