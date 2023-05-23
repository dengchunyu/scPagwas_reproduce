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

pdf(file="Seu_Hema_seurat_TSNE.pdf",height=10,width=10)
DimPlot(Seu_Hema_data,group.by="celltypes",reduction="tsne",pt.size=0.5,label = TRUE, repel=TRUE,label.size = 4)+ 
umap_theme()+ labs(x="TSNE",y="")+
ggtitle("BMMC")+
        scale_colour_manual(name = "celltypes", values = color_scanpy_viridis28) +
        theme(aspect.ratio=1)
dev.off()

```

### 2.tsne plot for different traits

Figure4CDE

```R

for(i in c("Lymphocytecount3","monocytecount","MeanCorpusVolume")){
 load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))

    all_fortify_can <- fortify.Seurat.tsne(Pagwas)
        p2<- ggplot() +
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
        print(p2)
        dev.off()
}

```

### 3.Barplot and dotplot for scPagwas TRS score

Figure4CDE

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
  
##########scPagwas500_scdrs.raw_score
 a1<- tapply(meta.data$scPagwas1000_scdrs.raw_score, meta.data$celltypes, function(x){
    mean(x)
  })
 a2<- tapply(meta.data$scPagwas500_scdrs.raw_score, meta.data$celltypes, function(x){
   mean(x)
 })
 a3<- tapply(meta.data$magma1000_scdrs.raw_score, meta.data$celltypes, function(x){
   mean(x)
 })
 a4<- tapply(meta.data$magma500_scdrs.raw_score, meta.data$celltypes, function(x){
   mean(x)
 })
 
  df<-data.frame(scPagwas1000=a1,
                 scPagwas500=a2,
                 magma1000=a3,
                 magma500=a4,
                 celltypes=names(a1))
  ggdf<-reshape2::melt( df,id.vars="celltypes")
  ggdf$celltypes<-factor(ggdf$celltypes,levels =celltyperank )
  ggdf$variable<-factor(ggdf$variable,levels =c("magma500","magma1000","scPagwas500","scPagwas1000") )
  p2 = ggplot(ggdf,
             aes(y=variable,
                 x = celltypes)) +
    geom_point(aes(colour=celltypes, size=value)) + 
    scale_color_manual(values = celltypecolor,guide = "none")+
    theme_light()#+
 pdf(paste0(i,"dotplot_score.pdf"),width =9,height =2.5)
 print(p2)
 dev.off()

##########cell p value percent plot
 meta.data$CellsrankPvalue<-p.adjust(meta.data$CellsrankPvalue,method = "bonferroni")
 a<- tapply(meta.data$CellsrankPvalue, meta.data$celltypes, function(x){
   sum(x<0.05)/length(x)
 })
 
 df<-data.frame(Cells.adjPvalue=a,celltypes=names(a))
 df$celltypes<-factor(df$celltypes,levels =celltyperank )
 df$scPagwas<-"Cells.adjPvalue"
 p3 = ggplot(df,
             aes(y=scPagwas,
                 x = celltypes)) +
   geom_point(aes(colour=celltypes, size=Cells.adjPvalue)) + 
   scale_color_manual(values = celltypecolor,guide = "none")+
   #scale_size(range=c(0, max.size)) +
   theme_light()#+
 #theme(legend.position="none")
 pdf(paste0(i,"dotplot_Cells.adjPvalue.pdf"),width =9,height = 1)
 print(p3)
 dev.off()
}
```

## visualization for pbmc

### 1. PBMC Lymphocytecount

**Supplementary Figure S7**

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 library(ggplot2)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)
library(ggtext)

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

### 2.PBMC Hemoglobinconcen

The names for TRS score is different for last version.

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 library(ggplot2)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)
library(ggtext)

load("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/Hemoglobinconcen_pbmc_scPagwas_singlecell.RData")

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
        ggtitle("Hemoglobinconcen pbmc")

        pdf(file="/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/pbmc.Umap.Hemoglobinconcen.scPagwas1000.pdf",width = 8, height = 8)
        print(p1)
        dev.off()
```

### 3.MeanCorpusVolume

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 library(ggplot2)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)
library(ggtext)

load("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/MeanCorpusVolume_pbmc_scPagwas_singlecell.RData")

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
        ggtitle("MeanCorpusVolume pbmc")

        pdf(file="/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/pbmc.Umap.MeanCorpusVolume.scPagwas1000.pdf",width = 8, height = 8)
print(p1)
dev.off()
```

4. ### Monocyte

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

color_scanpy_patient <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")

color_scanpy_viridis28 <- c("#D9DD6B","#ECEFA4","#D54C4C","#8D2828","#FDD2BF","#E98580","#DF5E5E","#492F10","#334257","#476072","#548CA8",
"#00A19D","#ECD662","#5D8233","#284E78","#3E215D","#835151","#F08FC0","#C6B4CE","#BB8760","#FFDADA","#3C5186",
"#558776","#E99497","#FFBD9B","#0A1D37","#01937C","#464660","#368B85")
color_scanpy_13 <- c("#FDD2BF","#DF5E5E","#492F10","#753422","#628395","#262A53","#ECD662","#5D8233","#284E78","#3E215D","#DF711B","#FFB740","#64C9CF")

color_scanpy_type3 <-c("#D9DD6B","#D54C4C","#548CA8")
color_scanpy_8 <- c("#FDD2BF","#753422","#297F87","#DF5E5E","#628395","#262A53","#ECD662","#5D8233")
```

