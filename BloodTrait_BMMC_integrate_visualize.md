# BMMC for Blood traits Integrate visualization

## visualization for BMMC

### 1.load the result and 

```R
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(stringr)
library(ggtext)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/IntegrateVisulize")
 library(scPagwas)
 suppressMessages(library(Seurat))
library(org.Hs.eg.db)
library(dplyr)
traits<-c("eosinophilcount","basophilcount","LymphocytePercent","Lymphocytecount","monocytecount","neutrophilcount","WhiteBloodCellcount","MeanCorpuscularHemoglobin","MeanCorpusVolume")
readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")

pdf(file="Seu_Hema_seurat_TSNE.pdf",height=10,width=10)
DimPlot(Pagwas,group.by="celltypes",reduction="tsne",pt.size=0.5,label = TRUE, repel=TRUE,label.size = 4)+ 
umap_theme()+ labs(x="TSNE",y="")+
ggtitle("BMMC")+
        scale_colour_manual(name = "celltypes", values = color_scanpy_viridis28) +
        theme(aspect.ratio=1)
dev.off()

```

### 2.tsne plot for different traits

```R
i<-"Lymphocytecount2"
i<-"monocytecount"
i<-"MeanCorpusVolume"
i<-"Hemoglobinconcen"

load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scPagwas1.8/",i,"_Hema_bmmc_scPagwas_v1.9.RData"))

    all_fortify_can <- fortify.Seurat.tsne(Pagwas)
    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = TSNE_1, y = TSNE_2,color =sclm_score), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle(i)

        pdf(file=paste0("TSNE.",i,".sclm_score.pdf"),width = 8, height = 8)
        print(p1)
        dev.off()

        p2<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = TSNE_1, y = TSNE_2,color =scPagwas.lmtopgenes.Score1), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D",
                             midpoint = median(all_fortify_can$scPagwas.lmtopgenes.Score1)
                              )+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle(i)

        pdf(file=paste0("TSNE.",i,".scPagwas.topgenes.Score1.pdf"),width = 8, height = 8)
        print(p2)
        dev.off()

```

### 3.Barplot for scPagwas TRS score

```R
setwd("E:/OneDrive/GWAS_Multiomics/Compare/IntegrateVisulize")
######################
traits2<-c("LymphocytePercent","MeanCorpuscularHemoglobin","MeanCorpusVolume","monocytecount","Plateletcount")

 celltyperank<-c("01_HSC", "05_CMP.LMPP","06_CLP.1","15_CLP.2","07_GMP",
      "02_Early.Eryth" , "03_Late.Eryth" , "04_Early.Baso","08_GMP.Neut", 
      "09_pDC","10_cDC" ,
     "11_CD14.Mono.1","12_CD14.Mono.2","13_CD16.Mono","14_Unk",
     "16_Pre.B","17_B","18_Plasma",
      "19_CD8.N","20_CD4.N1","21_CD4.N2","22_CD4.M","23_CD8.EM","24_CD8.CM",
      "25_NK",
      "26_Unk"
     )
celltypecolor<-c("FD5D5D","#FF8080","#FFC3C3","#FFCBCB",
    "#E26A2C","#FF6701","#FF8243","#FDA65D","#FFD07F",
    "#219F94","#C1DEAE",
    "#52734D","#91C788","#DDFFBC","#FEFFDE",
    "#F875AA","#FBACCC","#F1D1D0",
"#2F8F9D","#68A7AD","#3BACB6","#82DBD8","#B3E8E5","#CCF3EE",
    "#7C9473",
    "#AA8976"
)
for(i in traits){
    print(i)
    i<-"monocytecount"
 load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_meta.data.RData")) 

p <- ggplot(gg_gsva2, aes(x = ImmuneCell, y =ssGsva_score,fill=cluster))+
geom_boxplot(outlier.size = 0.6,alpha=0.6)+theme_classic()+
theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),legend.position = "top")+labs(x = "",y="",title=i)
    
    }

library(ggplot2)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)
setwd("E:/OneDrive/GWAS_Multiomics/Compare/IntegrateVisulize")
######################
traits2<-c("LymphocytePercent","MeanCorpuscularHemoglobin","MeanCorpusVolume","monocytecount","Plateletcount")

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

i<-"monocytecount"
i<-"MeanCorpusVolume"
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

```

### 4 Lymphocyte：

```R
load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas_v1.7.3.RData") 
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(stringr)
library(ggtext)
setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/IntegrateVisulize")
 library(scPagwas)
 suppressMessages(library(Seurat))

#Pagwas@misca
scPagwastop1000 <- names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation, decreasing = T),])[1:1000]
scPagwastop500 <- names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation, decreasing = T),])[1:500]

load(paste0(i,"_magma_genes.RData"))

magmatop1000<-intersect(magma_genes$symbol[1:1000],rownames(Pagwas))
 magmatop500<-intersect(magma_genes$symbol[1:500],rownames(Pagwas))

topgene<-list(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500)
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")


    Pagwas <- AddModuleScore(Pagwas, assay = "RNA", topgene, name = c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500"))
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas_v1.7.3.RData")
    
    all_fortify_can <- fortify.Seurat.tsne(Pagwas)

save(all_fortify_can,file="all_fortify_can.RData")
    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = TSNE_1, y = TSNE_2,color =scPagwas.topgenes.Score1), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D",
                             midpoint = median(all_fortify_can$scPagwas.topgenes.Score1))+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("Lymphocytecount")
ggsave(p1,file="E:/OneDrive/GWAS_Multiomics/Compare/IntegrateVisulize/TSNE.Lymphocytecount3.scPagwas1000.pdf",width = 8, height = 8)
        pdf(file="TSNE.Lymphocytecount3.scPagwas1000.pdf",width = 8, height = 8)
        print(p1)
        dev.off()

```

#### Lymphocyte run magma

```R
i<-"Lymphocytecount3"
gwas<-bigreadr::fread2(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_gwas_data.txt"))
gwas$N<-62076
print(head(gwas))
magma_Input1<-gwas[,c("rsid","p","N")]
magma_Input2<-gwas[,c("rsid","chrom","pos","pos")]
write.table(magma_Input2,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_",i,"Input2.txt"),sep="\t",row.names=F,quote=F,col.names=F)
colnames(magma_Input1)<-c("SNP","P","N")
write.table(magma_Input1,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_",i,"Input1.txt"),sep="\t",row.names=F,quote=F)

cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
for i in Lymphocytecount3 
do
./magma --annotate window=10,10 --snp-loc /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_${i}Input2.txt \
--gene-loc /share/pub/dengcy/Singlecell/COVID19/MAGMA/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded \
--out /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10_down
done

cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
for i in Lymphocytecount3 
do
./magma --bfile /share/pub/dengcy/Singlecell/COVID19/MAGMA/g1000_eur/g1000_eur \
--pval /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/magma_${i}Input1.txt ncol=3 \
--gene-annot /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10_down.genes.annot \
--out /share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10down
done

cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
for i in Lymphocytecount3
do
int="/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/${i}annotated_10kbup_10down.genes.raw"
cell_type="/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/Hema_top10.txt"
./magma --gene-results  $int --set-annot  $cell_type --out /share/pub/dengcy/GWAS_Multiomics/compare/magma/${i}_magma_Hema
done
##################
library(org.Hs.eg.db)
library(dplyr)
magma_genes<-read.table(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/magma/predata/",i,"annotated_10kbup_10down.genes.out"),header=T)
#save(magma_genes,file=paste0(i,"_magma_genes.RData"))
g2s=toTable(org.Hs.egSYMBOL)
magma_genes$gene_id<-magma_genes$GENE
magma_genes=merge(magma_genes,g2s,by="gene_id",all.x=T)
magma_genes<-na.omit(magma_genes)
magma_genes<-magma_genes[order(magma_genes$P,decreasing=F),]
save(magma_genes,file=paste0(i,"_magma_genes.RData"))

```

barp[lot]

```R
 celltyperank<-c("01_HSC", "05_CMP.LMPP","06_CLP.1","15_CLP.2","07_GMP",
      "02_Early.Eryth" , "03_Late.Eryth" , "04_Early.Baso","08_GMP.Neut", 
      "09_pDC","10_cDC" ,"11_CD14.Mono.1","12_CD14.Mono.2","13_CD16.Mono","14_Unk",
     "16_Pre.B","17_B","18_Plasma",
      "19_CD8.N","20_CD4.N1","21_CD4.N2","22_CD4.M","23_CD8.EM","24_CD8.CM",
      "25_NK",
      "26_Unk")
celltypecolor<-c("FD5D5D","#FF8080","#FFC3C3","#FFCBCB",
    "#E26A2C","#FF6701","#FF8243","#FDA65D","#FFD07F",
    "#219F94","#C1DEAE",
    "#52734D","#91C788","#DDFFBC","#FEFFDE",
    "#F875AA","#FBACCC","#F1D1D0",
"#2F8F9D","#68A7AD","#3BACB6","#82DBD8","#B3E8E5","#CCF3EE",
    "#7C9473",
    "#AA8976")

library(ggplot2)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)
setwd("E:/OneDrive/GWAS_Multiomics/Compare/IntegrateVisulize")
######################
 
celltypecolor<-c("#CB181D","#D7403E","#E3695F","#EF9280",
                 "#CD9B1D","#EEB422","#F1C148","#FFE557","#FFF3B0",
                 "#2171B5","#73A5D2",
                 "#7F2704","#A95426","#BE6A37","#D38148",
                 "#3F007D","#653698","#8C6DB3",
                 "#0E4F27","#2A653F","#39704C","#558664","#80A789","#BAD3BB",
                 "#763857",
                 "#555555"
)

  #load(paste0(i,"_meta.data.RData")) 
  head(all_fortify_can)
  all_fortify_can$celltypes<-factor(all_fortify_can$celltypes,levels =celltyperank)

  p <- ggplot(all_fortify_can, aes(x = celltypes, 
                             y =scPagwastop10001,
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
  
##########scPagwas500_scdrs.raw_score点图
 a1<- tapply(all_fortify_can$scPagwastop10001, all_fortify_can$celltypes, function(x){
    mean(x)
  })
 a2<- tapply(all_fortify_can$scPagwastop5003, all_fortify_can$celltypes, function(x){
   mean(x)
 })
 a3<- tapply(all_fortify_can$magmatop10002, all_fortify_can$celltypes, function(x){
   mean(x)
 })
 a4<- tapply(all_fortify_can$magmatop5004, all_fortify_can$celltypes, function(x){
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

##########cell p value 百分比图
 a<- tapply(all_fortify_can$Cells.lm.adjp, all_fortify_can$celltypes, function(x){
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

```

## visualization for pbmc

### 1. PBMC Lymphocytecount

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

```
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

Monocyte

```


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

## Integrate Pheatmap

change the version of result ,v1.7 to b1.9

```R
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(stringr)
library(ggtext)
#setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/IntegrateVisulize")
 library(scPagwas)
 suppressMessages(library(Seurat))
library(org.Hs.eg.db)
library(dplyr)
traits<-c("RBCcount","Plateletcount","basophilcount","eosinophilcount" ,"Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")
lapply(traits[c(1,3,4,6:12)],function(i){
   load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.7.RData"))   
Pagwas<-rerun_Pagwas(Single_data=Single_data,Pagwas=Pagwas,assay="RNA",n_topgenes=1000) 
  save(Pagwas,file=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
})
```

### Construct a score matrix

```R

load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))  

score_df<-lapply(traits,function(i){    load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))   
return(Pagwas$scPagwas.topgenes.Score1)
})

score_df<-as.data.frame(score_df)
colnames(score_df)<-traits
score_df$celltypes<-Pagwas$celltypes
save(score_df,file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/BMMC_scPagwas_topgene_score_df.RData")

load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/BMMC_scPagwas_score_df.RData")

mean_scoredf<-tapply(1:nrow(score_df),factor(score_df$celltypes),function(x){
    return(colMeans(score_df[x,1:12]))
})
ct<-names(mean_scoredf)
mean_scoredf <- Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),mean_scoredf)
rownames(mean_scoredf)<-ct
save(mean_scoredf,file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/mean_scoredf.RData")

load("E:/OneDrive/GWAS_Multiomics/Compare/Hema_test2/mean_scoredf.RData")
library('ComplexHeatmap')
library(circlize)

mean_scoredf2<-scale(mean_scoredf)
col_fun = colorRamp2(c(-0.1002435,0.6713905), c("#F8ECD1", "#85586F"))
col_fun(seq(-0.1002435,0.6713905))
#magma_pbmc_p

a<-c("10_cDC","09_pDC","26_Unk","16_Pre.B","18_Plasma",
     "17_B","14_Unk","13_CD16.Mono","12_CD14.Mono.2","11_CD14.Mono.1",
     "25_NK","24_CD8.CM","22_CD4.M",
     "23_CD8.EM","19_CD8.N","20_CD4.N1","21_CD4.N2",
     "06_CLP.1","15_CLP.2","05_CMP.LMPP","07_GMP",
     "01_HSC","04_Early.Baso","02_Early.Eryth",
     "03_Late.Eryth","08_GMP.Neut")

b<-c("Lymphocytecount3","LymphocytePercent","Hemoglobinconcen",
     "MeanCorpuscularHemoglobin","MeanCorpusVolume",
     "WhiteBloodCellcount","neutrophilcount","eosinophilcount",
     "basophilcount","monocytecount")
mean_scoredf2<-mean_scoredf2[t(a),t(b)]
colnames(mean_scoredf2)[1]<-"Lymphocytecount"
pdf("E:/OneDrive/GWAS_Multiomics/Compare/Hema_test2/Figure_pagwas_bmmc_score_heatmap.pdf",height =7,width =6)
Heatmap(as.matrix(mean_scoredf2),
        col = col_fun,
        #left_annotation = ha, 
        cluster_columns = F,
        cluster_rows = F,
        color_space="HLS",
        border=T,
        row_gap = unit(0.25, "mm"),
        show_parent_dend_line=T,
        name = "mean score",
        #row_order =order(rdf$types),
        show_row_names=T,
        show_column_names=T
        #row_split=rdf$phenotypes
)
dev.off()
```



## sub-function

```R

rerun_Pagwas<-function(Single_data=Single_data,Pagwas,assay="RNA",n_topgenes=1000){
     class(Pagwas)<-"list"  
    scPagwas_topgenes <- names(Pagwas$gene_heritability_correlation[order(Pagwas$gene_heritability_correlation, decreasing = T), ])[1:n_topgenes]
    Single_data <- Seurat::AddModuleScore(Single_data, assay = assay, list(scPagwas_topgenes), name = c("scPagwas.topgenes.Score"))
      Pagwas[c(
      "VariableFeatures", "merge_scexpr","snp_gene_df",
      "rawPathway_list","data_mat"
    )] <- NULL

    message("* Get rankPvalue for each single cell")
    CellScalepValue <- rankPvalue(datS = t(data.matrix(GetAssayData(Single_data, assay = assay)[scPagwas_topgenes, ])), pValueMethod = "scale")
     Pagwas[c("snp_gene_df",
               "Pathway_sclm_results","CellScalepValue",
               "scPagwas.topgenes.Score"
               )] <- NULL
      #Pagwas[c()] <- NULL
      scPagwas_pathway <- SeuratObject::CreateAssayObject(data = t(data.matrix(Pagwas$Pathway_single_results)))
      scPagwas_pca <- SeuratObject::CreateAssayObject(data = Pagwas$pca_scCell_mat)

      Single_data[["scPagwasPaHeritability"]] <- scPagwas_pathway
      Single_data[["scPagwasPaPca"]] <- scPagwas_pca

      Single_data$scPagwas.lm.score <- Pagwas$scPagwas_score[rownames(Pagwas$Celltype_anno)]
      Single_data$CellScalepValue <- CellScalepValue[rownames(Pagwas$Celltype_anno), "pValueHighScale"]
      Single_data$CellScaleqValue <- CellScalepValue[rownames(Pagwas$Celltype_anno), "qValueHighScale"]
      Pagwas[ c("scPagwas.score","Celltype_anno",
                "Pathway_single_results","pca_scCell_mat")] <- NULL
   
      Single_data@misc<-Pagwas
    return(Single_data)
}


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

scPagwas_score_filter <- function(scPagwas_score) {
  # scPagwas_score <- scPagwas_score$scPagwas_score
  # remove the NAN!
  if (sum(is.nan(scPagwas_score)) > 0) {
    scPagwas_score[is.nan(scPagwas_score)] <- 0
  }
  # remove the inf values!
  if (Inf %in% scPagwas_score) {
    scPagwas_score[which(scPagwas_score == Inf)] <- max(scPagwas_score[-which(scPagwas_score == Inf)], na.rm = TRUE)
  }
  if (-Inf %in% scPagwas_score) {
    scPagwas_score[which(scPagwas_score == -Inf)] <- min(scPagwas_score[-which(scPagwas_score == -Inf)], na.rm = TRUE)
  }
  lower_bound <- quantile(scPagwas_score, 0.05, na.rm = TRUE)
  upper_bound <- quantile(scPagwas_score, 0.95, na.rm = TRUE)

  lower_ind <- which(scPagwas_score < lower_bound)
  upper_ind <- which(scPagwas_score > upper_bound)
  scPagwas_score[lower_ind] <- lower_bound
  scPagwas_score[upper_ind] <- upper_bound

  return(scPagwas_score)
}

#################另一种方式计算score

sclm_anotherscore<-function(Pagwas=Pagwas){
 a<-matrix(Pagwas@misc$sclm_results,ncol = 1)
Pagwas$sclm_score2<-rowSums(apply(GetAssayData(Pagwas,assay = "scPagwasPaPca"),1,function(x){
   x*a
  }),na.rm = T)
  data_mat<-GetAssayData(Pagwas,assay = "RNA")  
    sparse_lmcor <- corSparse(
      X = t(as_matrix(data_mat)),
      Y = data.matrix(Pagwas$sclm_score2)
    )
  rownames(sparse_lmcor) <- rownames(data_mat)
  colnames(sparse_lmcor) <- "allsnp_gene_heritability_correlation"
  sparse_lmcor[is.nan(sparse_lmcor)] <- 0
  #allsnp_gene_heritability_correlation <- sparse_lmcor
  
scPagwaslm_topgenes <- names(sparse_lmcor[order(sparse_lmcor, decreasing = T), ])[1:1000]

Pagwas <- Seurat::AddModuleScore(Pagwas, assay = "RNA", list(scPagwaslm_topgenes),
                                      name = c("scPagwas.lm2topgenes.Score"))
Pagwas$scPagwas.lm2topgenes.Score1 <- scPagwas_score_filter(scPagwas_score = Pagwas$scPagwas.lm2topgenes.Score1)
 return(Pagwas)
}

```

```

corSparse <- function(X, Y) {
  # X <-as_matrix(X)
  n <- nrow(X)
  muX <- colMeans(X)

  stopifnot(nrow(X) == nrow(Y))

  muY <- colMeans(Y)
  covmat <- (as.matrix(crossprod(X, Y)) - n * tcrossprod(muX, muY))
  sdvecX <- sqrt((colSums(X^2) - n * muX^2))
  sdvecY <- sqrt((colSums(Y^2) - n * muY^2))
  cormat <- covmat / tcrossprod(sdvecX, sdvecY)
  # cormat[is.nan(cormat),1]<-0
  return(cormat)
}


#' the source code from RMTstat package
#'
WishartMaxPar <- (function() {
  mu <- function( n,p ) {
    n.sqrt <- sqrt( n )
    p.sqrt <- sqrt( p )
    res    <- ( n.sqrt + p.sqrt )^2
    res
  }

  sigma <- function( n,p ) {
    n.sqrt <- sqrt( n )
    p.sqrt <- sqrt( p )
    res    <- ( n.sqrt + p.sqrt )*( 1/n.sqrt + 1/p.sqrt )^( 1/3 )
    res
  }

  mu.real <- function( n,p ) {
    mu( n-1/2,p-1/2 )
  }

  sigma.real <- function( n,p ) {
    sigma( n-1/2,p-1/2 )
  }

  alpha <- function( n,p ) {
    1/( 1 + ( mu( n-1/2,p+1/2 )/mu( n+1/2,p-1/2 ) )
        * sqrt( sigma( n-1/2,p+1/2 )/sigma( n+1/2,p-1/2 ) ) )
  }

  mu.cplx <- function( n,p ) {
    a   <- alpha( n,p )
    res <- mu( n-1/2,p+1/2 )*a + mu( n+1/2,p-1/2 )*( 1-a )
    res
  }

  sigma.cplx <- function( n,p ) {
    a   <- alpha( n,p )
    res <- sigma( n-1/2,p+1/2 )*a + sigma( n+1/2,p-1/2 )*( 1-a )
    res
  }

  function( ndf, pdim, var=1, beta=1 ) {
    n <- ndf
    p <- pdim

    if( beta == 1 ) {
      m <- mu.real( n,p )
      s <- sigma.real( n,p )
    } else if( beta == 2 ) {
      m <- mu.cplx( n,p )
      s <- sigma.cplx( n,p )
    } else {
      stop( "`beta' must be 1 or 2, not `", beta, "'")
    }

    center <- var*( m/n )
    scale  <- var*( s/n )

    list( centering=center, scaling=scale )
  }
})()



```

