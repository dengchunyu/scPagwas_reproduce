# Lymphocytecount Realdata_groundtruth

**Figure 3DEF**

**Supplementary Figure S4  BC**

**Supplementary Figure S3**

## 一， BMMC

### 1.Select the example celltypes for compare

we need a classical Lymphocyte celltypes (NK) and monocyte to compare the result.

```R
library(ggpubr)
library(ggplot2)
library(ggtext)

suppressMessages(library(Seurat))

##################
#LYMPHCYTE real groudtruth, monocyte,dc,t
###############
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("Pagwas_groundtruth.monocyte.RData")

load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas_v1.10.0.RData")
Pagwas_groundtruth2<-Pagwas_groundtruth
Pagwas_groundtruth<-Pagwas[,colnames(Pagwas_groundtruth)]

table(Pagwas_groundtruth$celltypes)
#colors_celltypes=c("#F38BA0","#FFBCBC","#EDF6E5","#B5EAEA")
Pagwas_groundtruth$types<-"non_Lymphocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="19_CD8.N"]<-"Lymphocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="20_CD4.N1"]<-"Lymphocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="25_NK"]<-"Lymphocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="17_B"]<-"Lymphocyte"



a1<-df[df$types=="non_Lymphocyte",]
a1<-a1[order(a1$scPagwas.TRS.Score1,decreasing = T)[1:5000],]
a2<-df[df$types=="Lymphocyte",]
a2<-a2[order(a2$scPagwas.TRS.Score1,decreasing = F)[1:5000],]


groundtruth_samples<-c(rownames(a2),rownames(a1))

df<-df[groundtruth_samples,]
df$sig<-1
df$sig[df$CellScaleqValue>0.01]<-0
table(df[,c("types","sig")])

Pagwas_groundtruth<-Pagwas[,groundtruth_samples]
####################################
###functions
###########################
umap_theme <- function(){
  theme_grey() %+replace%
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key=element_blank())
}
####################################
###2.plot tsne plot
###########################
saveRDS(Pagwas_groundtruth,file="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_groundtruth.monocyte.rds")

all_fortify_can <- fortify.Seurat.umap(Pagwas_groundtruth)
all_fortify_can2 <- fortify.Seurat.umap(Pagwas_groundtruth2)
all_fortify_can$UMAP_1<-all_fortify_can2$UMAP_1
all_fortify_can$UMAP_2<-all_fortify_can2$UMAP_2

```

### 2.TSNE plot for groundtruth example celltypes

```R
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =scPagwas.TRS.Score1), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas.TRS.Score1")

pdf(file="Umap.real_Lymphocytecount_scPagwas.TRS.Score.pdf",width = 7, height =7)
print(p1)
dev.off()


p2<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =magma1000_scdrs.zscore), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("magma1000_scdrs.zscore")

pdf(file="Umap.real_Lymphocytecount_magma1000_scdrs.zscore.pdf",width = 7, height = 7)
print(p2)
dev.off()

p3<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =scPagwas1000_scdrs.zscore), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas1000_scdrs.zscore")

pdf(file="Umap.real_Lymphocytecount_scPagwas1000_scdrs.zscore.pdf",width = 6, height = 6)
print(p3)
dev.off()


```

### 2.1 umap for other methods

```R
library(ggpubr)
library(ggplot2)
library(ggtext)
suppressMessages(library(Seurat))
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_groundtruth.monocyte.RData")
load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas_v1.10.0.RData")
scPagwas_genes<-names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation,decreasing=T),])

topgene<-list(topgenes=scPagwas_genes[1:1000])
Pagwas_groundtruth <- Seurat::AddModuleScore(Pagwas_groundtruth,assay = "RNA", topgene,name="scPagwas")

a<-Pagwas[,colnames(Pagwas_groundtruth)]
Pagwas_groundtruth@meta.data<-a@meta.data

sl<-lapply(c('smultixcan', 'spredixcan' ,'TWAS'),function(i){
files<-file.path(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{i}/scDRS_result/bmmcmodel/Lymphocytecount3_scDRSresult.csv'))
scDRS_re<-read.csv(files)
return((scDRS_re$X98))
})

Pagwas_groundtruth@meta.data$smultixcan_scDRS<-sl[[1]]
Pagwas_groundtruth@meta.data$spredixcan_scDRS<-sl[[2]]
Pagwas_groundtruth@meta.data$TWAS_scDRS<-sl[[3]]

df<-read.csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/MAGMA_scDRSresult_bmmcmodel_Lymphocytecount.csv")
Pagwas_groundtruth@meta.data$MAGMAtop1000<-df$MAGMAtop100
all_fortify_can <- fortify.Seurat.umap(Pagwas_groundtruth)

#MAGMA
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =scPagwas100), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas Seurat")

pdf(file="Umap.realbmmc_Lymphocytecount3_scPagwas_seurat.pdf",width = 7, height =7)
print(p1)
dev.off()

#MAGMA
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =MAGMAtop1000), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("MAGMA scDRS")

pdf(file="Umap.realbmmc_Lymphocytecount3_MAGMA_scDRS.pdf",width = 7, height =7)
print(p1)
dev.off()

#scPagwas
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =scPagwas1000_scdrs.raw_score), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas scDRS")

pdf(file="Umap.realbmmc_Lymphocytecount3_scPagwas_scDRS.pdf",width = 7, height =7)
print(p1)
dev.off()
#

#smultixcan
#all_fortify_can <- fortify.Seurat.umap(Pagwas_groundtruth)

p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =smultixcan_scDRS), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("smultixcan scDRS")

pdf(file="Umap.realbmmc_Lymphocytecount3_smultixcan_scDRS.pdf",width = 7, height =7)
print(p1)
dev.off()

###spredixcan
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =spredixcan_scDRS), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("spredixcan scDRS")

pdf(file="Umap.realbmmc_Lymphocytecount3_spredixcan_scDRS.pdf",width = 7, height =7)
print(p1)
dev.off()
###
###TWAS
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =TWAS_scDRS), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("TWAS scDRS")

pdf(file="Umap.realbmmc_Lymphocytecount3_TWAS_scDRS.pdf",width = 7, height =7)
print(p1)
dev.off()
```



### 3. Trait associated genes for 10 to 1000

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("Pagwas_groundtruth_bmmc_lymphocyte_v10.RData")

load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas_v1.10.0.RData")
oa<-order(Pagwas@misc$gene_heritability_correlation,decreasing = T)
topgene<-lapply(1:100, function(i){
  g1<-rownames(Pagwas@misc$gene_heritability_correlation)[oa[1:(10*i)]]
})

names(topgene)<-paste0(1:100,"times_topgene")
Pagwas_groundtruth <- Seurat::AddModuleScore(Pagwas_groundtruth, 
                                 assay = "RNA", 
                                 topgene,
                                 name="scPagwas")

Pagwas_groundtruth$types<-"non_Lymphocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="19_CD8.N"]<-"Lymphocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="20_CD4.N1"]<-"Lymphocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="25_NK"]<-"Lymphocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="17_B"]<-"Lymphocyte"

percent1<-lapply(paste0("scPagwas",1:100), function(i){
  a<-Pagwas_groundtruth@meta.data$types[order(Pagwas_groundtruth@meta.data[,i],decreasing = T)[1:4800]]
  return(sum(a=="Lymphocyte")/4800)
})

i<-"Lymphocytecount3"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_magma_genes.RData"))

ab<-lapply(1:100, function(i){
  b<-paste(magma_genes$symbol[1:(10*i)],collapse=",")
  return(b)
})

b<-paste0("MAGMAtop",1:100)
a<-data.frame(genes=unlist(ab))
rownames(a)<-b
write.csv(a,file="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Lymphocytecount/magma_scRDS_10_1000.csv")

library(SeuratDisk)
library(Seurat)
a<-GetAssayData(Pagwas_groundtruth,assay = "RNA",slot="counts")
Pagwas_readgroundtruth <- CreateSeuratObject(counts = a,
                                             project = "scPagwas",
                                             min.cells = 3,
                                             min.features = 200)


DefaultAssay(Pagwas_readgroundtruth) <- "RNA"
SaveH5Seurat(Pagwas_readgroundtruth, "Pagwas_groundtruth_addata.h5seurat")
Convert("Pagwas_groundtruth_addata.h5seurat", dest="h5ad")

```

### 4.Run scDRS for top 10- 1000genes for magma

```python
DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/realgroundtruth"
H5AD_FILE = os.path.join(DATA_PATH, "Pagwas_groundtruth_addata.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
scdrs.preprocess(adata)
df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Lymphocytecount/magma_scRDS_10_1000.csv", index_col=0)
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Lymphocytecount/geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Lymphocytecount/geneset.gs")


for i in range(1,101):
 a = 'MAGMAtop' + str(i)
 gene_list  = df_gs[a][0]
 gene_weight  = df_gs[a][1]
 df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=20)
 df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Lymphocytecount/"+a+".scDRS.csv", sep=",", index=False)

```

### 5.R run plot for graduate nunber of genes

```R


df<-read.csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/MAGMA_scDRSresult_bmmcmodel_Lymphocytecount.csv")
 percent2<-lapply(paste0("MAGMAtop",1:100), function(i){
   a<-Pagwas_groundtruth$types[order(df[,i],decreasing = T)[1:4800]]
   return(sum(a=="Lymphocyte")/4800)
 })
 
###smultixcan
setwd("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/smultixcan/scDRS_result/bmmcmodel")
ad<-read.csv("Lymphocytecount3_scDRSresult.csv")
colnames(ad)<-paste0("smultixcantop",1:100)
#Pagwas_groundtruth@meta.data<-cbind(Pagwas_groundtruth@meta.data,ad)
percent3<-lapply(paste0("smultixcantop",1:100), function(i){
  a<-Pagwas_groundtruth@meta.data$types[order(ad[,i],decreasing = T)[1:4800]]
  return(sum(a=="Lymphocyte")/4800)
})
#spredixcan
setwd("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/spredixcan/scDRS_result/bmmcmodel")
ad<-read.csv("Lymphocytecount3_scDRSresult.csv")
colnames(ad)<-paste0("spredixcantop",1:100)
#Pagwas_groundtruth@meta.data<-cbind(Pagwas_groundtruth@meta.data,ad)
percent4<-lapply(paste0("spredixcantop",1:100), function(i){
  a<-Pagwas_groundtruth@meta.data$types[order(ad[,i],decreasing = T)[1:4800]]
  return(sum(a=="Lymphocyte")/4800)
})
#TWAS
setwd("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/TWAS/scDRS_result/bmmcmodel")
ad<-read.csv("Lymphocytecount3_scDRSresult.csv")
colnames(ad)<-paste0("TWAStop",1:100)
#Pagwas_groundtruth@meta.data<-cbind(Pagwas_groundtruth@meta.data,ad)
percent5<-lapply(paste0("TWAStop",1:100), function(i){
  a<-Pagwas_groundtruth@meta.data$types[order(ad[,i],decreasing = T)[1:4800]]
  return(sum(a=="Lymphocyte")/4800)
})
###########
#可视化

 library(ggplot2)
 library(ggpubr)
 gg1<-data.frame(Number=c(1:100 *10,1:100 *10,1:100 *10,1:100 *10,1:100 *10),
                 Percent=c(unlist(percent1),unlist(percent2),unlist(percent3),unlist(percent4),unlist(percent5)),

                 Types=c(rep("scPagwas",100),rep("MAGMA",100),rep("smultixcan",100),rep("spredixcan",100),rep("TWAS",100))
 )
gg1$Types<-factor(gg1$Types,levels=c("scPagwas","MAGMA","TWAS","spredixcan","smultixcan"))
 setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
 pdf("Lymphocyte_bmmc_real_percent.pdf")
 ggplot(gg1,aes(Number,Percent,group=Types, color=Types))+
   geom_point()+
   geom_line(position = position_dodge(0.1),cex=1)+
scale_color_manual(values=c('#e8505b','#f9d56e','#f3ecc2','#14b1ab','#d4b5b0'))+
scale_fill_manual(values=c('#e8505b','#f9d56e','#f3ecc2','#14b1ab','#d4b5b0'))+
   labs(x="Number of genes", y="Percent of Monocyte cells in top half cells") +
   theme_classic()
 dev.off()

save(Pagwas_groundtruth,file="Pagwas_groundtruth_bmmc_lymphocyte_v10.RData")
```

### 6.percent plot

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("Pagwas_groundtruth_bmmc_lymphocyte_v10.RData")
all_fortify_can<-fortify.Seurat.umap(Pagwas_groundtruth)
percent_plot<-function(df1){
        n<-nrow(df1)
        b<-rep("top",n)
        #b[1:(0.48*n)]<-1
        b[(0.48*n+1):n]<-"bottom"

        df1$rg<-b
        #df1$type
        print(table(data.frame(df1$types2,df1$rg)))
        t1<-table(data.frame(df1$types2,df1$rg))
        percent1<-t1[,2]/sum(t1[,2])
print(percent1)
       # df <- data.frame(
        #  group =c("monocyte","non_monocyte"),
        #  value = percent1)

       # p1<-ggdonutchart(df, "value", label = "group",
       #                  fill = "group", color = "white",
        #                 palette = c("#DF7861", "#D4E2D4") )
        return(p1)
}


df1<-all_fortify_can[order(all_fortify_can$scPagwas100,decreasing=T),]
p1<-percent_plot(df1)
#    Lymphocyte non_Lymphocyte
#   0.996666667    0.003333333


df1<-all_fortify_can[order(all_fortify_can$scPagwas1000_scdrs.zscore,decreasing=T),]
p1<-percent_plot(df1)
# monocyte non_monocyte
#      0.9836       0.0164
df1<-all_fortify_can[order(all_fortify_can$MAGMAtop1000,decreasing=T),]
p2<-percent_plot(df1)
#    Lymphocyte non_Lymphocyte
#     0.4235417      0.5764583
df1<-all_fortify_can[order(all_fortify_can$smultixcan_scDRS,decreasing=T),]
p3<-percent_plot(df1)
#    Lymphocyte non_Lymphocyte
 #     0.531875       0.468125
df1<-all_fortify_can[order(all_fortify_can$spredixcan_scDRS,decreasing=T),]
p4<-percent_plot(df1)
#Lymphocyte non_Lymphocyte
#          0.39           0.61
df1<-all_fortify_can[order(all_fortify_can$TWAS_scDRS,decreasing=T),]
p5<-percent_plot(df1)
#  Lymphocyte non_Lymphocyte
#     0.3147917      0.6852083
pdf("percent.methods_monocyte_realBMMC_scdrs.pdf",width = 10,height = 8)
ggpubr::ggarrange(p1,p2,p3,p4,p5,nrow = 2)
dev.off()

save(Pagwas_groundtruth,file="Pagwas_groundtruth_bmmc_lymphocyte_v10.RData")
```



## 二、PBMC

### 1.Select the example celltypes for compare

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
 
 load("/share/pub/dengcy/GWAS_Multiomics/compare/Pbmc/monocytecount_pbmc_scPagwas_singlecell.RData")
 
 celltypes<-c("CD14","CD4","NK_56hi")
 
 #colors_celltypes=c("#F38BA0","#FFBCBC","#EDF6E5","#B5EAEA")
 #Pagwas<-Pagwas[,Pagwas$initial_clustering %in% celltypes]
 df<-Pagwas@meta.data
 groundtruth_samples<-c(rownames(df)[sample(which(df$initial_clustering=="CD14"),5000)],
                        rownames(df)[sample(which(df$initial_clustering=="CD4"),3000)],
                        rownames(df)[sample(which(df$initial_clustering=="NK_56hi"),2000)])
 Pagwas_pbmc_groundtruth<-Pagwas[,groundtruth_samples]
 rm(Pagwas)
 Pagwas_pbmc_groundtruth$types<-"monocyte"
 Pagwas_pbmc_groundtruth$types[df$initial_clustering=="CD4"]<-"non_monocyte"
 Pagwas_pbmc_groundtruth$types[df$initial_clustering=="NK_56hi"]<-"non_monocyte"
 
 save(Pagwas_pbmc_groundtruth,file="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_pbmc_realgroundtruth.RData")

 load("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_pbmc_realgroundtruth.RData")
 
 Pagwas_pbmc_groundtruth <- FindVariableFeatures(Pagwas_pbmc_groundtruth,nfeatures = 3000)
 Pagwas_pbmc_groundtruth <- RunPCA(object = Pagwas_pbmc_groundtruth, assay = "RNA", npcs = 50)
 Pagwas_pbmc_groundtruth<-cluster_pca_umap(obj = Pagwas_pbmc_groundtruth,reduction="pca",cluster_res = 0.3)
 ######计算scPagwas
 Idents(Pagwas_pbmc_groundtruth)<-as.vector(Pagwas_pbmc_groundtruth$initial_clustering)
 
 saveRDS(Pagwas_pbmc_groundtruth,file="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_pbmc_realgroundtruth.rds")
 
 i<-"Lymphocytecount3"
 Pagwas_pbmc_groundtruth<-scPagwas_main(Pagwas = NULL,
                       gwas_data=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_gwas_data.txt"),
                       Single_data =Pagwas_pbmc_groundtruth,
                       output.prefix="real_pbmc_Lymphocytecount",
                       output.dirs="real_pbmc_Lymphocytecount",
                       block_annotation = block_annotation,
                       assay="RNA",
                       Pathway_list=Genes_by_pathway_kegg,
                       chrom_ld = chrom_ld,
                       singlecell=T,
                       seruat_return=T,
                       celltype=F,
                       ncores = 1)
 save(Pagwas_pbmc_groundtruth,file = "Pagwas_pbmc_groundtruth_Lymphocytecount.RData")
 
 load("Pagwas_pbmc_groundtruth_Lymphocytecount.RData")
 color_scanpy_patient <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")
 pdf("real_PBMC_groundtruth_monocyte.pdf",width = 6,height = 6)
 DimPlot(Pagwas_pbmc_groundtruth,reduction ="umap",
         group.by = "initial_clustering",pt.size=0.3,
         label = F, repel=TRUE)+ 
   umap_theme()+labs(x="TSNE",y="")+
   scale_colour_manual(name = "celltypes", values =color_scanpy_patient[c(2:9,1)]) +
   theme(aspect.ratio=1)
 dev.off()


```

### 2. Get top genes for magma and scPagwas for scDRS

```R
 ######################
 i<-"Lymphocytecount3"
 load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_magma_genes.RData"))
 
 scPagwas_genes<-names(Pagwas_pbmc_groundtruth@misc$gene_heritability_correlation[order(Pagwas_pbmc_groundtruth@misc$gene_heritability_correlation,decreasing=T),])
 
 magmatop1000<-intersect(magma_genes$symbol,rownames(Pagwas_pbmc_groundtruth))[1:1000]
 scPagwastop1000<-scPagwas_genes[1:1000]
 magmatop500<-intersect(magma_genes$symbol,rownames(Pagwas_pbmc_groundtruth))[1:500]
 scPagwastop500<-scPagwas_genes[1:500]
 topgene<-list(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500)
 names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
 
 a<-data.frame(genes=c(paste(scPagwastop1000,collapse=","),paste(magmatop1000,collapse=","),paste(scPagwastop500,collapse=","),paste(magmatop500,collapse=",")))
 rownames(a)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
 write.csv(a,file=paste0("scPagwas.magma.PBMC.topgenes.",i,".csv"))

```

PYTHON for scDRS for top 1000 genes

```python
df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/scPagwas.magma.PBMC.topgenes.Lymphocytecount3.csv", index_col=0)
df_gs = df_gs.loc[["scPagwastop1000","magmatop1000","scPagwastop500","magmatop500"],:]
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/geneset.gs")
 
gene_list  = df_gs['magmatop1000'][0]
gene_weight  = df_gs['magmatop1000'][1]
df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=20)
df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Lympho.pbmc.1000magma.df_model.csv", sep=",", index=False)
 

gene_list  = df_gs['scPagwastop1000'][0]
gene_weight  = df_gs['scPagwastop1000'][1]
df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=20)
df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Lympho.pbmc.1000scPagwas.df_model2.csv", sep=",", index=False)

```

### 3.TSNE plot for groundtruth example celltypes

```R
 
 library(ggpubr)
 library(ggplot2)
 library(ggtext)
 all_fortify_can <- fortify.Seurat.umap(Pagwas_pbmc_groundtruth)
 p1<- ggplot() +
   geom_point(data = all_fortify_can,
              aes(x = UMAP_1, y = UMAP_2,color =scPagwas.TRS.Score1), size = 0.3, alpha = 1) +
   umap_theme() +
   scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
   theme(aspect.ratio=1) +
   theme(legend.text=element_markdown(size=14),
         legend.title=element_text(size=14)) +
   guides(colour = guide_legend(override.aes = list(size=3)))  +
   ggtitle("scPagwas.TRS.Score1")
 
 pdf(file="Umap.pbmc_Lymphocytecount_scPagwas.TRS.Score.pdf",width = 7, height =7)
 print(p1)
 dev.off()
 
 
 p1<- ggplot() +
   geom_point(data = all_fortify_can,
              aes(x = UMAP_1, y = UMAP_2,color =scPagwas.TRS.Score1), size = 0.3, alpha = 1) +
   umap_theme() +
   scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
   theme(aspect.ratio=1) +
   theme(legend.text=element_markdown(size=14),
         legend.title=element_text(size=14)) +
   guides(colour = guide_legend(override.aes = list(size=3)))  +
   ggtitle("scPagwas.TRS.Score1")
 
 pdf(file="Umap.pbmc_Lymphocytecount_scPagwas.TRS.Score.pdf",width = 7, height =7)
 print(p1)
 dev.off()

```

### 3.1umap for other methods

```R
library(ggpubr)
library(ggplot2)
library(ggtext)
suppressMessages(library(Seurat))
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("Pagwas_pbmc_groundtruth_Lymphocytecount.RData")

sl<-lapply(c('smultixcan', 'spredixcan' ,'TWAS'),function(i){
files<-file.path(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{i}/scDRS_result/pbmcmodel/Lymphocytecount3_scDRSresult.csv'))
scDRS_re<-read.csv(files)
return((scDRS_re$X98))
})

Pagwas_pbmc_groundtruth@meta.data$smultixcan_scDRS<-sl[[1]]
Pagwas_pbmc_groundtruth@meta.data$spredixcan_scDRS<-sl[[2]]
Pagwas_pbmc_groundtruth@meta.data$TWAS_scDRS<-sl[[3]]

#smultixcan
all_fortify_can <- fortify.Seurat.umap(Pagwas_pbmc_groundtruth)

p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =smultixcan_scDRS), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("smultixcan scDRS")

pdf(file="Umap.pbmcreal_Lymphocytecount3_smultixcan_scDRS.pdf",width = 7, height =7)
print(p1)
dev.off()

###spredixcan
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =spredixcan_scDRS), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("smultixcan scDRS")

pdf(file="Umap.pbmcreal_Lymphocytecount3_spredixcan_scDRS.pdf",width = 7, height =7)
print(p1)
dev.off()
###
###TWAS
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =TWAS_scDRS), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("smultixcan scDRS")

pdf(file="Umap.pbmcreal_Lymphocytecount3_TWAS_scDRS.pdf",width = 7, height =7)
print(p1)
dev.off()
```



### 4.Trait associated genes for 10 to 1000

```R
i<-"Lymphocytecount3"
 load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_magma_genes.RData"))
 
 ab<-lapply(1:100, function(i){
   b<-paste(magma_genes$symbol[1:(10*i)],collapse=",")
   return(b)
 })
 b<-paste0("MAGMAtop",1:100)
 a<-data.frame(genes=unlist(ab))
 rownames(a)<-b
 write.csv(a,file="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/PBMC_Lymphocytecount3_magma_scRDS_10_1000.csv")
 
 setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
 load("Pagwas_pbmc_groundtruth_Lymphocytecount.RData")
 library(SeuratDisk)
 library(Seurat)
 a<-GetAssayData(Pagwas_pbmc_groundtruth,assay = "RNA",slot="counts")
 Pagwas_pbmc_groundtruth <- CreateSeuratObject(counts = a,
                                               project = "scPagwas",
                                               min.cells = 3,
                                               min.features = 200)
 DefaultAssay(Pagwas_pbmc_groundtruth) <- "RNA"
 SaveH5Seurat(Pagwas_pbmc_groundtruth, "Pagwas_pbmc_groundtruth.h5seurat")
 Convert("Pagwas_pbmc_groundtruth.h5seurat", dest="h5ad")
```

### 5.scDRS python

```python
 DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/"
 H5AD_FILE = os.path.join(DATA_PATH, "Pagwas_pbmc_groundtruth.h5ad")
 adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
 scdrs.preprocess(adata)
 
df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/PBMC_Lymphocytecount3_magma_scRDS_10_1000.csv", index_col=0)
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/geneset.gs")
 
for i in range(1,101):
 a = 'MAGMAtop' + str(i)
 gene_list  = df_gs[a][0]
 gene_weight  = df_gs[a][1]
 df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=20)
 df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/PBMC/Lymphocytecount."+a+".scDRS.csv", sep=",", index=False)

```

### 6.R run plot for graduate nunber of genes

```R

load("Pagwas_pbmc_groundtruth_Lymphocytecount.RData")
#Pagwas_pbmc_groundtruth@meta.data
oa<-order(Pagwas_pbmc_groundtruth@misc$gene_heritability_correlation,decreasing = T)
 topgene<-lapply(1:100, function(i){
   g1<-rownames(Pagwas_pbmc_groundtruth@misc$gene_heritability_correlation)[oa[1:(10*i)]]
 })
 names(topgene)<-paste0(1:100,"times_topgene")
 Pagwas_pbmc_groundtruth <- Seurat::AddModuleScore(Pagwas_pbmc_groundtruth, 
                                                   assay = "RNA", 
                                                   topgene,
                                                   name="scPagwas")
 
 Pagwas_pbmc_groundtruth$types<-"Lymphocyte"
 Pagwas_pbmc_groundtruth$types[Pagwas_pbmc_groundtruth$initial_clustering=="CD14"]<-"non_Lymphocyte"
 
 percent1<-lapply(paste0("scPagwas",1:100), function(i){
   a<-Pagwas_pbmc_groundtruth@meta.data$types[order(Pagwas_pbmc_groundtruth@meta.data[,i],decreasing = T)[1:5000]]
   return(sum(a=="Lymphocyte")/5000)
 })
 
 
df<-read.csv("MAGMA_scDRSresult_pbmcmodel_Lymphocytecount3.csv")
percent2<-lapply(paste0("MAGMAtop",1:100), function(i){
   a<-Pagwas_pbmc_groundtruth@meta.data$types[order(df[,i],decreasing = T)[1:5000]]
   return(sum(a=="Lymphocyte")/5000)
 })
###smultixcan
setwd("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/smultixcan/scDRS_result/pbmcmodel")
ad<-read.csv("Lymphocytecount3_scDRSresult.csv")
colnames(ad)<-paste0("smultixcantop",1:100)
#Pagwas_pbmc_groundtruth@meta.data<-cbind(Pagwas_pbmc_groundtruth@meta.data,ad)
percent3<-lapply(paste0("smultixcantop",1:100), function(i){
  a<-Pagwas_pbmc_groundtruth@meta.data$types[order(ad[,i],decreasing = T)[1:5000]]
  return(sum(a=="Lymphocyte")/5000)
})
#spredixcan
setwd("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/spredixcan/scDRS_result/pbmcmodel")
ad<-read.csv("Lymphocytecount3_scDRSresult.csv")
colnames(ad)<-paste0("spredixcantop",1:100)
#Pagwas_pbmc_groundtruth@meta.data<-cbind(Pagwas_pbmc_groundtruth@meta.data,ad)
percent4<-lapply(paste0("spredixcantop",1:100), function(i){
  a<-Pagwas_pbmc_groundtruth@meta.data$types[order(ad[,i],decreasing = T)[1:5000]]
  return(sum(a=="Lymphocyte")/5000)
})
#TWAS
setwd("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/TWAS/scDRS_result/pbmcmodel")
ad<-read.csv("Lymphocytecount3_scDRSresult.csv")
colnames(ad)<-paste0("TWAStop",1:100)
#Pagwas_pbmc_groundtruth@meta.data<-cbind(Pagwas_pbmc_groundtruth@meta.data,ad)
percent5<-lapply(paste0("TWAStop",1:100), function(i){
  a<-Pagwas_pbmc_groundtruth@meta.data$types[order(ad[,i],decreasing = T)[1:5000]]
  return(sum(a=="Lymphocyte")/5000)
}) 
 
 library(ggplot2)
 library(ggpubr)
 gg1<-data.frame(Number=c(1:100 *10,1:100 *10,1:100 *10,1:100 *10,1:100 *10),
                 Percent=c(unlist(percent1),unlist(percent2),unlist(percent3),unlist(percent4),unlist(percent5)),

                 Types=c(rep("scPagwas",100),rep("MAGMA",100),rep("smultixcan",100),rep("spredixcan",100),rep("TWAS",100))
 )
gg1$Types<-factor(gg1$Types,levels=c("scPagwas","MAGMA","TWAS","spredixcan","smultixcan"))
 setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
 pdf("Lymphocyte_realPBMC_percent.pdf")
 ggplot(gg1,aes(Number,Percent,group=Types, color=Types))+
   geom_point()+
   geom_line(position = position_dodge(0.1),cex=1)+
scale_color_manual(values=c('#e8505b','#f9d56e','#f3ecc2','#14b1ab','#d4b5b0'))+
scale_fill_manual(values=c('#e8505b','#f9d56e','#f3ecc2','#14b1ab','#d4b5b0'))+
   labs(x="Number of genes", y="Percent of Lymphocyte cells in top half cells") +
   theme_classic()
 dev.off()

```

### 6.Rankplot

```R
library('ComplexHeatmap')
library(circlize)
####
colors_celltypes=c("#DF7861","#8CC0DE")
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("Pagwas_groundtruth.monocyte.RData")
load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/monocytecount_bmmc_Pagwas1.10.0.RData")
Pagwas_groundtruth$scPagwas.TRS.Score1<-Pagwas@meta.data[colnames(Pagwas_groundtruth),"scPagwas.TRS.Score1"]

Pagwas_groundtruth$types<-"monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="09_pDC"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="10_cDC"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="19_CD8.N"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="20_CD4.N1"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="25_NK"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="17_B"]<-"non_monocyte"
all_fortify_can<-fortify.Seurat.umap(Pagwas_groundtruth)
all_fortify_can<-all_fortify_can[order(all_fortify_can$scPagwas.TRS.Score1,decreasing=T),]

pdf("BMMCmodel_monocyte_scPagwas_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "scPagwas rank"))
dev.off()


df<-read.csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/MAGMAtop100.scDRS.csv")
Pagwas_groundtruth$MAGMA<-df$zscore

all_fortify_can<-fortify.Seurat.umap(Pagwas_groundtruth)
all_fortify_can<-all_fortify_can[order(all_fortify_can$MAGMA,decreasing=T),]

pdf("BMMCmodel_monocyte_magma_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "MAGMA rank"))
dev.off()

###smultixcan

ad<-read.csv("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/smultixcan/scDRS_result/bmmcmodel/monocytecount_scDRSresult.csv")
colnames(ad)<-paste0("smultixcantop",1:100)
Pagwas_groundtruth@meta.data$Smultixcan<-ad$smultixcantop100
all_fortify_can<-fortify.Seurat.umap(Pagwas_groundtruth)
all_fortify_can<-all_fortify_can[order(all_fortify_can$Smultixcan,decreasing=T),]
pdf("BMMCmodel_monocyte_Smultixcan_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "Smultixcan rank"))
dev.off()

#spredixcan
ad<-read.csv("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/spredixcan/scDRS_result/bmmcmodel/monocytecount_scDRSresult.csv")
colnames(ad)<-paste0("spredixcantop",1:100)
Pagwas_groundtruth@meta.data$Spredixcan<-ad$spredixcantop100
all_fortify_can<-fortify.Seurat.umap(Pagwas_groundtruth)
all_fortify_can<-all_fortify_can[order(all_fortify_can$Spredixcan,decreasing=T),]
pdf("BMMCmodel_monocyte_Spredixcan_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "Spredixcan rank"))
dev.off()


#TWAS
ad<-read.csv("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/TWAS/scDRS_result/bmmcmodel/monocytecount_scDRSresult.csv")
colnames(ad)<-paste0("TWAStop",1:100)
Pagwas_groundtruth@meta.data$TWAS<-ad$TWAStop100

all_fortify_can<-fortify.Seurat.umap(Pagwas_groundtruth)
all_fortify_can<-all_fortify_can[order(all_fortify_can$TWAS,decreasing=T),]
pdf("BMMCmodel_monocyte_TWAS_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "TWAS rank"))
dev.off()
save(Pagwas_groundtruth,file='Pagwas_groundtruth_bmmc_monocyte_v10.RData')
```



# Monocytecount Realdata_groundtruth

## 一， BMMC

### 1.Select the example celltypes for compare

```R
library(ggpubr)
library(ggplot2)
library(ggtext)

suppressMessages(library(Seurat))

##################
#GET EXAMPLE
###############
setwd("D:/OneDrive/GWAS_Multiomics/realgroundtruth")
load("D:/OneDrive/GWAS_Multiomics/Compare/Hema_test2/monocytecount_bmmc_Pagwas1.10.0.RData")
table(Pagwas$celltypes)
meta.data<-Pagwas@meta.data
celltypes<-c("11_CD14.Mono.1","12_CD14.Mono.2","13_CD16.Mono","09_pDC","10_cDC",
             "17_B","19_CD8.N","20_CD4.N1","25_NK")

#colors_celltypes=c("#F38BA0","#FFBCBC","#EDF6E5","#B5EAEA")
df<-meta.data[meta.data$celltypes %in% celltypes,]
df$types<-"monocyte"
df$types[df$celltypes=="09_pDC"]<-"non_monocyte"
df$types[df$celltypes=="10_cDC"]<-"non_monocyte"
df$types[df$celltypes=="19_CD8.N"]<-"non_monocyte"
df$types[df$celltypes=="20_CD4.N1"]<-"non_monocyte"
df$types[df$celltypes=="25_NK"]<-"non_monocyte"
df$types[df$celltypes=="17_B"]<-"non_monocyte"

groundtruth_samples<-c(rownames(df)[sample(which(df$celltypes=="09_pDC"),100)],
                       rownames(df)[sample(which(df$celltypes=="10_cDC"),100)],
                       rownames(df)[sample(which(df$celltypes=="19_CD8.N"),1500)],
                       rownames(df)[sample(which(df$celltypes=="17_B"),1000)],
                       rownames(df)[sample(which(df$celltypes=="20_CD4.N1"),1500)],
                       rownames(df)[sample(which(df$celltypes=="25_NK"),800)],
  rownames(df)[sample(which(df$types=="monocyte"),5000)])

df<-df[groundtruth_samples,]
df$sig<-1
df$sig[df$CellScaleqValue>0.01]<-0
table(df[,c("types","sig")])

Pagwas_groundtruth<-Pagwas[,groundtruth_samples]

Pagwas_groundtruth <- FindVariableFeatures(Pagwas_groundtruth,nfeatures = 3000)
Pagwas_groundtruth <- RunPCA(object = Pagwas_groundtruth, assay = "RNA", npcs = 50)

##subfunction:
cluster_pca_umap <- function(obj,assay=NULL, reduction,cluster_res = 0.3){
  #obj2 <- RunPCA(obj, assay = "SCT", reduction = "harmony",verbose = F)
  obj2 <- RunTSNE(object = obj,assay = assay, reduction = reduction, dims = 1:50,check_duplicates = FALSE)
  obj2 <- RunUMAP(object = obj2, assay =assay, reduction = reduction, dims = 1:50,check_duplicates = FALSE)
  obj2 <- FindNeighbors(object=obj2, assay = assay, reduction = reduction, dims = 1:50)
  obj2 <- FindClusters(object=obj2, resolution = cluster_res)
  return(obj2)
}
umap_theme <- function(){
  theme_grey() %+replace%
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.1),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key=element_blank())
}
fortify.Seurat.umap <- function(x) {
  xy1 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "umap")
  )
  colnames(xy1) <- c("UMAP_1", "UMAP_2")
  xy1$UMAP_1 <- as.numeric(xy1$UMAP_1)
  xy1$UMAP_2 <- as.numeric(xy1$UMAP_2)
  
  return(cbind(xy1, as.data.frame(x@meta.data)))
}

load("Pagwas_groundtruth.monocyte.RData")
Pagwas_groundtruth<-cluster_pca_umap(obj = Pagwas_groundtruth,reduction="pca",cluster_res = 0.3)
save(Pagwas_groundtruth,file = "Pagwas_groundtruth.monocyte.RData")


color_scanpy_patient <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")
pdf("groundtruth_monocyte_test8.13.pdf",width = 6,height = 6)
 DimPlot(Pagwas_groundtruth,reduction ="umap",
                         group.by = "celltypes",pt.size=0.2,
                         label = F, repel=TRUE)+ 
  umap_theme()+labs(x="TSNE",y="")+
  scale_colour_manual(name = "celltypes", values =color_scanpy_patient[c(2:9,1)]) +
  theme(aspect.ratio=1)
dev.off()

```

### 2.TSNE plot for groundtruth example celltypes

```R
all_fortify_can <- fortify.Seurat.umap(Pagwas_groundtruth)
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =Cluster100), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas.TRS.Score1")

pdf(file="Umap.monocytecount_scPagwas.TRS.Score.pdf",width = 6, height = 6)
print(p1)
dev.off()


p2<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =MAGMAtop100.scdrs), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("magma1000_scdrs.zscore")

pdf(file="Umap.monocytecount_magma1000_scdrs.zscore.pdf",width = 6, height = 6)
print(p2)
dev.off()

p3<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =scPagwas1000_scdrs.zscore), size = 0.2, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas1000_scdrs.zscore")

pdf(file="Umap.monocytecount_scPagwas1000_scdrs.zscore.pdf",width = 6, height = 6)
print(p3)
dev.off()

```



### 2.1umap for other methods

```
library(ggpubr)
library(ggplot2)
library(ggtext)
suppressMessages(library(Seurat))
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_groundtruth.monocyte.RData")

scDRS_re<-read.csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/ThreeMethods_scDRSresult_bmmcmodel.csv")

Pagwas_groundtruth@meta.data$smultixcan_top100_monocytecount<-scDRS_re$smultixcan_top100_monocytecount
Pagwas_groundtruth@meta.data$spredixcan_top100_monocytecount<-scDRS_re$spredixcan_top100_monocytecount
Pagwas_groundtruth@meta.data$TWAS_top100_monocytecount<-scDRS_re$TWAS_top100_monocytecount

df<-read.csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/mono.bmmc.1000magma.df_model.csv")
Pagwas_groundtruth$magma_monocytecount_scDRS<-df$norm_score
df<-read.csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/mono.bmmc.1000scPagwas.df_model.csv")
Pagwas_groundtruth$scPagwas_monocytecount_scDRS<-df$norm_score

#smultixcan
all_fortify_can <- fortify.Seurat.umap(Pagwas_groundtruth)
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =magma_monocytecount_scDRS), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("magma scDRS")

pdf(file="Umap.real_monocytecount_magma_scDRS.pdf",width = 7, height =7)
print(p1)
dev.off()

p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =scPagwas_monocytecount_scDRS), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scPagwas scDRS")

pdf(file="Umap.real_monocytecount_scPagwas_scDRS.pdf",width = 7, height =7)
print(p1)
dev.off()

p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =smultixcan_top100_monocytecount), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("smultixcan scDRS")

pdf(file="Umap.real_monocytecount_smultixcan_scDRS.pdf",width = 7, height =7)
print(p1)
dev.off()

###spredixcan
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =spredixcan_top100_monocytecount), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("spredixcan scDRS")

pdf(file="Umap.real_monocytecount_spredixcan_scDRS.pdf",width = 7, height =7)
print(p1)
dev.off()
###
###TWAS
p1<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =TWAS_top100_monocytecount), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("TWAS scDRS")

pdf(file="Umap.real_monocytecount_TWAS_scDRS.pdf",width = 7, height =7)
print(p1)
dev.off()
```



### 3,Heatmap and ranked plot for scPagwas and magma

```R
library('ComplexHeatmap')
library(circlize)
####
colors_celltypes=c("#DF7861","#8CC0DE")

all_fortify_can<-all_fortify_can[order(all_fortify_can$scPagwas.TRS.Score1,decreasing=T),]
all_fortify_can$celltype<-as.vector(all_fortify_can$celltype)
pdf("rankplot_scPagwas.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "scPagwas rank"))
dev.off()

all_fortify_can<-all_fortify_can[order(all_fortify_can$MAGMAtop100.scdrs,decreasing=T),]
all_fortify_can$types<-as.vector(all_fortify_can$types)
pdf("rankplot_magma.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "MAGMA rank"))
dev.off()

############
#proportion
##############

df1<-all_fortify_can[order(all_fortify_can$scPagwas.TRS.Score1,decreasing=T),]
#a<-scPagwas_magma_genelist[[i]]
n<-nrow(df1)
b<-rep(1,n)
b[1:(0.5*n)]<-1
b[(0.5*n+1):(1*n)]<-2

df1$rg<-b

print(table(data.frame(df1$types,df1$rg)))
t1<-table(data.frame(df1$types,df1$rg))
percent1<-t1[,1]/sum(t1[,1])
#   Correct         non 
#0.8844 0.1156

a <- data.frame(
  group = names(percent1),
  value = percent1)

p1<-ggdonutchart(a, "value", label = rev(a$value),
                 fill = "group", color = "white",
                 palette = c("#DF7861","#8CC0DE") )

t1<-table(data.frame(df1$types,df1$rg))
percent1<-t1[,2]/sum(t1[,2])

a <- data.frame(
  group = names(percent1),
  value = percent1)

p2<-ggdonutchart(a, "value", label = rev(a$value),
                 fill = "group", color = "white",
                 palette = c("#DF7861","#8CC0DE"))

####################magma
df1<-all_fortify_can[order(all_fortify_can$MAGMAtop100.scdrs,decreasing=T),]
#a<-scPagwas_magma_genelist[[i]]
n<-nrow(df1)
b<-rep(1,n)
b[1:(0.5*n)]<-1
b[(0.5*n+1):(1*n)]<-2

df1$rg<-b

print(table(data.frame(df1$types,df1$rg)))
t1<-table(data.frame(df1$types,df1$rg))
percent1<-t1[,1]/sum(t1[,1])
a <- data.frame(
  group = names(percent1),
  value = percent1)
p3<-ggdonutchart(a, "value", label = rev(a$value),
                 fill = "group", color = "white",
                 palette = c("#DF7861","#8CC0DE") )

#t1<-table(data.frame(df1$types,df1$rg))
percent1<-t1[,2]/sum(t1[,2])
a <- data.frame(
  group = names(percent1),
  value = percent1)
p4<-ggdonutchart(a, "value", label = rev(a$value),
                 fill = "group", color = "white",
                 palette = c("#DF7861","#8CC0DE") )

pdf("percent_monocyte.pdf",width =8,height = 8)
ggpubr::ggarrange(p1,p2,p3,p4,nrow = 2,ncol =1)
dev.off()


```

### 4.Trait associated genes for 10 to 1000

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/monocytecount_bmmc_Pagwas1.10.0.RData")
load("Pagwas_groundtruth.monocyte.RData")

#table(Pagwas$celltypes)
meta.data<-Pagwas@meta.data

oa<-order(Pagwas@misc$gene_heritability_correlation,decreasing = T)

topgene<-lapply(1:100, function(i){
  g1<-rownames(Pagwas@misc$gene_heritability_correlation)[oa[1:(10*i)]]
})
names(topgene)<-paste0(1:100,"times_topgene")
Pagwas_groundtruth <- Seurat::AddModuleScore(Pagwas_groundtruth, 
                                 assay = "RNA", 
                                 topgene,name="scPagwas")

i<-"monocytecount"

load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_magma_genes.RData"))

#load("Pagwas_groundtruth.monocyte.RData")

#Pagwas_groundtruth<-Pagwas[,colnames(Pagwas_groundtruth)]
Pagwas_groundtruth$types<-"monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="09_pDC"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="10_cDC"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="19_CD8.N"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="20_CD4.N1"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="25_NK"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="17_B"]<-"non_monocyte"

percent1<-lapply(paste0("scPagwas",1:100), function(i){
  a<-Pagwas_groundtruth@meta.data$types[order(Pagwas_groundtruth@meta.data[,i],decreasing = T)[1:5000]]
  return(sum(a=="monocyte")/5000)
})

#####scDRS

ab<-lapply(1:100, function(i){
  b<-paste(magma_genes$symbol[1:(10*i)],collapse=",")
  return(b)
})

b<-paste0("MAGMAtop",1:100)
a<-data.frame(genes=unlist(ab))
rownames(a)<-b
write.csv(a,file="magma_scRDS_10_1000.csv")

########——————————————————————————————————————————————
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("Pagwas_groundtruth.monocyte.RData")
library(SeuratDisk)
library(Seurat)
a<-GetAssayData(Pagwas_groundtruth,assay = "RNA",slot="counts")
Pagwas_readgroundtruth <- CreateSeuratObject(counts = a,
                           project = "scPagwas",
                           min.cells = 3,
                           min.features = 200)


DefaultAssay(Pagwas_readgroundtruth) <- "RNA"
SaveH5Seurat(Pagwas_readgroundtruth, "Pagwas_groundtruth_addata.h5seurat")
Convert("Pagwas_groundtruth_addata.h5seurat", dest="h5ad")
```

### scDRS for magma and scPagwas

```python

DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/"
H5AD_FILE = os.path.join(DATA_PATH, "Pagwas_groundtruth_addata.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
#sc.tl.pca(adata, svd_solver="arpack")
#sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
scdrs.preprocess(adata)

i='monocytecount'
df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_scRDS"+i+".csv", index_col=0)
df_gs = df_gs.loc[["scPagwastop1000","magmatop1000"],:]
df_gs=  df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".geneset.gs")

df_score = scdrs.score_cell(
             data=adata,
             gene_list=df_gs['magmatop1000'][0],
             gene_weight=df_gs['magmatop1000'][1],
             ctrl_match_key="mean_var",
             n_ctrl=200,
             weight_opt="vs",
             return_ctrl_raw_score=False,
             return_ctrl_norm_score=True,
             verbose=False)

df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/mono.bmmc.1000magma.df_model.csv", sep=",")

df_score = scdrs.score_cell(
             data=adata,
             gene_list=df_gs['scPagwastop1000'][0],
             gene_weight=df_gs['scPagwastop1000'][1],
             ctrl_match_key="mean_var",
             n_ctrl=200,
             weight_opt="vs",
             return_ctrl_raw_score=False,
             return_ctrl_norm_score=True,
             verbose=False)

df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/mono.bmmc.1000scPagwas.df_model.csv", sep=",")
```



### 5.R run plot for graduate nunber of genes

```R


Pagwas_groundtruth$types<-"monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="09_pDC"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="10_cDC"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="19_CD8.N"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="20_CD4.N1"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="25_NK"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="17_B"]<-"non_monocyte"


percent1<-lapply(paste0("Cluster",1:100), function(i){
  a<-Pagwas_groundtruth@meta.data$types[order(Pagwas_groundtruth@meta.data[,i],decreasing = T)[1:5000]]
  return(sum(a=="monocyte")/5000)
})

df<-read.csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/MAGMA_scDRSresult_bmmcmodel.csv")

 percent2<-lapply(paste0("MAGMAtop",1:100), function(i){
   a<-df_list$types[order(df[,i],decreasing = T)[1:5000]]
   return(sum(a=="monocyte")/5000)
 })
 
###smultixcan

ad<-read.csv("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/smultixcan/scDRS_result/bmmcmodel/monocytecount_scDRSresult.csv")
colnames(ad)<-paste0("smultixcantop",1:100)
#Pagwas_groundtruth@meta.data<-cbind(Pagwas_groundtruth@meta.data,ad)
percent3<-lapply(paste0("smultixcantop",1:100), function(i){
  a<-Pagwas_groundtruth@meta.data$types[order(ad[,i],decreasing = T)[1:5000]]
  return(sum(a=="monocyte")/5000)
})
#spredixcan
ad<-read.csv("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/spredixcan/scDRS_result/bmmcmodel/monocytecount_scDRSresult.csv")
colnames(ad)<-paste0("spredixcantop",1:100)
#Pagwas_groundtruth@meta.data<-cbind(Pagwas_groundtruth@meta.data,ad)
percent4<-lapply(paste0("spredixcantop",1:100), function(i){
  a<-Pagwas_groundtruth@meta.data$types[order(ad[,i],decreasing = T)[1:5000]]
  return(sum(a=="monocyte")/5000)
})
#TWAS
ad<-read.csv("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/TWAS/scDRS_result/bmmcmodel/monocytecount_scDRSresult.csv")
colnames(ad)<-paste0("TWAStop",1:100)
#Pagwas_groundtruth@meta.data<-cbind(Pagwas_groundtruth@meta.data,ad)
percent5<-lapply(paste0("TWAStop",1:100), function(i){
  a<-Pagwas_groundtruth@meta.data$types[order(ad[,i],decreasing = T)[1:5000]]
  return(sum(a=="monocyte")/5000)
})
###########
#可视化

 library(ggplot2)
 library(ggpubr)
 gg1<-data.frame(Number=c(1:100 *10,1:100 *10,1:100 *10,1:100 *10,1:100 *10),
                 Percent=c(unlist(percent1),unlist(percent2),unlist(percent3),unlist(percent4),unlist(percent5)),

                 Types=c(rep("scPagwas",100),rep("MAGMA",100),rep("smultixcan",100),rep("spredixcan",100),rep("TWAS",100))
 )
gg1$Types<-factor(gg1$Types,levels=c("scPagwas","MAGMA","TWAS","spredixcan","smultixcan"))
 setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
 pdf("Monocyte_bmmc_real_percent.pdf")
 ggplot(gg1,aes(Number,Percent,group=Types, color=Types))+
   geom_point()+
   geom_line(position = position_dodge(0.1),cex=1)+
scale_color_manual(values=c('#e8505b','#f9d56e','#f3ecc2','#14b1ab','#d4b5b0'))+
scale_fill_manual(values=c('#e8505b','#f9d56e','#f3ecc2','#14b1ab','#d4b5b0'))+
   labs(x="Number of genes", y="Percent of Monocyte cells in top half cells") +
   theme_classic()
 dev.off()


```

### 6.RANK PLOT

```R
library('ComplexHeatmap')
library(circlize)
####
colors_celltypes=c("#DF7861","#8CC0DE")
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("Pagwas_groundtruth.monocyte.RData")
load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/monocytecount_bmmc_Pagwas1.10.0.RData")
Pagwas_groundtruth$scPagwas.TRS.Score1<-Pagwas@meta.data[colnames(Pagwas_groundtruth),"scPagwas.TRS.Score1"]
scDRS_re<-read.csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/ThreeMethods_scDRSresult_bmmcmodel.csv")

Pagwas_groundtruth@meta.data$smultixcan_top100_monocytecount<-scDRS_re$smultixcan_top100_monocytecount
Pagwas_groundtruth@meta.data$spredixcan_top100_monocytecount<-scDRS_re$spredixcan_top100_monocytecount
Pagwas_groundtruth@meta.data$TWAS_top100_monocytecount<-scDRS_re$TWAS_top100_monocytecount

Pagwas_groundtruth$types<-"monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="09_pDC"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="10_cDC"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="19_CD8.N"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="20_CD4.N1"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="25_NK"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="17_B"]<-"non_monocyte"
all_fortify_can<-fortify.Seurat.umap(Pagwas_groundtruth)

all_fortify_can<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_scDRS,decreasing=T),]

pdf("BMMCmodel_monocyte_scPagwas_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "scPagwas rank"))
dev.off()

all_fortify_can<-all_fortify_can[order(all_fortify_can$Cluster100,decreasing=T),]

pdf("BMMCmodel_monocyte_scPagwas_seurat_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "scPagwas seurat rank"))
dev.off()


all_fortify_can<-all_fortify_can[order(all_fortify_can$magma_monocytecount_scDRS,decreasing=T),]

pdf("BMMCmodel_monocyte_magma_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "MAGMA rank"))
dev.off()

###smultixcan
all_fortify_can<-all_fortify_can[order(all_fortify_can$smultixcan_top100_monocytecount,decreasing=T),]
pdf("BMMCmodel_monocyte_Smultixcan_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "Smultixcan rank"))
dev.off()

#spredixcan

all_fortify_can<-all_fortify_can[order(all_fortify_can$spredixcan_top100_monocytecount,decreasing=T),]
pdf("BMMCmodel_monocyte_Spredixcan_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "Spredixcan rank"))
dev.off()


#TWAS
all_fortify_can<-all_fortify_can[order(all_fortify_can$TWAS_top100_monocytecount,decreasing=T),]
pdf("BMMCmodel_monocyte_TWAS_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "TWAS rank"))
dev.off()


all_fortify_can<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_AUCell,decreasing=T),]

pdf("BMMCmodel_monocyte_scPagwas_AUCell_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "scPagwas rank"))
dev.off()

all_fortify_can<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_Vision,decreasing=T),]

pdf("BMMCmodel_monocyte_scPagwas_Vision_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "scPagwas rank"))
dev.off()

save(Pagwas_groundtruth,file='Pagwas_groundtruth_bmmc_monocyte_v10.RData')
```

### 7.percent plot

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load('Pagwas_groundtruth_bmmc_monocyte_v10.RData')

all_fortify_can<-fortify.Seurat.umap(Pagwas_groundtruth)
percent_plot<-function(df1){
        n<-nrow(df1)
        b<-rep("top",n)
        #b[1:(0.5*n)]<-1
        b[(0.5*n+1):n]<-"bottom"

        df1$rg<-b
        #df1$type
        print(table(data.frame(df1$type,df1$rg)))
        t1<-table(data.frame(df1$type,df1$rg))
        percent1<-t1[,2]/sum(t1[,2])
print(percent1)
        df <- data.frame(
          group =c("monocyte","non_monocyte"),
          value = percent1)

        p1<-ggdonutchart(df, "value", label = "group",
                         fill = "group", color = "white",
                         palette = c("#DF7861", "#D4E2D4") )
        return(p1)
}


df1<-all_fortify_can[order(all_fortify_can$Cluster100,decreasing=T),]
p1<-percent_plot(df1)
#    monocyte non_monocyte
#       0.982        0.018


df1<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_scDRS,decreasing=T),]
p1<-percent_plot(df1)
# monocyte non_monocyte
#      0.9836       0.0164

df1<-all_fortify_can[order(all_fortify_can$magma_monocytecount_scDRS,decreasing=T),]
p2<-percent_plot(df1)
#    monocyte non_monocyte
#      0.7874       0.2126

df1<-all_fortify_can[order(all_fortify_can$smultixcan_top100_monocytecount,decreasing=T),]
p3<-percent_plot(df1)
#     monocyte non_monocyte
#      0.6704       0.3296

df1<-all_fortify_can[order(all_fortify_can$spredixcan_top100_monocytecount,decreasing=T),]
p4<-percent_plot(df1)
#monocyte non_monocyte
#      0.6174       0.3826
df1<-all_fortify_can[order(all_fortify_can$TWAS_top100_monocytecount,decreasing=T),]
p5<-percent_plot(df1)
#  monocyte non_monocyte
 #     0.7104       0.2896
pdf("percent.methods_monocyte_realBMMC_scdrs.pdf",width = 10,height = 8)
ggpubr::ggarrange(p1,p2,p3,p4,p5,nrow = 2)
dev.off()
```

### 8.Compare with other score methods

2022/11/25日重新计算magma AUCell结果检查

```
files<-"/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/mono.1000magma.df_modeldata1.csv"

magma_genes<-read.table("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/monocytecountmagmatop1000.txt",header=T)

a<-c("MAGMA1000","NA",magma_genes$x)
a<-matrix(a,nrow=1)
out<-"/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/monocytecountmagmatop1000.genes.mgt"
write.table(a,file=out,append = T,col.names = F,row.names = F,quote = F,sep = "\t")


  auccell = cal_PAS(seurat_object = Pagwas_groundtruth,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/monocytecountmagmatop1000.genes.mgt")

  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
  Pagwas_groundtruth$magma_monocytecount_AUCell<-auccell_df$MAGMA1000 

df1<-all_fortify_can[order(all_fortify_can$magma_monocytecount_AUCell,decreasing=T),]
p<-percent_plot(df1)

save(Pagwas_groundtruth,file="Pagwas_groundtruth_bmmc_monocyte_v10.RData")
```



```R
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
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")

load('Pagwas_groundtruth_bmmc_monocyte_v10.RData')
##输出mgt文件

gwass<-c('monocytecount')
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

counts = load_counts()
se_oj = CreateSeuratObject(counts)
se_oj = cal_PAS(seurat_object = se_oj,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = 'log',
              species = 'mouse', 
              pathway='kegg')

DefaultAssay(object = Pagwas) <- "RNA"

methods<-c('smultixcan', 'spredixcan' ,'TWAS')
gwas='monocytecount'
m_l<-lapply(methods,function(method){
  gmt<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
  DefaultAssay(object = Pagwas_groundtruth) <- "RNA"
  auccell = cal_PAS(seurat_object = Pagwas_groundtruth,
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
Pagwas_groundtruth@meta.data<-cbind(Pagwas_groundtruth@meta.data,auc_df)

###scPagwas magma
auccell = cal_PAS(seurat_object = Pagwas_groundtruth,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
gmt_file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scPagwas.magma.Model.monocytecount.gmt")

auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
Pagwas_groundtruth$magma_monocytecount_AUCell<-auccell_df$magmatop1000 
Pagwas_groundtruth$scPagwas_monocytecount_AUCell<-auccell_df$scPagwastop1000

####################
######Vision
##################
m_l<-lapply(methods,function(method){
  gmt<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
  DefaultAssay(object = Pagwas_groundtruth) <- "RNA"

  auccell = cal_PAS(seurat_object = Pagwas_groundtruth,
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
Pagwas_groundtruth@meta.data<-cbind(Pagwas_groundtruth@meta.data,df)

###scPagwas magma
auccell = cal_PAS(seurat_object = Pagwas_groundtruth,
                       tool = 'Vision',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
gmt_file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scPagwas.magma.Model.monocytecount.gmt")

auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
Pagwas_groundtruth$magma_monocytecount_Vision<-auccell_df$magmatop1000 
Pagwas_groundtruth$scPagwas_monocytecount_Vision<-auccell_df$scPagwastop1000
save(Pagwas_groundtruth,file="Pagwas_groundtruth_bmmc_monocyte_v10.RData")

#1.umap plot
all_fortify_can <- fortify.Seurat.umap(Pagwas_groundtruth)

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas_monocytecount_AUCell), size = 0.3, alpha = 1) +
        umap_theme() +
         scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("scPagwas + AUCell")

pdf(file="./Umap.monocytecount_BMMCgroundtruth.scPagwas.AUCell.pdf",width = 6, height = 6)
print(p1)
dev.off()

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas_monocytecount_Vision), size = 0.3, alpha = 1) +
        umap_theme() +
         scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("scPagwas + Vision")

pdf(file="./Umap.monocytecount_BMMCgroundtruth.scPagwas.Vision.pdf",width = 6, height = 6)
print(p1)
dev.off()

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =magma_monocytecount_AUCell), size = 0.3, alpha = 1) +
        umap_theme() +
         scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("MAGMA + AUCell")

pdf(file="./Umap.monocytecount_BMMCgroundtruth.MAGMA.AUCell.pdf",width = 8, height = 8)
print(p1)
dev.off()

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =magma_monocytecount_Vision), size = 0.3, alpha = 1) +
        umap_theme() +
         scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("Vision + AUCell")

pdf(file="./Umap.monocytecount_BMMCgroundtruth.MAGMA.Vision.pdf",width = 8, height = 8)
print(p1)
dev.off()

df1<-all_fortify_can[order(all_fortify_can$Cluster100,decreasing=T),]
p<-percent_plot(df1)
#  monocyte non_monocyte
#       0.982        0.018
df1<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_scDRS,decreasing=T),]
p<-percent_plot(df1)
# monocyte non_monocyte
#      0.9836       0.0164
df1<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_AUCell,decreasing=T),]
p<-percent_plot(df1)
#monocyte non_monocyte
#      0.9836       0.0164
df1<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_Vision,decreasing=T),]
p<-percent_plot(df1)
# monocyte non_monocyte
#       0.983        0.017
df1<-all_fortify_can[order(all_fortify_can$magma_monocytecount_Vision,decreasing=T),]
p<-percent_plot(df1)
#   monocyte non_monocyte
#      0.9124       0.0876
df1<-all_fortify_can[order(all_fortify_can$magma_monocytecount_AUCell,decreasing=T),]
p<-percent_plot(df1)
#  monocyte non_monocyte
#      0.9004       0.0996

save(Pagwas_groundtruth,file="Pagwas_groundtruth_bmmc_monocyte_v10.RData")

colors_celltypes=c("#DF7861","#8CC0DE")
########绘制umap图
library('ComplexHeatmap')
all_fortify_can<-all_fortify_can[order(all_fortify_can$magma_monocytecount_Vision,decreasing=T),]
#library()
pdf("BMMCmodel_monocyte_magma_Vision_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "MAGMA Vision"))
dev.off()

all_fortify_can<-all_fortify_can[order(all_fortify_can$magma_monocytecount_AUCell,decreasing=T),]
#library()
pdf("BMMCmodel_monocyte_magma_AUCell_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "MAGMA AUCell"))
dev.off()
```



## 二、PBMC

### 1.Run Pagwas and get top1000 genes

```R
library(scPagwas)
library(ggplot2)
suppressMessages(library(Seurat))
suppressMessages(library("dplyr"))
i<-"monocytecount"
Pagwas_pbmc_groundtruth<-scPagwas_main(Pagwas = NULL,
                                       gwas_data=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"),
                                       Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_pbmc_realgroundtruth.rds",
                                       output.prefix="real_pbmc_monocytecount",
                                       output.dirs="real_pbmc_monocytecount",
                                       block_annotation = block_annotation,
                                       assay="RNA",
                                       Pathway_list=Genes_by_pathway_kegg,
                                       chrom_ld = chrom_ld,
                                       singlecell=T,
                                       seruat_return=T,
                                       celltype=F,
                                       ncores = 1)
save(Pagwas_pbmc_groundtruth,file = "/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_pbmc_groundtruth_monocytecount.RData")

load("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_pbmc_groundtruth_monocytecount.RData")
i<-"monocytecount"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_magma_genes.RData"))

scPagwas_genes<-names(Pagwas_pbmc_groundtruth@misc$gene_heritability_correlation[order(Pagwas_pbmc_groundtruth@misc$gene_heritability_correlation,decreasing=T),])

magmatop1000<-intersect(magma_genes$symbol,rownames(Pagwas_pbmc_groundtruth))[1:1000]
scPagwastop1000<-scPagwas_genes[1:1000]
magmatop500<-intersect(magma_genes$symbol,rownames(Pagwas_pbmc_groundtruth))[1:500]
scPagwastop500<-scPagwas_genes[1:500]
topgene<-list(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500)
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")

a<-data.frame(genes=c(paste(scPagwastop1000,collapse=","),paste(magmatop1000,collapse=","),paste(scPagwastop500,collapse=","),paste(magmatop500,collapse=",")))
rownames(a)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
write.csv(a,file=paste0("scPagwas.magma.PBMC.topgenes.",i,".csv"))

```



### 3.Get top genes for 10 to 1000

```R
i<-"monocytecount"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_magma_genes.RData"))

ab<-lapply(1:100, function(i){
  b<-paste(magma_genes$symbol[1:(10*i)],collapse=",")
  return(b)
})
b<-paste0("MAGMAtop",1:100)
a<-data.frame(genes=unlist(ab))
rownames(a)<-b
write.csv(a,file="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/PBMC_monocyte_magma_scRDS_10_1000.csv")


setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("Pagwas_pbmc_groundtruth_monocytecount.RData")
library(SeuratDisk)
library(Seurat)
a<-GetAssayData(Pagwas_pbmc_groundtruth,assay = "RNA",slot="counts")
Pagwas_pbmc_groundtruth <- CreateSeuratObject(counts = a,
                                             project = "scPagwas",
                                             min.cells = 3,
                                             min.features = 200)
DefaultAssay(Pagwas_pbmc_groundtruth) <- "RNA"
SaveH5Seurat(Pagwas_pbmc_groundtruth, "Pagwas_pbmc_groundtruth.h5seurat")
Convert("Pagwas_pbmc_groundtruth.h5seurat", dest="h5ad") 

```

### 4.Run scDRS for top 10-1000 genes

```python
##scDRS
# load adata
DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/"
H5AD_FILE = os.path.join(DATA_PATH, "Pagwas_pbmc_groundtruth.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
scdrs.preprocess(adata)

df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/PBMC_monocyte_magma_scRDS_10_1000.csv", index_col=0)
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/geneset.gs")

for i in range(1,101):
 a = 'MAGMAtop' + str(i)
 gene_list  = df_gs[a][0]
 gene_weight  = df_gs[a][1]
 df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=1000)
 df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/PBMC/"+a+".scDRS.csv", sep=",", index=False)
#######r
```

### 5.R run plot for graduate nunber of genes

```R
library(scPagwas)
library(ggplot2)
suppressMessages(library(Seurat))
suppressMessages(library("dplyr"))
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load('Pagwas_pbmc_groundtruth_monocytecount.RData')
oa<-order(Pagwas_pbmc_groundtruth@misc$gene_heritability_correlation,decreasing = T)
 topgene<-lapply(1:100, function(i){
   g1<-rownames(Pagwas_pbmc_groundtruth@misc$gene_heritability_correlation)[oa[1:(10*i)]]
 })
 names(topgene)<-paste0(1:100,"times_topgene")
 Pagwas_pbmc_groundtruth <- Seurat::AddModuleScore(Pagwas_pbmc_groundtruth, 
                                  assay = "RNA", 
                                  topgene,
                                  name="scPagwas")
 
 Pagwas_pbmc_groundtruth$types<-"non_monocyte"
 Pagwas_pbmc_groundtruth$types[Pagwas_pbmc_groundtruth$initial_clustering=="CD14"]<-"monocyte"
 
 percent1<-lapply(paste0("scPagwas",1:100), function(i){
   a<-Pagwas_pbmc_groundtruth@meta.data$types[order(Pagwas_pbmc_groundtruth@meta.data[,i],decreasing = T)[1:5000]]
   return(sum(a=="monocyte")/5000)
 })
 
 
df<-read.csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/MAGMA_scDRSresult_pbmcmodel.csv")
 percent2<-lapply(paste0("MAGMAtop",1:100), function(i){
   a<-Pagwas_pbmc_groundtruth@meta.data$types[order(df[,i],decreasing = T)[1:5000]]
   return(sum(a=="monocyte")/5000)
 })
 
###smultixcan
setwd("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/smultixcan/scDRS_result/pbmcmodel")
ad<-read.csv("monocytecount_scDRSresult.csv")
colnames(ad)<-paste0("smultixcantop",1:100)
#Pagwas_pbmc_groundtruth@meta.data<-cbind(Pagwas_pbmc_groundtruth@meta.data,ad)
percent3<-lapply(paste0("smultixcantop",1:100), function(i){
  a<-Pagwas_pbmc_groundtruth@meta.data$types[order(ad[,i],decreasing = T)[1:5000]]
  return(sum(a=="monocyte")/5000)
})
#spredixcan
setwd("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/spredixcan/scDRS_result/pbmcmodel")
ad<-read.csv("monocytecount_scDRSresult.csv")
colnames(ad)<-paste0("spredixcantop",1:100)
#Pagwas_pbmc_groundtruth@meta.data<-cbind(Pagwas_pbmc_groundtruth@meta.data,ad)
percent4<-lapply(paste0("spredixcantop",1:100), function(i){
  a<-Pagwas_pbmc_groundtruth@meta.data$types[order(ad[,i],decreasing = T)[1:5000]]
  return(sum(a=="monocyte")/5000)
})
#TWAS
setwd("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/TWAS/scDRS_result/pbmcmodel")
ad<-read.csv("monocytecount_scDRSresult.csv")
colnames(ad)<-paste0("TWAStop",1:100)
#Pagwas_pbmc_groundtruth@meta.data<-cbind(Pagwas_pbmc_groundtruth@meta.data,ad)
percent5<-lapply(paste0("TWAStop",1:100), function(i){
  a<-Pagwas_pbmc_groundtruth@meta.data$types[order(ad[,i],decreasing = T)[1:5000]]
  return(sum(a=="monocyte")/5000)
}) 
 
 library(ggplot2)
 library(ggpubr)
 gg1<-data.frame(Number=c(1:100 *10,1:100 *10,1:100 *10,1:100 *10,1:100 *10),
                 Percent=c(unlist(percent1),unlist(percent2),unlist(percent3),unlist(percent4),unlist(percent5)),

                 Types=c(rep("scPagwas",100),rep("MAGMA",100),rep("smultixcan",100),rep("spredixcan",100),rep("TWAS",100))
 )
gg1$Types<-factor(gg1$Types,levels=c("scPagwas","MAGMA","TWAS","spredixcan","smultixcan"))
 setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
 pdf("monocyte_realPBMC_percent.pdf")
 ggplot(gg1,aes(Number,Percent,group=Types, color=Types))+
   geom_point()+
   geom_line(position = position_dodge(0.1),cex=1)+
scale_color_manual(values=c('#e8505b','#f9d56e','#f3ecc2','#14b1ab','#d4b5b0'))+
scale_fill_manual(values=c('#e8505b','#f9d56e','#f3ecc2','#14b1ab','#d4b5b0'))+
   labs(x="Number of genes", y="Percent of monocyte cells in top half cells") +
   theme_classic()
 dev.off()

```

### 6.Rankplot

```
library('ComplexHeatmap')
library(circlize)
####
colors_celltypes=c("#DF7861","#8CC0DE")
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("Pagwas_pbmc_groundtruth_monocytecount.RData")

Pagwas_pbmc_groundtruth$types<-"non_monocyte"
Pagwas_pbmc_groundtruth$types[Pagwas_pbmc_groundtruth$initial_clustering=="CD14"]<-"monocyte"
 
all_fortify_can<-fortify.Seurat.umap(Pagwas_pbmc_groundtruth)
all_fortify_can<-all_fortify_can[order(all_fortify_can$scPagwas.TRS.Score1,decreasing=T),]

pdf("PBMCmodel_monocyte_scPagwas_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "scPagwas rank"))
dev.off()


df<-read.csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/PBMC/MAGMAtop100.scDRS.csv")
Pagwas_pbmc_groundtruth$MAGMA<-df$zscore

all_fortify_can<-fortify.Seurat.umap(Pagwas_pbmc_groundtruth)
all_fortify_can<-all_fortify_can[order(all_fortify_can$MAGMA,decreasing=T),]

pdf("PBMCmodel_monocyte_magma_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "MAGMA rank"))
dev.off()

###smultixcan

ad<-read.csv("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/smultixcan/scDRS_result/pbmcmodel/monocytecount_scDRSresult.csv")
colnames(ad)<-paste0("smultixcantop",1:100)
Pagwas_pbmc_groundtruth@meta.data$Smultixcan<-ad$smultixcantop100
all_fortify_can<-fortify.Seurat.umap(Pagwas_pbmc_groundtruth)
all_fortify_can<-all_fortify_can[order(all_fortify_can$Smultixcan,decreasing=T),]
pdf("PBMCmodel_monocyte_Smultixcan_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "Smultixcan rank"))
dev.off()

#spredixcan
ad<-read.csv("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/spredixcan/scDRS_result/pbmcmodel/monocytecount_scDRSresult.csv")
colnames(ad)<-paste0("spredixcantop",1:100)
Pagwas_pbmc_groundtruth@meta.data$Spredixcan<-ad$spredixcantop100
all_fortify_can<-fortify.Seurat.umap(Pagwas_pbmc_groundtruth)
all_fortify_can<-all_fortify_can[order(all_fortify_can$Spredixcan,decreasing=T),]
pdf("PBMCmodel_monocyte_Spredixcan_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "Spredixcan rank"))
dev.off()


#TWAS
ad<-read.csv("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/TWAS/scDRS_result/pbmcmodel/monocytecount_scDRSresult.csv")
colnames(ad)<-paste0("TWAStop",1:100)
Pagwas_pbmc_groundtruth@meta.data$TWAS<-ad$TWAStop100

all_fortify_can<-fortify.Seurat.umap(Pagwas_pbmc_groundtruth)
all_fortify_can<-all_fortify_can[order(all_fortify_can$TWAS,decreasing=T),]
pdf("PBMCmodel_monocyte_TWAS_rankplot.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "TWAS rank"))
dev.off()
save(Pagwas_pbmc_groundtruth,file='Pagwas_pbmc_groundtruth_monocytecount.RData')
```

## subfunction

```
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


```

