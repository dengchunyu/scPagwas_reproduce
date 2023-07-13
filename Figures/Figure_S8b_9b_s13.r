### 2.TSNE plot for groundtruth example celltypes
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


### 2.1umap for other methods

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


### 3,Heatmap and ranked plot for scPagwas and magma
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




### 5.R run plot for graduate nunber of genes

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

### 6.RANK PLOT

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


### 7.percent plot

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

##################################
#PBMC real groundtruth

### 5.R run plot for graduate nunber of genes

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
