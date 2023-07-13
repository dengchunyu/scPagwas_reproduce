##figure3a
library(SeuratObject)
library(Seurat)
sim_data=readRDS("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds")
Idents(sim_data)<-sim_data$celltype
library(ggpubr)
library(ggplot2)
color_seruat_patient <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")
pdf("model_groundtruth_monocyte.pdf",width = 6,height = 6)
 DimPlot(sim_data,reduction ="umap",
                         group.by = "celltype",pt.size=0.3,
                         label = F, repel=TRUE)+ 
  umap_theme()+
  scale_colour_manual(name = "celltype", values =color_seruat_patient[c(2:9,1)]) +
  theme(aspect.ratio=1)
dev.off()



##figure3b
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

pdf(file="Umap.monocytecount_model_scPagwas.TRS.Score.pdf",width = 6, height = 6)
print(p1)
dev.off()

#figure3d
p2<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =scDRS.1000magma), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scDRS.1000magma")

pdf(file="Umap.monocytecount_model_scDRS.1000magma.pdf",width = 6, height = 6)
print(p2)
dev.off()
#figure3c
p3<- ggplot() +
  geom_point(data = all_fortify_can,
             aes(x = UMAP_1, y = UMAP_2,color =scDRS.1000scPagwas), size = 0.3, alpha = 1) +
  umap_theme() +
  scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
  theme(aspect.ratio=1) +
  theme(legend.text=element_markdown(size=14),
        legend.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))  +
  ggtitle("scDRS.1000scPagwas")

pdf(file="Umap.monocytecount_model_scDRS.1000scPagwas.pdf",width = 6, height = 6)
print(p3)
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
#    monocyte non_monocyte
#       0.959        0.041

#magma
df1<-all_fortify_can[order(all_fortify_can$scDRS.1000magma,decreasing=T),]
#a<-scPagwas_magma_genelist[[i]]
n<-nrow(df1)
b<-rep(1,n)
b[1:(0.5*n)]<-1
b[(0.5*n+1):(1*n)]<-2

df1$rg<-b
print(table(data.frame(df1$types,df1$rg)))
t1<-table(data.frame(df1$types,df1$rg))
percent1<-t1[,1]/sum(t1[,1])
# monocyte non_monocyte
#        0.94         0.06
a <- data.frame(
  group = names(percent1),
  value = percent1)

################################################################
#figures13
###################
library(ggplot2)
library(ggpubr)
library(ggtext)
#################
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_output")
oa<-order(Pagwas@misc$gene_heritability_correlation,decreasing = T)
topgene<-lapply(1:100, function(i){
  g1<-rownames(Pagwas@misc$gene_heritability_correlation)[oa[1:(10*i)]]
})
names(topgene)<-paste0(1:100,"times_topgene")
Pagwas <- Seurat::AddModuleScore(Pagwas, 
                                 assay = "RNA", 
                                 topgene,
                                name="scPagwas")

Pagwas$types<-"non_monocyte"
Pagwas$types[Pagwas$celltype=="monocytes"]<-"monocyte"

percent<-lapply(paste0("scPagwas",1:100), function(i){
  a<-Pagwas@meta.data$types[order(Pagwas@meta.data[,i],decreasing = T)[1:1000]]
  return(sum(a=="monocyte")/1000)
})
#################
df_list<-lapply(1:100,function(i){
   df<-read.csv(paste0("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scDRS/MAGMAtop",i,".scDRS.csv"))
   return(df$zscore)
 })
 df_list<-as.data.frame(df_list)
colnames(df_list)<-paste0("MAGMAtop",1:100)
save(df_list,file = "/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scDRS/magma_df_list.RData")

df_list$types<-Pagwas$types
percent2<-lapply(paste0("MAGMAtop",1:100), function(i){
  a<-df_list$types[order(df_list[,i],decreasing = T)[1:1000]]
  return(sum(a=="monocyte")/1000)
})

library(ggplot2)
library(ggpubr)
gg1<-data.frame(Number=c(1:100 *10,1:100 *10),
           Percent=c(unlist(percent),unlist(percent2)),
           Types=c(rep("scPagwas",100),rep("MAGMA",100))
           )

pdf("Monocyte_model_percent.pdf")
ggplot(gg1,aes(Number,Percent,group=Types, color=Types))+
  geom_point()+
  geom_line(position = position_dodge(0.1),cex=1)+
  labs(x="Number of genes", y="Percent of Lymphocyte cells in top half cells") +
  theme_classic()
dev.off()
################

df_list<-lapply(c("1000scPagwas","1000magma"),function(i){
 df<-read.csv(paste0("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Mono",i,".df_model2.csv"))
 return(df$norm_score)
})

df_list<-as.data.frame(df_list)
Pagwas$scDRS.1000scPagwas<-df_list[[1]]
Pagwas$scDRS.1000magma<-df_list[[2]]
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")
all_fortify_can <- fortify.Seurat.umap(Pagwas)
all_fortify_can$types<-"non_monocyte"
all_fortify_can$types[all_fortify_can$celltype=="monocytes"]<-"monocyte"


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

#########################
#FigureS8
######################

setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata.RData")
library(scPagwas)
suppressMessages(library(Seurat))
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggtext)

df<-read.csv('/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/mono.1000scPagwas.df_modeldata1.csv')
Pagwas@meta.data$scPagwas_monocytecount_scDRS<-df$norm_score
df<-read.csv('/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/mono.1000magma.df_modeldata1.csv')
Pagwas@meta.data$magma_monocytecount_scDRS<-df$raw_score

scDRS_re<-read.csv('/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0/ThreeMethods_scDRSresult_model1.csv')
Pagwas@meta.data$smultixcan_scDRS<-scDRS_re$smultixcan_top100_monocytecount
Pagwas@meta.data$spredixcan_scDRS<-scDRS_re$spredixcan_top100_monocytecount
Pagwas@meta.data$TWAS_scDRS<-scDRS_re$TWAS_top100_monocytecount
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata.RData")
```

### 3.UMAP plot
FigureS8a
```R
###############umap
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
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata.RData")

all_fortify_can <- fortify.Seurat.umap(Pagwas)

#'smultixcan',
    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =smultixcan_scDRS), size = 0.3, alpha = 1) +
        umap_theme() +
         scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("Smultixcan + scDRS")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0/Umap.monocytecount_model1.smultixcan_scDRS_score.pdf",width = 6, height = 6)
print(p1)
dev.off()

 #'spredixcan' 
    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =spredixcan_scDRS), size = 0.3, alpha = 1) +
        umap_theme() +
         scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("Spredixcan + scDRS")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0/Umap.monocytecount_model1.spredixcan_scDRS_score.pdf",width = 6, height = 6)
print(p1)
dev.off()
#,'TWAS

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =TWAS_scDRS), size = 0.3, alpha = 1) +
        umap_theme() +
         scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("TWAS + scDRS")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0/Umap.monocytecount_model1.TWAS_scDRS_score.pdf",width = 6, height = 6)
print(p1)
dev.off()

save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata.RData")

save(all_fortify_can,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0/all_fortify_can.RData")

######################
#rank plot for figure3 and S8
###################

library('ComplexHeatmap')
library(circlize)
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_output")
####
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata.RData")
colors_celltypes=c("#DF7861","#8CC0DE")

all_fortify_can <- fortify.Seurat.umap(Pagwas)

all_fortify_can<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_scDRS,decreasing=T),]
#all_fortify_can$celltype<-as.vector(all_fortify_can$type)
pdf("rankplot_scPgawas.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$type),
              col =rev(colors_celltypes),
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "magma rank"))
dev.off()


all_fortify_can<-all_fortify_can[order(all_fortify_can$scDRS.1000magma,decreasing=T),]
#all_fortify_can$celltype<-as.vector(all_fortify_can$type)
pdf("rankplot_magma.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$type),
              col =rev(colors_celltypes),
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "magma rank"))
dev.off()


all_fortify_can<-all_fortify_can[order(all_fortify_can$smultixcan_scDRS,decreasing=T),]
#all_fortify_can$celltype<-as.vector(all_fortify_can$type)
pdf("rankplot_smultixcan.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$type),
              col =rev(colors_celltypes),
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "Smultixcan rank"))
dev.off()

all_fortify_can<-all_fortify_can[order(all_fortify_can$spredixcan_scDRS,decreasing=T),]
all_fortify_can$type<-as.vector(all_fortify_can$type)
pdf("rankplot_spredixcan.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$type),
              col =rev(colors_celltypes),
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "Spredixcan rank"))
dev.off()
#TWAS_scDRS
all_fortify_can<-all_fortify_can[order(all_fortify_can$TWAS_scDRS,decreasing=T),]
all_fortify_can$type<-as.vector(all_fortify_can$type)
pdf("rankplot_TWAS.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$type),
              col =rev(colors_celltypes),
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "TWAS rank"))
dev.off()
# Percent
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
          group =c("Non_Monocyte","Monocyte"),
          value = percent1)
}

df1<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_scDRS,decreasing=T),]
p1<-percent_plot(df1)
#       0     1
#0.053 0.947
df1<-all_fortify_can[order(all_fortify_can$magma_monocytecount_scDRS,decreasing=T),]
p2<-percent_plot(df1)
#     0    1
#0.06 0.94
df1<-all_fortify_can[order(all_fortify_can$smultixcan_scDRS,decreasing=T),]
p3<-percent_plot(df1)
#   0   1
#0.1 0.9
df1<-all_fortify_can[order(all_fortify_can$spredixcan_scDRS,decreasing=T),]
p4<-percent_plot(df1)
#0     1
#0.161 0.839
df1<-all_fortify_can[order(all_fortify_can$TWAS_scDRS,decreasing=T),]
p5<-percent_plot(df1)
#     0     1
#0.371 0.629

#############
#FigureS8c
################

#umap plot
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_output")
all_fortify_can <- fortify.Seurat.umap(Pagwas)
save(all_fortify_can,file="all_fortify_can.RData")
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

pdf(file="./Umap.monocytecount_model1.scPagwas.AUCell.pdf",width = 6, height = 6)
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

pdf(file="./Umap.monocytecount_model1.scPagwas.Vision.pdf",width = 6, height = 6)
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

pdf(file="./Umap.monocytecount_model1.magma.AUCell.pdf",width =7, height =7)
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
        ggtitle("MAGMA + Vision")

pdf(file="./Umap.monocytecount_model1.magma.Vision.pdf",width =7, height =7)
print(p1)
dev.off()

####rank plot
library('ComplexHeatmap')
library(circlize)
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_output")
####
colors_celltypes=c("#DF7861","#8CC0DE")

all_fortify_can<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_Vision,decreasing=T),]
#all_fortify_can$celltype<-as.vector(all_fortify_can$type)
pdf("rankplot_Vision.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$type),
              col =rev(colors_celltypes),
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "Vision"))
dev.off()

all_fortify_can<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_AUCell,decreasing=T),]
#all_fortify_can$celltype<-as.vector(all_fortify_can$type)
pdf("rankplot_AUCell.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$type),
              col =rev(colors_celltypes),
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "AUCell"))
dev.off()

#scDRS.1000scPagwas
all_fortify_can<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_scDRS,decreasing=T),]
#all_fortify_can$celltype<-as.vector(all_fortify_can$type)
pdf("rankplot_scDRS.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$type),
              col =rev(colors_celltypes),
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "scDRS"))
dev.off()

#seruat

all_fortify_can<-all_fortify_can[order(all_fortify_can$scPagwas.TRS.Score1,decreasing=T),]
#all_fortify_can$celltype<-as.vector(all_fortify_can$type)
pdf("rankplot_seurat.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$type),
              col =rev(colors_celltypes),
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "Seurat"))
dev.off()

all_fortify_can<-all_fortify_can[order(all_fortify_can$magma_monocytecount_Vision,decreasing=T),]
#all_fortify_can$celltype<-as.vector(all_fortify_can$type)
pdf("magma_rankplot_Vision.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$type),
              col =rev(colors_celltypes),
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "magma Vision"))
dev.off()
all_fortify_can<-all_fortify_can[order(all_fortify_can$magma_monocytecount_AUCell,decreasing=T),]
#all_fortify_can$celltype<-as.vector(all_fortify_can$type)
pdf("magma_rankplot_AUCell.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$type),
              col =rev(colors_celltypes),
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "magma AUCell"))
dev.off()
###########
#percent plot
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
          group =c("Non_Monocyte","Monocyte"),
          value = percent1)
}


df1<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_scDRS,decreasing=T),]
p1<-percent_plot(df1)
#   0     1
#0.053 0.947

df1<-all_fortify_can[order(all_fortify_can$scPagwas.TRS.Score1,decreasing=T),]
p2<-percent_plot(df1)
#     0     1
#0.041 0.959
df1<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_Vision,decreasing=T),]
p3<-percent_plot(df1)
#     0     1
#0.041 0.959
df1<-all_fortify_can[order(all_fortify_can$scPagwas_monocytecount_AUCell,decreasing=T),]
p4<-percent_plot(df1)
#    0     1
#0.038 0.962

df1<-all_fortify_can[order(all_fortify_can$magma_monocytecount_Vision,decreasing=T),]
p1<-percent_plot(df1)
#   0    1
#0.05 0.95
df1<-all_fortify_can[order(all_fortify_can$magma_monocytecount_AUCell,decreasing=T),]
p1<-percent_plot(df1)
#     0     1
#0.039 0.961
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata.RData")


##############################
#Lymphocytecount FigureS9
#############################
### Umap
# FigureS9A
scDRS_re<-read.csv('/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount31000magma.df_model.csv')
Pagwas$magma_scDRS_Lymphocytecount3.score <- scDRS_re$norm_score

scDRS_re<-read.csv('/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount31000scPagwas.df_model.csv')
Pagwas$scPagwas_scDRS_Lymphocytecount3.score <- scDRS_re$norm_score

load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas.RData")

genes<-names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation, decreasing = T), ])[1:1000]

load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Pagwas_Lymphocytecount_modeldata.RData")

Pagwas <- Seurat::AddModuleScore(Pagwas, assay = "RNA", list(genes), name = "scPagwas_seurat_score")

save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Pagwas_Lymphocytecount_modeldata.RData")

#all_fortify_can<-all_fortify_can[,c(1:2,6,10,129:134,138:142)]

all_fortify_can <- fortify.Seurat.umap(Pagwas)

save(all_fortify_can,file="Lymphocytecount3_model1_all_fortify_can.RData")
   p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas_seurat_score1), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("scPagwas.TRS.seurat.Lymphocytecount3")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Umap.Lymphocytecount3_model1.scPagwas_seurat_score.pdf",width =6, height = 6)
print(p1)
dev.off()

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =magma_scDRS_Lymphocytecount3.score), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("magma1000_scDRS_score")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Umap.Lymphocytecount3_model1.magma1000_scDRS_score.pdf",width = 8, height = 6)
print(p1)
dev.off()

p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas_scDRS_Lymphocytecount3.score ), size = 0.3, alpha = 1) +
        umap_theme() +
      scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("scPagwas_scDRS_Lymphocytecount3.score")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Umap.Lymphocytecount3_model1.scPagwas1000_scDRS_score.pdf",width = 8, height = 6)
print(p1)
dev.off()


    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =smultixcan_scDRS_Lymphocytecount3), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("smultixcan_scDRS_Lymphocytecount3")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Umap.Lymphocytecount3_model1.smultixcan1000_scDRS_score.pdf",width = 8, height = 6)
print(p1)
dev.off()

   p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =spredixcan_scDRS_Lymphocytecount3), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("spredixcan_scDRS_Lymphocytecount3")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Umap.Lymphocytecount3_model1.spredixcan1000_scDRS_score.pdf",width = 8, height = 6)
print(p1)
dev.off()

   p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =TWAS_scDRS_Lymphocytecount3), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("TWAS_scDRS_Lymphocytecount3")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Umap.Lymphocytecount3_model1.TWAS1000_scDRS_score.pdf",width =8, height = 6)
print(p1)
dev.off()

## sub-functions
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



