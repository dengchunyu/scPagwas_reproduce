
#######--------------------------------------------------------------

#Add_Module_score for effector marker genes
library(tidyverse)
library(RColorBrewer)
library(Seurat)
library(SeuratData)

#FeaturePlot
#features <- c("CSNK2B")
#FeaturePlot(liverCells, features = features)

##The effector marker genes for calculating molecular signature scores for Navie CD8+T cells
effector_marker_genes<-c("PRDM1","PRF1","GZMB", "GNLY", "GZMA","IFNG", "FASLG")

#RidgePlot(liverCells, features = unique(pbc_risk_genes), ncol = 2)
#DotPlot(liverCells, features = unique(pbc_risk_genes)) + RotatedAxis()
#DotPlot(pbmc, features = unique(pbc_risk_genes)) + RotatedAxis()

#?Dotplot()
#?DiscretePalette()

scPagwas_naive_CD8T <- AddModuleScore(scPagwas_naive_CD8T,
                                      features = list(effector_marker_genes),
                                      name="effector_gene_scores")
W2<-FeaturePlot(scPagwas_naive_CD8T, reduction = "tsne",
                features = "effector_gene_scores1", label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



#Plotting the dot figures for the association between TRS and effector gene scores
library(ggplot2)
library(Seurat)
library(dplyr)

EGS <- as.matrix(scPagwas_naive_CD8T$effector_gene_scores1)
TRS <-as.matrix(scPagwas_naive_CD8T$scPagwas_TRS)
data_single_cells <- cbind(EGS,TRS)
colnames(data_single_cells) <- c("effector_gene_scores1","scPagwas_TRS")
data_single_cells <- as.data.frame(data_single_cells)

#取bins
binning <- data_single_cells %>% mutate(rank=ntile(data_single_cells$effector_gene_scores1,5))

xx <- aggregate(binning,by = list(binning$rank), mean)




##ggplot
ggplot(binning,aes(x=as.factor(rank),y=scPagwas_TRS))+
     geom_boxplot()+
     scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","#999999", "#E69F00"))+
    theme(legend.position="right")+
   labs(title="Plot of length  per dose",x="Dose (mg)", y = "Length")+
   geom_dotplot(binaxis='y', stackdir='center', stackratio=1.0, dotsize=0.5)

#######--------------------------------------------------------------
