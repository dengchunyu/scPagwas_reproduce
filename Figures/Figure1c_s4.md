

# scPagwas data density and sparsity

## Density for Pathway PCA and single cell 

Supplementary Figure S4

```R
load("E:/OneDrive/GWAS_Multiomics/compare/Hema_test2/Lymphocytecount3_Hema_bmmc_scPagwas.RData")
pca_scCell_mat <- GetAssayData(Pagwas,assay = "scPagwasPaPca")
median<- apply(pca_scCell_mat,1,median)
gg_pca<-data.frame(median=median,pathway=rownames(pca_scCell_mat))
library(reshape2)

pdf("E:/OneDrive/GWAS_Multiomics/Distribution_test/Pathway_gene_distribution.pdf",width = 4,height = 6)
ggplot(gg_pca, aes(x=median))+
  geom_density()+
  #scale_color_brewer()+
  theme_classic()+
labs(title="pathway")+
 theme(legend.position="none")
dev.off()

gg_scgene<-as.data.frame(Pagwas$data_mat[,1:10])
gg_scgene$gene<-rownames(gg_scgene)
#library(reshape2)
gg_scgene<-reshape2::melt(gg_scgene,id.vars="gene")

pdf("E:/OneDrive/GWAS_Multiomics/Distribution_test/singlegene_distribution.pdf",width = 4,height = 6)
ggplot(gg_scgene, aes(x=value,color=variable))+
  geom_density()+
  scale_color_brewer()+
  theme_classic()+
labs(title="gene")+
 theme(legend.position="none")
dev.off()
```

## Variance for pathway and gene

**Figure 1C**

```R
library("ggplot2")
library(ggpubr)
library(scPagwas)
suppressMessages(library(Seurat))

setwd("/share/pub/dengcy/GWAS_Multiomics/Distribution_test")
load("/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/scPgwas/scpagwas_HCA_Migraine.RData")
Pagwas$Pathway_ld_gwas_data<-NULL
as_matrix <- function(mat) {
  tmp <- matrix(data = 0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  row_pos <- mat@i + 1
  col_pos <- findInterval(seq(mat@x) - 1, mat@p[-1]) + 1
  val <- mat@x
  for (i in seq_along(val)) {
    tmp[row_pos[i], col_pos[i]] <- val[i]
  }
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

##
data_mat <- Pagwas$data_mat
#data_mat <- as_matrix(data_mat)
data_mat2 <-Pagwas$pca_scCell_mat
paths<-rownames(data_mat2)

x_list<-lapply(paths,function(i){
    gene<- Genes_by_pathway_kegg[[i]]
    x2<-data_mat[intersect(gene,rownames(data_mat)),]
    x2<- as_matrix(x2)
     x2 <- apply(x2, 2, function(x) (x - min(x)) / (max(x) - min(x)))
    return(x2)
})
data_mat2<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),x_list)

x_list<-lapply(paths,function(i){
    gene<- Genes_by_pathway_kegg[[i]]
    x2<-data_mat[intersect(gene,rownames(data_mat)),]
    return(x2)
})
data_mat<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),x_list)
              
pa_var_mean_df<-data.frame(
 pathway=rownames(data_mat2),
 var_x_Pa=unlist(apply(data_mat2,1,var)),
 mean_x_Pa=unlist(apply(data_mat2,1,mean)))

 gene_sparsity<-apply(data.matrix(1:nrow(data_mat)),1,function(x){
  a<-data_mat[x,]
  return(sum(a==0)/length(a))
})
 gene_var_mean_df<-data.frame(  
     gene=rownames(data_mat),
     var_x_ge=unlist(apply(data.matrix(1:nrow(data_mat)),1,function(x){
  a<-data_mat[x,]
  return(var(a))
})),
 mean_x_ge=unlist(apply(data.matrix(1:nrow(data_mat)),1,function(x){
  a<-data_mat[x,]
  return(mean(a))
}))

pa_var_mean_df$mean_x_Pa<-  log10(abs(pa_var_mean_df$mean_x_Pa))
gene_var_mean_df$mean_x_ge<-  log10(abs(gene_var_mean_df$mean_x_ge))

pdf("HCA_smoothScatter_pathway.pdf")
smoothScatter(pa_var_mean_df$mean_x_Pa,pa_var_mean_df$var_x_Pa, 
              colramp = colorRampPalette(c("white", "#B20600")),
              transformation =  function(x) x^.2,
              col = "black",
              ylim=c(0,3),
              xlab = "snp-pathway log10(x)", 
              ylab = "Variance")
lines(lowess(pa_var_mean_df$mean_x_Pa,pa_var_mean_df$var_x_Pa),col='#B20600',lwd=1,lty=2)
dev.off()

pdf("HCA_smoothScatter_genes.pdf")
smoothScatter(gene_var_mean_df$mean_x_ge,gene_var_mean_df$var_x_ge,
              colramp = colorRampPalette(c("white", "#398AB9")),
              col = "black",
              transformation =  function(x) x^.4,
              ylim=c(0,3),
              xlab = "snp-genes log10(x)", 
              ylab = "Variance")
lines(lowess(gene_var_mean_df$mean_x_ge, gene_var_mean_df$var_x_ge),col='#243A73',lwd=1.5,lty=1)
dev.off()
```

## scRNA-seq and pathway data sparsity proportion plot
Supplementary Figure S4
```R
library("ggplot2")
library(ggpubr)
library(scPagwas)
suppressMessages(library(Seurat))

setwd("/share/pub/dengcy/GWAS_Multiomics/Distribution_test")
load("/share/pub/dengcy/GWAS_Multiomics/banchmarkresult/scPgwas/scpagwas_HCA_Migraine.RData")
Pagwas$Pathway_ld_gwas_data<-NULL
as_matrix <- function(mat) {
  tmp <- matrix(data = 0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  row_pos <- mat@i + 1
  col_pos <- findInterval(seq(mat@x) - 1, mat@p[-1]) + 1
  val <- mat@x
  for (i in seq_along(val)) {
    tmp[row_pos[i], col_pos[i]] <- val[i]
  }
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}


data_mat <- Pagwas$data_mat
data_mat2 <-Pagwas$pca_scCell_mat

rm(Pagwas)
gene_sparsity<-apply(data.matrix(1:nrow(data_mat)),1,function(x){
  a<-data_mat[x,]
  return(sum(a==0)/length(a))
})

cellsingle_sparsity<-apply(data.matrix(1:ncol(data_mat),nrow=1),2,function(x){
  a<-data_mat[,x]  
  return(sum(a==0)/length(a))
})
gene_sparsity2<-round(gene_sparsity,3)
cellsingle_sparsity2<-round(cellsingle_sparsity,3)
Pathway_sparsity<-apply(data_mat2,1,function(x){
  x<-round(x,2)
  return(sum(x==0)/length(x))
})

cellPathway_sparsity<-apply(data_mat2,2,function(x){
  x<-round(x,2)
  return(sum(x==0)/length(x))
})
Pathway_sparsity2<-round(Pathway_sparsity,4)
cellPathway_sparsity2<-round(cellPathway_sparsity,4)

df<- data.frame(gene_sparsity2)
#pdf("gene_sparsity_distribution.pdf",width = 4,height = 4)
p1<-ggplot(df, aes(x=gene_sparsity2))+
  geom_density(fill="#3A5BA0",color="#3A5BA0",alpha=0.8)+
  scale_color_brewer()+
  theme_classic()+
  labs(title="",x="Sparsity of gene(%)")+
  theme(legend.position="none")
#dev.off()


df2<- data.frame(Pathway_sparsity2)

#pdf("Pathway_sparsity_distribution.pdf",width = 4,height = 4)
p2<-ggplot(df2, aes(x=Pathway_sparsity2))+
  geom_density(fill="#FF5B00",color="#FF5B00",alpha=0.8)+
  scale_color_brewer()+
  theme_classic()+
  labs(title="",x="Sparsity of pathway(%)")+
  theme(legend.position="none")
#dev.off()


df3<- data.frame(cellsingle_sparsity2)
#pdf("cellsingle_sparsity_distribution.pdf",width = 4,height = 4)
p3<-ggplot(df3, aes(x=cellsingle_sparsity2))+
  geom_density(fill="#3A5BA0",color="#3A5BA0",alpha=0.8)+
  scale_color_brewer()+
  theme_classic()+
  labs(title="",x="Sparsity of cell in gene profile(%)")+
  theme(legend.position="none")
#dev.off()


df4<- data.frame(cellPathway_sparsity2)

#pdf("cellPathway_sparsity_distribution.pdf",width = 4,height =4)
p4<-ggplot(df4, aes(x=cellPathway_sparsity2))+
  geom_density(fill="#FF5B00",color="#FF5B00",alpha=0.8)+
  scale_color_brewer()+
  theme_classic()+
  labs(title="",x="Sparsity of cell in pathway profile(%)")+
  theme(legend.position="none")

  median<- apply(data_mat2,1,median)
#gg_pca<-as.data.frame(pca_scCell_mat[,sample(1:ncol(pca_scCell_mat),100)])
gg_pca<-data.frame(median=median,pathway=rownames(data_mat2))
#gg_pca$pathway<-rownames(gg_pca)

#pdf("Pathway_gene_distribution.pdf",width = 4,height = 6)
p5<-ggplot(gg_pca, aes(x=median))+
  geom_density(fill="#FF5B00",color="#FF5B00",alpha=0.8)+
  theme_classic()+
labs(title="pathway",x="Median expression")+
 theme(legend.position="none")
#dev.off()

median2<- apply(data_mat,1,median)
gg_scgene<-data.frame(median=median2,gene=rownames(data_mat))
#library(reshape2)
#gg_scgene<-reshape2::melt(gg_scgene,id.vars="gene")

#pdf("singlegene_distribution.pdf",width = 4,height = 6)
p6<-ggplot(gg_scgene, aes(x=median))+
  geom_density(fill="#3A5BA0",color="#3A5BA0",alpha=0.8)+
  theme_classic()+
labs(title="gene",x="Median expression")+
 theme(legend.position="none")
#dev.off()
#dev.off()
pdf("HCA_cellPathway_sparsity_distribution.pdf",width =7,height =5)
patchwork::wrap_plots(plotlist = list(p1,p3,p6,p2,p4,p5), ncol = 3)
dev.off()

```

