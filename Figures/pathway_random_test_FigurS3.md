### 2.1.1Building pathway data randomly from all gene expression profiles based on KEGG data can study the interactions between genes and pathways.


```R
library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
#Genes_by_pathway_kegg
modeldata<-readRDS("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds")
genes_all<-rownames(modeldata)
set.seed(1234)
Genes_by_pathway_kegg_random<-lapply(1:100,function(i){
    pl<-lapply(Genes_by_pathway_kegg,function(j){
        sample(genes_all,length(j))
    })
return(pl)
})
save(Genes_by_pathway_kegg_random,file="Genes_by_pathway_kegg_random.RData")
```

scPagwas

```R
library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result")
load("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Genes_by_pathway_kegg_random.RData")
for(i in 1:100){
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_random_kegg", 
                     output.prefix=paste0("kegg_",i),
                     iters_singlecell = 10,
                     block_annotation = block_annotation,
                     Pathway_list=Genes_by_pathway_kegg_random[[i]],
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file=paste0("model_random_kegg_",i,".RData"))
}

```

To calculate the accuracy of the results, the proportion of mononuclear cells among the top 50% of genes is calculated, as well as the number of significantly enriched pathways among the top 1000 genes.

```R
library(scPagwas)
library(Seurat)
library(WebGestaltR)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result")
percent_list<-c()
refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
pa_num<-list()
outputDirectory <- getwd()
for(i in 1:100){
    load(paste0("model_random_kegg_",i,".RData"))
    df<-Pagwas@meta.data
    df<-df[order(df$scPagwas.TRS.Score1,decreasing = T),]
    df<-df[1:round(nrow(df)/2),]
    df<-table(df$celltype)
    percent_list[i]<-df["monocytes"]/sum(df)
    scPagwas_topgenes <- names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation, decreasing = T), ])[1:1000]
    result<-WebGestaltR(enrichMethod="ORA", organism="hsapiens",
        enrichDatabase="geneontology_Biological_Process_noRedundant", 
        interestGene=scPagwas_topgenes,
        interestGeneType="genesymbol", 
        referenceGeneFile=refFile,
        referenceGeneType="genesymbol", 
        isOutput=FALSE,
        projectName=paste0("random_",i))
        a1<-sum(result$FDR<0.01)
        a2<-sum(result$FDR<0.001)
        a3<-sum(result$FDR<0.0001)
        a4<-sum(result$FDR<0.00001)
        pa_num[[i]]<-c(a1,a2,a3,a4)
}
save(percent_list,file="random_kegg_percent_list.RData")
save(pa_num,file="random_kegg_pa_num.RData")

```
### 2.1.2 visual

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result")
load("random_kegg_percent_list.RData")
load("random_kegg_pa_num.RData")
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(ggstatsplot)
num_l<-unlist(lapply(pa_num,function(i) i[1]))
df<-data.frame(percent=percent_list,pa_num=num_l)
library(ggstatsplot)
p<-ggscatterstats(df, 
                x = percent, y = pa_num,
               type = "pearson",
               centrality.para = "mean",    
               margins = "both",
               xfill = "#009E73",
               yfill = "#D55E00",
               xlab = "Percent of monocytes in top 50% cells (Sensitivity)",
          ylab = "Number of significant pathways for top 1000 genes")

pdf("random_kegg_percent.pdf",width=8,height=8)
p
dev.off()
```