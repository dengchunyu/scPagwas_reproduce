## Increase the random proportion in pathways from 0.05 to 0.2 at an increment of 0.05 each time and calculate the results.


```R
modeldata <- readRDS("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds")
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result2")
library(scPagwas)
percent_list<-c(0.05,0.1,0.15,0.2)
genes_all<-rownames(modeldata)
genes_other <- genes_all[!(genes_all %in% unique(unlist(Genes_by_pathway_kegg)))]
percent_list<-c(0.15,0.25)
for(i in percent_list){
Genes_by_pathway_kegg_random<-lapply(1:20,function(i){
    pl<-lapply(Genes_by_pathway_kegg,function(j){
        num <- floor(i*length(j))
        other <- sample(genes_other,num)
        j<-c(j,other)
    })
return(pl)
})
save(Genes_by_pathway_kegg_random,file=paste0("Genes_by_pathway_kegg_random_",i,".RData"))
}

```

## scPagwas for model data

```R
library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result2")
Args <- commandArgs(T)
#gwas='finngen_r7_c3_breast_exallc'
i = print(Args[1])
percent_list<-c(0.05,0.1,0.15,0.2)
#for(i in percent_list){
i<-0.15
    load(paste0("Genes_by_pathway_kegg_random_",i,".RData"))
    TRS_score<-list()
    gPas_score<-list()
    for(j in 1:20){
        Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs=paste0("random_kegg_",i,"_",j),
output.prefix="",
            block_annotation = block_annotation,
            Pathway_list=Genes_by_pathway_kegg_random[[j]],
            chrom_ld = chrom_ld,
            iters_singlecell = 10,
            singlecell=T,
            celltype=F
        )
        TRS_score[j]<-Pagwas@meta.data$scPagwas.TRS.Score1
        gPas_score[j]<-Pagwas@meta.data$scPagwas.gPAS.score
    }
    save(TRS_score,gPas_score,file=paste0("model_random_kegg_TRS_",i,".RData"))
#}
```


```shell
#!/bin/bash
#SBATCH -e model9.err
#SBATCH -o model9.out
#SBATCH -J model9
#SBATCH -w in007
#SBATCH --mem=50000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate R4.2
cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result2
Rscript 1.r 0.1

```

## plot

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result2")
library(scPagwas)
library(Seurat)
library(WebGestaltR)
refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
outputDirectory <- getwd()
#percent_list<-seq(0.1,1,0.1)
Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds"
scdata<-readRDS(Single_data)
anno<-scdata$celltype
percent_list<-c(0.05,0.1,0.15,0.2)
random_kegg_sensitivity <- lapply(percent_list,function(i){
    load(paste0("model_random_kegg_TRS_",i,".RData"))
    pcl<-lapply(TRS_score,function(x){
      df<-data.frame(TRS=x,celltype=anno)
    df<-df[order(df$TRS,decreasing = T),]
    df<-df[1:round(nrow(df)/2),]
    df<-table(df$celltype)
    pc<-df["monocytes"]/sum(df)
    return(pc)
    })
return(unlist(pcl))
})
random_kegg_sensitivity<-as.data.frame(do.call(rbind,random_kegg_sensitivity))
random_kegg_sensitivity<-t(as.matrix(random_kegg_sensitivity))
random_kegg_sensitivity<-as.data.frame(random_kegg_sensitivity)
colnames(random_kegg_sensitivity)<- paste0("addnoise",percent_list)
random_kegg_sensitivity$RandomTimes<-paste0("times",1:20)
library(reshape2)
random_kegg_sensitivity<-melt(random_kegg_sensitivity,id.vars = "RandomTimes")
random_kegg_sensitivity$noise_percent<-rep(percent_list,20)
save(random_kegg_sensitivity,file="random_kegg_sensitivity.RData")

setwd("D:/OneDrive/GWAS_Multiomics/Manuscripts/Revise_comments/pathway_add/Random_kegg_result2")
load("random_kegg_sensitivity.RData")
library(ggplot2)
library(ggpubr)
p<-ggboxplot(random_kegg_sensitivity, x = "variable", y = "value",color="#60281E",add = "jitter",size = 0.5,palette = "jco",ylab = "Sensitivity",xlab = "Noise percent",ggtheme = theme_bw())

pdf("random_kegg_noise_sensitivity.pdf")
p
dev.off()
png("random_kegg_noise_sensitivity.png")
p
dev.off()
```