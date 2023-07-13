# Covid19 traits for PBMC scRNA-seq data for scPagwas

### Severe

> covid
> An object of class Seurat 
> 37985 features across 104484 samples within 2 assays 
> Active assay: SCT (17803 features, 0 variable features)
>  1 other assay present: RNA



### moderate

> An object of class Seurat 
> 37985 features across 204508 samples within 2 assays 
> Active assay: SCT (17803 features, 0 variable features)
>  1 other assay present: RNA



### mild

> An object of class Seurat 
> 37985 features across 122686 samples within 2 assays 
> Active assay: SCT (17803 features, 0 variable features)
>  1 other assay present: RNA



### Normal

> An object of class Seurat 
> 37985 features across 37775 samples within 2 assays 
> Active assay: SCT (17803 features, 0 variable features)
>  1 other assay present: RNA

## scPagwas

The scPagwas and other results for Figure5

### 1. severe

#### Run celltypes scPagwas result

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 setwd("/share/pub/dengcy/GWAS_Multiomics/test/covid19/test5.6")
 for(i in c("severe","Normal","moderate","mild")){
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                       assay="RNA",
                     block_annotation = block_annotation,
                     Single_data =paste0("/share2/pub/jiangdp/jiangdp/COVID/data/",i,"_all.rds"),
                     singlecell=F,
                      output.prefix=i,
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
Pagwas$Pathway_ld_gwas_data<-NULL
save(Pagwas,file=paste0("scpagwas_",i,".RData"))
}
#single cell result for severe
load("scpagwas_severe.RData")
Pagwas<-scPagwas_main(Pagwas = Pagwas,
                     gwas_data ="/share2/pub/jiangdp/jiangdp/COVID/data/COVID19_GWAS.txt",
                       assay="RNA",
                     block_annotation = block_annotation,
                     singlecell=T,
                      output.prefix="severe",
                     celltype=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld)
save(Pagwas,file="scpagwas_severe.RData")

Bootstrap_P_Barplot(p_results=Pagwas$bootstrap_results$bp_value[-1],
                                p_names=rownames(Pagwas$bootstrap_results)[-1],
                                title = "severe",
                                figurenames = "barplot_severe_kegg.pdf",
                                width = 5,
                                height = 7,
                                do_plot=F)

```

#### MAGMA

```R
gwas<-bigreadr::fread2("/share/pub/dengcy/GWAS_Multiomics/test/covid19/COVID19_GWAS.txt")
write.table(gwas,file="/share/pub/dengcy/GWAS_Multiomics/test/covid19/COVID19_GWAS.txt",quote=F)

gwas$N <- 
magma_Input1<-gwas[,c("rsid","pvalue","N")]
magma_Input2<-gwas[,c("rsid","chrom","pos","pos")]
write.table(magma_Input2,file=paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/magma_Input2.txt"),sep="\t",row.names=F,quote=F,col.names=F)
colnames(magma_Input1)<-c("SNP","P","N")
write.table(magma_Input1,file=paste0("/share/pub/dengcy/GWAS_Multiomics/test/covid19/magma_Input1.txt"),sep="\t",row.names=F,quote=F)
```

```shell
cd /share/pub/dengcy/Singlecell/COVID19/MAGMA
./magma --annotate window=10,10 --snp-loc /share/pub/dengcy/GWAS_Multiomics/test/covid19/magma_Input2.txt \
--gene-loc /share/pub/dengcy/Singlecell/COVID19/MAGMA/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded \
--out /share/pub/dengcy/GWAS_Multiomics/test/covid19/annotated_10kbup_10_down

./magma --bfile /share/pub/dengcy/Singlecell/COVID19/MAGMA/g1000_eur/g1000_eur \
--pval /share/pub/dengcy/GWAS_Multiomics/test/covid19/magma_Input1.txt ncol=3 \
--gene-annot /share/pub/dengcy/GWAS_Multiomics/test/covid19/annotated_10kbup_10_down.genes.annot \
--out /share/pub/dengcy/GWAS_Multiomics/test/covid19/annotated_10kbup_10down
```
