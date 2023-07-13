## HCA single celll data

### dowload from GEO
E:/OneDrive/GWAS_Multiomics/HCLdata/GSE134355_RAW.tar

```R

library("stringr") 
library("readr")
library(Seurat)

setwd("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata")
files<-list.files("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata/HCAadult")

Seurat_list<-lapply(files,function(x){
    count<-read_table(x)
    gene<-unlist(count[,1])
    count<-as(as(count[,-1],"matrix"),"dgCMatrix")
    rownames(count)<-gene
    a<- strsplit(x,split = "_",fixed=T)[[1]][2]
    a<-str_replace_all(a,"-","_")
    Seurat_object <- CreateSeuratObject(
               counts = count, 
               min.cells = 3, 
               min.features = 200)
    Idents(Seurat_object)<-rep(a,ncol(Seurat_object))
    return(Seurat_object)
                })
    
Seurat_data<-merge(Seurat_list[[1]],Seurat_list[[2]])
for(i in 3:length(files)){
    Seurat_data<-merge(Seurat_data,Seurat_list[[i]])
}
rm(Seurat_list)
Seurat_data <- NormalizeData(Seurat_data, normalization.method = "LogNormalize", scale.factor = 10000)
saveRDS(Seurat_data,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata/HCA_tissues.rds")
```

scPagwas

```R
library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 #Input pathway gene list, you can construct with youself.
  suppressMessages(library(Seurat))
 HCA_tissues<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata/HCA_tissues.rds")
 Idents(HCA_tissues)<-unlist(lapply(as.vector(Idents(HCA_tissues)),function(x){
  return(substr(x,1,nchar(x)-1))
 }))

b<-as.vector(Idents(HCA_tissues))
b[which(b=="Adult_Kidney4_")]<-"Adult_Kidney"
b[which(b=="Adult_Liver1_")]<-"Adult_Liver"
b[which(b=="Adult_Liver4_")]<-"Adult_Liver"
b[which(b=="Adult_Lung3_")]<-"Adult_Lung"
b[which(b=="Adult_Peripheral_Blood3_")]<-"Adult_Peripheral_Blood"
b[which(b=="Adult_Peripheral_Blood4_")]<-"Adult_Peripheral_Blood"
b[which(b=="Adult_Stomach3_")]<-"Adult_Stomach"
b[which(b=="Adult_Transverse_Colon2_")]<-"Adult_Transverse_Colon"
table(b)
Idents(HCA_tissues)<-b
saveRDS(HCA_tissues,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata/HCA_tissues.rds")

#############
 suppressMessages(library(Seurat))
  library(scPagwas)
  
 HCA_tissues<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata/HCA_tissues.rds")
HCA_tissues_scexpr <- Seurat::AverageExpression(HCA_tissues,assays="RNA")[["RNA"]]
save(HCA_tissues_scexpr,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCA_tissues_scexpr.RData")
rm(HCA_tissues_scexpr)

Pagwas <- list();
class(Pagwas) <- 'Pagwas'
Pagwas <- Single_data_input(Pagwas=Pagwas,
                                assay="RNA",
Single_data=/share/pub/dengcy/GWAS_Multiomics/singlecelldata/HCAdata/HCA_tissues.rds,
                                Pathway_list=Genes_by_pathway_kegg)
  
  Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,n.cores=1,
                                 Pathway_list=Genes_by_pathway_kegg
                                 )
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Pagwas_HCA_tissues.RData")
```

## mouse_brain 
```R
 library(scPagwas)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
 #Input pathway gene list, you can construct with youself.
 data(Genes_by_pathway_kegg)
 #gene annotation files.
 data(block_annotation)
 #LD data
 data(chrom_ld)
 ####先计算平均表达
   mouse_brain_scexpr <- Seurat::AverageExpression(mouse_brain_seu2,assays="RNA")[["RNA"]]
save(mouse_brain_scexpr,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/mouse_brain_scexpr.RData")
 Pagwas <- list();
  class(Pagwas) <- 'Pagwas'
    Pagwas <- Single_data_input(Pagwas=Pagwas,
                                assay="RNA",
                                nfeatures =NULL,
                                Single_data=mouse_brain_seu2,
                                Pathway_list=Genes_by_pathway_kegg,
                                min_clustercells=10)
  
  Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,n.cores=1,
                                 Pathway_list=Genes_by_pathway_kegg
                                 )
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Pagwas_mouse_brain.RData")

```
### Gtex tissue data

```R
  Gtex_tissue_seu<-NormalizeData(Gtex_tissue_seu,assay = "RNA")
save(Gtex_tissue_seu,file="/share/pub/dengcy/GWAS_Multiomics/pagwas/GtexRnaseq/Gtex_tissue_seu.RData")

 suppressMessages(library(Seurat))
 library(scPagwas)
 data(Genes_by_pathway_kegg)
 data(block_annotation)
 data(chrom_ld)
 load("/share/pub/dengcy/GWAS_Multiomics/pagwas/GtexRnaseq/Gtex_tissue_seu.RData")

 Pagwas <- list();
 class(Pagwas) <- 'Pagwas'
 Pagwas <- Single_data_input(Pagwas=Pagwas,
                                assay="RNA",
                                nfeatures =NULL,
                                Single_data=Gtex_tissue_seu,
                                Pathway_list=Genes_by_pathway_kegg)
  Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,
                                 Pathway_list=Genes_by_pathway_kegg
                                 )
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Pagwas_GTEX_tissues.RData")
```
