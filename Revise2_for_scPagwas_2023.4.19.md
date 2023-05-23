# 第二个reviewer的问题

## 1. 加入三个新的通路分析结果

```R
remove.packages("scPagwas")
install.packages("/share/pub/dengcy/software/scPagwas_1.1.1.tar.gz",repos=NULL,type="source")
#genes.by.reactome.pathway.RData
library(scPagwas)
library(Seurat)
#length(genes.by.reactome.pathway)
#[1] 1615
set.seed(123)
reduce_genes.by.reactome.pathway<-genes.by.reactome.pathway
#删除小于50大于300个基因的list
reduce_genes.by.reactome.pathway<-reduce_genes.by.reactome.pathway[sapply(reduce_genes.by.reactome.pathway,length)>50]
reduce_genes.by.reactome.pathway<-reduce_genes.by.reactome.pathway[sapply(reduce_genes.by.reactome.pathway,length)<300]

reduce_genes.by.reactome.pathway<-scPagwas::reduce_pathway(
  pathway_seed=names(reduce_genes.by.reactome.pathway)[sample(1:length(reduce_genes.by.reactome.pathway),50)],
                                                 pathway_list=reduce_genes.by.reactome.pathway,
                                                 remove_proporion=0.8)
length(reduce_genes.by.reactome.pathway)
length(unique(unlist(reduce_genes.by.reactome.pathway)))
save(reduce_genes.by.reactome.pathway,file="reduce_genes.by.reactome.pathway.RData")

setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
load("reduce_genes.by.reactome.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.prefix="reactome", 
                     output.dirs="model_reduce_reactome",
                     block_annotation = block_annotation,
                     Pathway_list=reduce_genes.by.reactome.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_model_reduce_reactome.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.prefix="reactome", 
                     output.dirs="model_reactome",
                     block_annotation = block_annotation,
                     Pathway_list=genes.by.reactome.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_model_reactome.RData")
#genes.by.regulatory.pathway.RData

set.seed(123)
reduce_genes.by.regulatory.pathway<-genes.by.regulatory.pathway
#删除genes.by.regulatory.pathway中少于50个基因的list
reduce_genes.by.regulatory.pathway<-reduce_genes.by.regulatory.pathway[sapply(reduce_genes.by.regulatory.pathway,length)>50]
#删除大于300个基因的list
reduce_genes.by.regulatory.pathway<-reduce_genes.by.regulatory.pathway[sapply(reduce_genes.by.regulatory.pathway,length)<300]

for(i in 1:10){
reduce_genes.by.regulatory.pathway<-scPagwas::reduce_pathway(
  pathway_seed=names(reduce_genes.by.regulatory.pathway)[sample(1:length(reduce_genes.by.regulatory.pathway),50)],
                                                 pathway_list=reduce_genes.by.regulatory.pathway,
                                                 remove_proporion=0.3)
}
length(reduce_genes.by.regulatory.pathway)
length(unique(unlist(reduce_genes.by.regulatory.pathway)))
save(reduce_genes.by.regulatory.pathway,file="reduce_genes.by.regulatory.pathway.RData")

setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
load("reduce_genes.by.regulatory.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_reduce_regulatory", 
                     output.prefix="regulatory",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=reduce_genes.by.regulatory.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_model_reduce_regulatory.RData")

Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_regulatory", 
                     output.prefix="regulatory",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=genes.by.regulatory.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_model_regulatory.RData")

#genes.by.tft.pathway.RData
reduce_genes.by.tft.pathway<-genes.by.tft.pathway
set.seed(123)
#删除genes.by.tft.pathway中少于50个基因的list
reduce_genes.by.tft.pathway<-reduce_genes.by.tft.pathway[sapply(reduce_genes.by.tft.pathway,length)>50]
#删除大于300个基因的list
reduce_genes.by.tft.pathway<-reduce_genes.by.tft.pathway[sapply(reduce_genes.by.tft.pathway,length)<300]

reduce_genes.by.tft.pathway<-scPagwas::reduce_pathway(
  pathway_seed=names(reduce_genes.by.tft.pathway)[sample(1:length(reduce_genes.by.tft.pathway),50)],
                                                 pathway_list=reduce_genes.by.tft.pathway,
                                                 remove_proporion=0.5)


length(reduce_genes.by.tft.pathway)
length(unique(unlist(reduce_genes.by.tft.pathway)))
save(reduce_genes.by.tft.pathway,file="reduce_genes.by.tft.pathway.RData")

setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
load("reduce_genes.by.tft.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_reduce_tft", 
                     output.prefix="tft",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=reduce_genes.by.tft.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_model_reduce_tft.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_tft", 
                     output.prefix="tft",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=genes.by.tft.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_model_tft.RData")


#genes.by.tft.pathway.RData
reduce_genes.by.tft.pathway<-genes.by.tft.pathway
set.seed(123)
#删除genes.by.tft.pathway中少于50个基因的list
reduce_genes.by.tft.pathway<-reduce_genes.by.tft.pathway[sapply(reduce_genes.by.tft.pathway,length)>50]
#删除大于300个基因的list
reduce_genes.by.tft.pathway<-reduce_genes.by.tft.pathway[sapply(reduce_genes.by.tft.pathway,length)<300]

reduce_genes.by.tft.pathway<-scPagwas::reduce_pathway(
  pathway_seed=names(reduce_genes.by.tft.pathway)[sample(1:length(reduce_genes.by.tft.pathway),50)],
                                                 pathway_list=reduce_genes.by.tft.pathway,
                                                 remove_proporion=0.5)

#genes.by.gobp.pathway
reduce_genes.by.gobp.pathway<-genes.by.gobp.pathway
set.seed(123)
#删除genes.by.gobp.pathway中少于50个基因的list
reduce_genes.by.gobp.pathway<-reduce_genes.by.gobp.pathway[sapply(reduce_genes.by.gobp.pathway,length)>50]
#删除大于300个基因的list
reduce_genes.by.gobp.pathway<-reduce_genes.by.gobp.pathway[sapply(reduce_genes.by.gobp.pathway,length)<300]
reduce_genes.by.gobp.pathway<-scPagwas::reduce_pathway(
  pathway_seed=names(reduce_genes.by.gobp.pathway)[sample(1:length(reduce_genes.by.gobp.pathway),50)],
                                                 pathway_list=reduce_genes.by.gobp.pathway,
                                                 remove_proporion=0.5)
save(reduce_genes.by.gobp.pathway,file="reduce_genes.by.gobp.pathway.RData")

setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
load("reduce_genes.by.gobp.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_reduce_gobp", 
                     output.prefix="gobp",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=reduce_genes.by.gobp.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_model_reduce_gobp.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_gobp", 
                     output.prefix="gobp",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=genes.by.gobp.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_model_gobp.RData")

#genes.by.immunologic.pathway.RData
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
load("reduce_genes.by.immunologic.pathway.RData")
reduce_genes.by.immunologic.pathway<-as.list(reduce_genes.by.immunologic.pathway)
class(reduce_genes.by.immunologic.pathway)
#将class为array的list转换为list
reduce_genes.by.immunologic.pathway<-lapply(reduce_genes.by.immunologic.pathway,as.character)
save(reduce_genes.by.immunologic.pathway,file="reduce_genes.by.immunologic.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_reduce_immunologic", 
                     output.prefix="immunologic",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=reduce_genes.by.immunologic.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_model_reduce_immunologic.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_immunologic", 
                     output.prefix="immunologic",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=genes.by.immunologic.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_model_immunologic.RData")

#genes.by.immunesigdb.pathway.RData
reduce_genes.by.immunesigdb.pathway<-genes.by.immunesigdb.pathway
set.seed(1234)
#删除genes.by.immunesigdb.pathway中少于50个基因的list
reduce_genes.by.immunesigdb.pathway<-reduce_genes.by.immunesigdb.pathway[sapply(reduce_genes.by.immunesigdb.pathway,length)>50]
#删除大于300个基因的list
reduce_genes.by.immunesigdb.pathway<-reduce_genes.by.immunesigdb.pathway[sapply(reduce_genes.by.immunesigdb.pathway,length)<300]
#去冗余
for(i in 1:50){
reduce_genes.by.immunesigdb.pathway<-scPagwas::reduce_pathway(
  pathway_seed=names(reduce_genes.by.immunesigdb.pathway)[sample(1:length(reduce_genes.by.immunesigdb.pathway),200)],
                                                 pathway_list=reduce_genes.by.immunesigdb.pathway,
                                                 remove_proporion=0.1)
print(length(reduce_genes.by.immunesigdb.pathway))
}

length(reduce_genes.by.immunesigdb.pathway)
length(unique(unlist(reduce_genes.by.immunesigdb.pathway)))
save(reduce_genes.by.immunesigdb.pathway,file="reduce_genes.by.immunesigdb.pathway.RData")
###immunesigdb
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
load("reduce_genes.by.immunesigdb.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_reduce_immunesigdb", 
                     output.prefix="immunesigdb",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=reduce_genes.by.immunesigdb.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_model_reduce_immunesigdb.RData")

Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_immunesigdb", 
                     output.prefix="immunesigdb",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=genes.by.immunesigdb.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T)
save(Pagwas,file="Pagwas_model_immunesigdb.RData")

##reduce_genes.by.celltype.pathway.RData
library(scPagwas)
library(Seurat)
reduce_genes.by.celltype.pathway<-genes.by.celltype.pathway
set.seed(1234)
#删除genes.by.celltype.pathway中少于50个基因的list
reduce_genes.by.celltype.pathway<-reduce_genes.by.celltype.pathway[sapply(reduce_genes.by.celltype.pathway,length)>50]
#删除大于300个基因的list
reduce_genes.by.celltype.pathway<-reduce_genes.by.celltype.pathway[sapply(reduce_genes.by.celltype.pathway,length)<300]
#去冗余
for(i in 1:10){
reduce_genes.by.celltype.pathway<-scPagwas::reduce_pathway(
  pathway_seed=names(reduce_genes.by.celltype.pathway)[sample(1:length(reduce_genes.by.celltype.pathway),100)],
                                                 pathway_list=reduce_genes.by.celltype.pathway,
                                                 remove_proporion=0.5)

print(length(reduce_genes.by.celltype.pathway))
}
save(reduce_genes.by.celltype.pathway,file="reduce_genes.by.celltype.pathway.RData")
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
load("reduce_genes.by.celltype.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_reduce_celltype", 
                     output.prefix="celltype",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=reduce_genes.by.celltype.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_model_reduce_celltype.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_celltype", 
                     output.prefix="celltype",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=genes.by.celltype.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T)
save(Pagwas,file="Pagwas_model_celltype.RData")

#kegg
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
reduce_Genes_by_pathway_kegg<-Genes_by_pathway_kegg
set.seed(1234)
#删除genes.by.celltype.pathway中少于50个基因的list
reduce_Genes_by_pathway_kegg<-reduce_Genes_by_pathway_kegg[sapply(reduce_Genes_by_pathway_kegg,length)>50]
#删除大于300个基因的list
reduce_Genes_by_pathway_kegg<-reduce_Genes_by_pathway_kegg[sapply(reduce_Genes_by_pathway_kegg,length)<300]
#去冗余
#for(i in 1:10){
reduce_Genes_by_pathway_kegg<-scPagwas::reduce_pathway(
  pathway_seed=names(reduce_Genes_by_pathway_kegg)[sample(1:length(reduce_Genes_by_pathway_kegg),10)],
                                                 pathway_list=reduce_Genes_by_pathway_kegg,
                                                 remove_proporion=0.5)

print(length(reduce_Genes_by_pathway_kegg))
#}
save(reduce_Genes_by_pathway_kegg,file="reduce_Genes_by_pathway_kegg.RData")
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
load("reduce_Genes_by_pathway_kegg.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_reduce_kegg", 
                     output.prefix="kegg",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=reduce_Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_model_reduce_kegg.RData")

#################################
##real data
###################################
load("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_groundtruth.monocyte.RData")
Pagwas_groundtruth@meta.data<-Pagwas_groundtruth@meta.data[,c("celltypes","seurat_clusters")]
saveRDS(Pagwas_groundtruth,file="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds")
library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
load("reduce_genes.by.reactome.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_reduce_reactome", 
                     output.prefix="reactome",
                    iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=reduce_genes.by.reactome.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_reduce_reactome.RData")

Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_reactome", 
                     output.prefix="reactome",
                    iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=genes.by.reactome.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_reactome.RData")

library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
load("reduce_genes.by.regulatory.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_reduce_regulatory", 
                     output.prefix="regulatory",
                     block_annotation = block_annotation,
                                          iters_singlecell = 500,
                     Pathway_list=reduce_genes.by.regulatory.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_reduce_regulatory.RData")

Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_regulatory", 
                     output.prefix="regulatory",
                     block_annotation = block_annotation,
                                          iters_singlecell = 500,
                     Pathway_list=genes.by.regulatory.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_regulatory.RData")

library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
load("reduce_genes.by.tft.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_reduce_tft", 
                     output.prefix="tft",
                     block_annotation = block_annotation,
                                          iters_singlecell = 500,
                     Pathway_list=reduce_genes.by.tft.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_reduce_tft.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_tft", 
                     output.prefix="tft",
                     block_annotation = block_annotation,
                                          iters_singlecell = 500,
                     Pathway_list=genes.by.tft.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_tft.RData")

library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
load("reduce_genes.by.gobp.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_reduce_gobp", 
                     output.prefix="gobp",
                     block_annotation = block_annotation,
                                          iters_singlecell = 500,
                     Pathway_list=reduce_genes.by.gobp.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_reduce_gobp.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_gobp", 
                     output.prefix="gobp",
                     block_annotation = block_annotation,
                                          iters_singlecell = 500,
                     Pathway_list=genes.by.gobp.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_gobp.RData")
#immunologic
load("reduce_genes.by.immunologic.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_reduce_immunologic", 
                     output.prefix="immunologic",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=reduce_genes.by.immunologic.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_reduce_immunologic.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_immunologic", 
                     output.prefix="immunologic",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=genes.by.immunologic.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_immunologic.RData")
#immunesigdb
load("reduce_genes.by.immunesigdb.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_reduce_immunesigdb", 
                     output.prefix="immunesigdb",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=reduce_genes.by.immunesigdb.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_reduce_immunesigdb.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_immunesigdb", 
                     output.prefix="immunesigdb",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=genes.by.immunesigdb.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_immunesigdb.RData")

#celltype
load("reduce_genes.by.celltype.pathway.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_reduce_celltype", 
                     output.prefix="celltype",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=reduce_genes.by.celltype.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_reduce_celltype.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_celltype", 
                     output.prefix="celltype",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=genes.by.celltype.pathway,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_celltype.RData")
#kegg
load("reduce_Genes_by_pathway_kegg.RData")
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds",
                     output.dirs="real_reduce_kegg", 
                     output.prefix="kegg",
                     iters_singlecell = 500,
                     block_annotation = block_annotation,
                     Pathway_list=reduce_Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file="Pagwas_real_reduce_kegg.RData")

```

```shell
#!/bin/bash
#SBATCH -e model1.err
#SBATCH -o model1.out
#SBATCH -J model1
#SBATCH -w in008
#SBATCH --mem=100000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate R4.2
Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/reactome_model.r
Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/reg_model.r
Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/gp_model.r
```

### 1.2 以上通路结果的可视化

```R
library(scPagwas)
library(ggplot2)
library(dplyr)
library(Seurat)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(ggrepel)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
##model数据对没有reduce的通路结果可视化
p_list<-list()
for(i in c("kegg","reactome","celltype","gobp","tft")){
  load(paste0("Pagwas_model_",i,".RData"))
  print(i)
  #输出可以换行的文本
  reduce_percent <- length(unlist(Pagwas@misc$Pathway_list))/length(unique(unlist(Pagwas@misc$Pathway_list)))

  df<-Pagwas@meta.data
  #将df按照scPagwas.TRS.Score1的列从大到小排序，获得celltype列中前50%的数据中monocytes元素的占比
  df<-df[order(df$scPagwas.TRS.Score1,decreasing = T),]
  df<-df[1:round(nrow(df)/2),]
  #计算monocytes在celltypes列中的占比
  df<-table(df$celltype)
  Sensitivity<- df["monocytes"]/sum(df)
  #保留小数点后两位
  reduce_percent<-round(reduce_percent,2)
  ti<-paste0("Pathway: ",i,"\n","Pathway numbers: ",length(Pagwas@misc$Pathway_list),"\n","Duplication times: ",reduce_percent,"\n","Sensitivity: ",Sensitivity)

  all_fortify_can <- fortify.Seurat.umap(Pagwas)
  p1<- ggplot() +
      geom_point(data = all_fortify_can,
                  aes(x = UMAP_1, y = UMAP_2,color =scPagwas.TRS.Score1), size = 0.2, alpha = 1) +
      umap_theme() +
    scale_colour_gradient2(name="TRS",low="#8479E1",mid="#F7F5F2",high="#FD5D5D")+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))  +
      ggtitle(ti)
  p_list[[i]]<-p1
}
p<-ggarrange(plotlist = p_list, ncol = 5, nrow =1)
ggsave("pathway_model.png",p,width=20,height=6)
ggsave("pathway_model.pdf",p,width=20,height=6)

#low = "#000957", high = "#EBE645"
##model数据对reduce后的通路结果可视化
p_list<-list()
for(i in c("reduce_kegg","reduce_reactome","reduce_celltype","reduce_gobp","reduce_tft")){
  load(paste0("Pagwas_model_",i,".RData"))

reduce_percent <- length(unlist(Pagwas@misc$Pathway_list))/length(unique(unlist(Pagwas@misc$Pathway_list)))
df<-Pagwas@meta.data
    #将df按照scPagwas.TRS.Score1的列从大到小排序，获得celltype列中前50%的数据中monocytes元素的占比
    df<-df[order(df$scPagwas.TRS.Score1,decreasing = T),]
    df<-df[1:round(nrow(df)/2),]
    #计算monocytes在celltypes列中的占比
    df<-table(df$celltype)
    Sensitivity<- df["monocytes"]/sum(df)
#保留小数点后两位
reduce_percent<-round(reduce_percent,2)
  ti<-paste0("Pathway: ",i,"\n","Pathway numbers: ",length(Pagwas@misc$Pathway_list),"\n","Duplication times: ",reduce_percent,"\n","Sensitivity: ",Sensitivity)
  all_fortify_can <- fortify.Seurat.umap(Pagwas)
p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas.TRS.Score1), size = 0.2, alpha = 1) +
        umap_theme() +
      scale_colour_gradient2(name="TRS",low="#8479E1",mid="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle(ti)
p_list[[i]]<-p1
}
p<-ggarrange(plotlist = p_list,  ncol = 5, nrow =1)
ggsave("pathway_model_reduce.png",p,width=20,height=6)
ggsave("pathway_model_reduce.pdf",p,width=20,height=6)
```

真实数据的可视化

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
p_list<-list()
for(i in c("kegg","reactome","celltype","gobp","tft")){
  load(paste0("Pagwas_real_",i,".RData"))
#输出可以换行的文本
reduce_percent <- length(unlist(Pagwas@misc$Pathway_list))/length(unique(unlist(Pagwas@misc$Pathway_list)))
df<-Pagwas@meta.data
    #将df按照scPagwas.TRS.Score1的列从大到小排序，获得celltype列中前50%的数据中monocytes元素的占比
    df<-df[order(df$scPagwas.TRS.Score1,decreasing = T),]
    df<-df[1:round(nrow(df)/2),]
    #计算monocytes在celltypes列中的占比
    df<-table(df$celltype)
    Sensitivity<- sum(df["11_CD14.Mono.1"],df["12_CD14.Mono.2"],df["13_CD16.Mono"])/sum(df)
#保留小数点后两位
reduce_percent<-round(reduce_percent,2)
  ti<-paste0("Pathway: ",i,"\n","Pathway numbers: ",length(Pagwas@misc$Pathway_list),"\n","Duplication times: ",reduce_percent,"\n","Sensitivity: ",Sensitivity)

  all_fortify_can <- fortify.Seurat.umap(Pagwas)
p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas.TRS.Score1), size = 0.2, alpha = 1) +
        umap_theme() +
      scale_color_gradient2(name="TRS",low="#8479E1",mid="#F7F5F2",high="#FD5D5D") +
        theme(aspect.ratio=1) + 
        #在图中加入一段文字
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle(ti)
p_list[[i]]<-p1
}
p<-ggarrange(plotlist = p_list, ncol = 5, nrow =1)
ggsave("pathway_real.png",p,width=20,height=6)
ggsave("pathway_real.pdf",p,width=20,height=6)

##model数据对reduce后的通路结果可视化
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")

p_list<-list()
for(i in c("reduce_kegg","reduce_reactome","reduce_celltype","reduce_gobp","reduce_tft")){
    load(paste0("Pagwas_real_",i,".RData"))
    
reduce_percent <- length(unlist(Pagwas@misc$Pathway_list))/length(unique(unlist(Pagwas@misc$Pathway_list)))
df<-Pagwas@meta.data
    #将df按照scPagwas.TRS.Score1的列从大到小排序，获得celltype列中前50%的数据中monocytes元素的占比
    df<-df[order(df$scPagwas.TRS.Score1,decreasing = T),]
    df<-df[1:round(nrow(df)/2),]
    #计算monocytes在celltypes列中的占比
    df<-table(df$celltypes)
    Sensitivity<- sum(df["11_CD14.Mono.1"],df["12_CD14.Mono.2"],df["13_CD16.Mono"])/sum(df)
#保留小数点后两位
reduce_percent<-round(reduce_percent,2)
  ti<-paste0("Pathway: ",i,"\n","Pathway numbers: ",length(Pagwas@misc$Pathway_list),"\n","Duplication times: ",reduce_percent,"\n","Sensitivity: ",Sensitivity)
    all_fortify_can <- fortify.Seurat.umap(Pagwas)
p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas.TRS.Score1), size = 0.2, alpha = 1) +
        umap_theme() +
      scale_color_gradient2(name="TRS",low="#8479E1",mid="#F7F5F2",high="#FD5D5D") +
        theme(aspect.ratio=1) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle(ti)
p_list[[i]]<-p1
}
p<-ggarrange(plotlist = p_list, ncol = 5, nrow =1)
ggsave("pathway_real_reduce.png",p,width=20,height=6)
ggsave("pathway_real_reduce.pdf",p,width=20,height=6)

```

### 1.3 以上通路结果的量化

model数据

```R
library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
percent_list<-c()
for(i in c("kegg","reactome","celltype","gobp","tft","regulatory","immunologic","immunesigdb","reduce_kegg","reduce_reactome","reduce_celltype","reduce_gobp","reduce_tft","reduce_regulatory","reduce_immunologic","reduce_immunesigdb")){
    load(paste0("Pagwas_model_",i,".RData"))
    df<-Pagwas@meta.data
    #将df按照scPagwas.TRS.Score1的列从大到小排序，获得celltype列中前50%的数据中monocytes元素的占比
    df<-df[order(df$scPagwas.TRS.Score1,decreasing = T),]
    df<-df[1:round(nrow(df)/2),]
    #计算monocytes在celltypes列中的占比
    df<-table(df$celltype)
    percent_list[i]<-df["monocytes"]/sum(df)
    }
save(percent_list,file="Pagwas_model_percent_list.RData")
```

real数据
  
```R
library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
percent_list<-c()
for(i in c("kegg","reactome","celltype","gobp","tft","regulatory","immunologic","immunesigdb","reduce_kegg","reduce_reactome","reduce_celltype","reduce_gobp","reduce_tft","reduce_regulatory","reduce_immunologic","reduce_immunesigdb")){
    load(paste0("Pagwas_real_",i,".RData"))
    df<-Pagwas@meta.data
    #将df按照scPagwas.TRS.Score1的列从大到小排序，获得celltype列中前50%的数据中monocytes元素的占比
    df<-df[order(df$scPagwas.TRS.Score1,decreasing = T),]
    df<-df[1:round(nrow(df)/2),]
    #计算monocytes在celltypes列中的占比
    df<-table(df$celltypes)
    percent_list[i]<- sum(df["11_CD14.Mono.1"],df["12_CD14.Mono.2"],df["13_CD16.Mono"])/sum(df)
    }
save(percent_list,file="Pagwas_real_percent_list.RData")
```

## 2. 通路的测试

### 2.1 不同通路是否具有较高的相关性
```shell
#安装rpy2.robjects
pip install rpy2
```

```python
#读取Rdata文件
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
r = robjects.r
#设置文件夹
r['setwd']("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
r['load']("Genes_by_pathway_kegg.RData")
Genes_by_pathway_kegg = r['Genes_by_pathway_kegg']
r['load']("genes.by.celltype.pathway.RData")
genes_by_celltype_pathway = r['genes.by.celltype.pathway']
r['load']("genes.by.gobp.pathway.RData")
genes_by_gobp_pathway = r['genes.by.gobp.pathway']
r['load']("genes.by.reactome.pathway.RData")
genes_by_reactome_pathway = r['genes.by.reactome.pathway']
r['load']("genes.by.regulatory.pathway.RData")
genes_by_regulatory_pathway = r['genes.by.regulatory.pathway']
r['load']("genes.by.tft.pathway.RData")
genes_by_tft_pathway = r['genes.by.tft.pathway']
r['load']("genes.by.hallmark.pathway.RData")
genes_by_hallmark_pathway = r['genes.by.hallmark.pathway']
r['load']("genes.by.immunesigdb.pathway.RData")
genes_by_immunesigdb_pathway = r['genes.by.immunesigdb.pathway']
r['load']("genes.by.immunologic.pathway.RData")
genes_by_immunologic_pathway = r['genes.by.immunologic.pathway']

r['load']("reduce_Genes_by_pathway_kegg.RData")
reduce_Genes_by_pathway_kegg = r['reduce_Genes_by_pathway_kegg']
r['load']("reduce_genes.by.celltype.pathway.RData")
reduce_genes_by_celltype_pathway = r['reduce_genes.by.celltype.pathway']
r['load']("reduce_genes.by.gobp.pathway.RData")
reduce_genes_by_gobp_pathway = r['reduce_genes.by.gobp.pathway']
r['load']("reduce_genes.by.reactome.pathway.RData")
reduce_genes_by_reactome_pathway = r['reduce_genes.by.reactome.pathway']
r['load']("reduce_genes.by.regulatory.pathway.RData")
reduce_genes_by_regulatory_pathway = r['reduce_genes.by.regulatory.pathway']
r['load']("reduce_genes.by.tft.pathway.RData")
reduce_genes_by_tft_pathway = r['reduce_genes.by.tft.pathway']
r['load']("reduce_genes.by.immunesigdb.pathway.RData")
reduce_genes_by_immunesigdb_pathway = r['reduce_genes.by.immunesigdb.pathway']
r['load']("reduce_genes.by.immunologic.pathway.RData")
reduce_genes_by_immunologic_pathway = r['reduce_genes.by.immunologic.pathway']
```

```python
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import os
#计算通路list中两个通路的共有基因数占第一个通路的比例
def get_overlap(pathway1,pathway2):
  pathway1=list(set(pathway1))
  pathway2=list(set(pathway2))
  overlap=len(set(pathway1).intersection(set(pathway2)))
  return overlap/len(pathway1)
#1.循环计算通路list中两两通路的共有基因数占第一个通路的比例
pahtway_list=[Genes_by_pathway_kegg,genes_by_celltype_pathway,genes_by_gobp_pathway,genes_by_reactome_pathway,genes_by_regulatory_pathway,genes_by_tft_pathway,genes_by_hallmark_pathway,genes_by_immunesigdb_pathway,genes_by_immunologic_pathway]
#设置输出图片的地址
os.chdir("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
#查看pahtway_list的格式
type(pahtway_list)
#给pahtway_list中的每个通路list命名
pathway_name=["Genes_by_pathway_kegg","genes_by_celltype_pathway","genes_by_gobp_pathway","genes_by_reactome_pathway","genes_by_regulatory_pathway","genes_by_tft_pathway","genes_by_hallmark_pathway","genes_by_immunesigdb_pathway","genes_by_immunologic_pathway"]
for i in range(len(pahtway_list)):
  pahtway_list[i].name=pathway_name[i]
#循环计算通路list中两两通路的共有基因数占第一个通路的比例
for m in pahtway_list:
    pathway_l=m
    overlap=np.zeros((len(pathway_l),len(pathway_l)))
    for i in range(len(pathway_l)):
      for j in range(len(pathway_l)):
        overlap[i,j]=get_overlap(pathway_l[i],pathway_l[j])
    plt.figure(figsize=(10,10))
    sns.heatmap(overlap,cmap="YlGnBu")
    plt.savefig(pathway_l.name+"_overlap.pdf",dpi=300)
    plt.savefig(pathway_l.name+"_overlap.png",dpi=300)
###reduce pathway
pathway_list2=[reduce_Genes_by_pathway_kegg,reduce_genes_by_celltype_pathway,reduce_genes_by_gobp_pathway,reduce_genes_by_reactome_pathway,reduce_genes_by_regulatory_pathway,reduce_genes_by_tft_pathway,reduce_genes_by_hallmark_pathway,reduce_genes_by_immunesigdb_pathway,reduce_genes_by_immunologic_pathway]
pathway_name=["reduce_Genes_by_pathway_kegg","reduce_genes_by_celltype_pathway","reduce_genes_by_gobp_pathway","reduce_genes_by_reactome_pathway","reduce_genes_by_regulatory_pathway","reduce_genes_by_tft_pathway","reduce_genes_by_hallmark_pathway","reduce_genes_by_immunesigdb_pathway","reduce_genes_by_immunologic_pathway"]
for i in range(len(pathway_list2)):
  pathway_list2[i].name=pathway_name[i]

os.chdir("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
#循环计算通路list中两两通路的共有基因数占第一个通路的比例
for m in pathway_list2:
    pathway_l=m
    overlap=np.zeros((len(pathway_l),len(pathway_l)))
    for i in range(len(pathway_l)):
      for j in range(len(pathway_l)):
        overlap[i,j]=get_overlap(pathway_l[i],pathway_l[j])
    plt.figure(figsize=(10,10))
    sns.heatmap(overlap,cmap="YlGnBu")
    plt.savefig(pathway_l.name+"_overlap.pdf",dpi=300)
    plt.savefig(pathway_l.name+"_overlap.png",dpi=300)
```
