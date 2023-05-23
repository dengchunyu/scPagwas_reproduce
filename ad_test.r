#读取之前的结果数据
ad2<-readRDS("/share/pub/dengcy/GWAS_Multiomics/ad_test/pagwas_ad_GSE160936.rds")
#基于单细胞矩阵构造单细胞数据
library(Seurat)
library(scPagwas)
Single_data<-Seurat::CreateSeuratObject(
  counts=ad2$data_mat,
  assay = "RNA",
  meta.data=ad2$Celltype_anno
)
Idents(Single_data)<-Single_data$annotation
saveRDS(Single_data,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/GSE160936.rds")
Pagwas<-scPagwas_main(
    Pagwas = NULL,
                     gwas_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/AD_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/ad_test/GSE160936.rds",
            block_annotation = block_annotation,
            Pathway_list=Genes_by_pathway_kegg,
            chrom_ld = chrom_ld,
            iters_singlecell = 100,
output.dirs="GSE160936",
output.prefix="ad2",
            singlecell=T,
            celltype=T
)
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/ad_test/pagwas_ad_GSE160936.rds")