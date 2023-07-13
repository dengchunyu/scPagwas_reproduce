
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/CorrectBGp")
load("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_groundtruth_bmmc_monocyte_v10.RData")
all_fortify_can <- fortify.Seurat.umap(Pagwas_groundtruth)
all_fortify_can$CorrectBG_adj_p <- correct_pdf$adj_p
all_fortify_can$CorrectBG_z <- correct_pdf$pooled_z
plots_sigp2 <- ggplot() +
    geom_point(
    data = all_fortify_can[all_fortify_can$CorrectBG_adj_p > p_thre, ],
    aes(x = UMAP_1, y = UMAP_2), size = size, alpha = 0.8,
    color = "#E4DCCF"
    ) +
    umap_theme() +
    # new_scale_color() +
    geom_point(
    data = all_fortify_can[all_fortify_can$CorrectBG_adj_p <= p_thre, ],
    aes(x = UMAP_1, y = UMAP_2), color = "#EA5455", size =size
    ) +
    umap_theme() +
    # new_scale_color() +
    ggtitle(paste0("CorrectBG adjp value <", p_thre, " significant cells"))

df$CorrectBG_adj_p<-ifelse(df$CorrectBG_adj_p<=0.05,1,0)
df$celltypes<-ifelse(df$celltypes %in% c("12_CD14.Mono.2","11_CD14.Mono.1","13_CD16.Mono"),1,0)
#CorrectBG_adj_p的假阳性率
sum(df$CorrectBG_adj_p==1&df$celltypes==0)/sum(df$CorrectBG_adj_p==1)

##cell types

library(scPagwas)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/CelltypeP")
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")

Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data =Pagwas,
                     output.prefix="celltype", 
                     output.dirs="model_monocytecount_celltype",
                     block_annotation = block_annotation,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     singlecell=F,
                     celltype=T
)
save(Pagwas,file="Pagwas_monocytecount_modeldata_celltype.RData")
Pagwas$bootstrap_results$bp_value<-p.adjust(Pagwas$bootstrap_results$bp_value,method = "fdr")
pdf("model_monocytecount_celltype.pdf")
Bootstrap_estimate_Plot(bootstrap_results=Pagwas$bootstrap_results,
                        width = 9,
                        height = 6,
                        do_plot=T)
dev.off()
df<-Pagwas$bootstrap_results[-1,]
write.table(df,file="model_monocytecount_celltypeP.txt",sep="\t",quote=F,row.names = F)



