# Groudtruth data for simulated data

We used scDesign2  to simulate a ground truth scRNA-seq dataset containing five cell types of monocytes,  DC, B, NK, T cells for assessing the performance of scPagwas on identifying monocyte count trait-relevant individual cells. DC, a differentiated cell from monocyte [28930664], was setted as non-trait-relevant cell type which can be a confounding factor for monocyte is a powerful celltype for bone marrow and PBMC. In the model-fitting step, we first fitted a multivariate generative model to a real dataset via the fluorescence activated cell sorting (FACS)-sorted bulk hematopoietic populations download from GEO database(GSE107011, https://www.ncbi.nlm.nih.gov/geo/). As there were five sorted cell types, we split the datasets into five subsets depending on the cell type, and fitted a cell type-specific model to each subset. In the data-generation step, we generated a synthetic scRNA-seq data from the fitted model to represent genetic positive cell populations (monocytes) and negative cell populations (DC, B, NK, Tcells) for monocyte count trait. Finally, we obtained 2000 cells synthetic scRNA-seq data  with a cell propotion of 0.5 (monocytes), 0.05 (DC) , 0.2 (B cells), 0.05 (NK) and 0.2 (T cells).



## 1.Construct the simulated data

export OPENBLAS_NUM_THREADS=1

### Get the bulk data for simulating

```R
library(scDesign2)
library(copula)    
library(Rtsne)
library(plyr)     
library(reshape2) 
library(gridExtra)
library(ggpubr)
library(cowplot)
library(ggplot2); 
library(stringr)
theme_set(theme_bw());

set.seed(42)
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")

bulk_tpm<-read.delim("/share/pub/dengcy/GWAS_Multiomics/pbmc_bulkdata/GSE107011_Processed_data_TPM.txt",header=T)
colnames(bulk_tpm)
Mono<-c("DZQV_C_mono","DZQV_I_mono","DZQV_NC_mono" ,"X925L_C_mono" ,"X925L_NC_mono","X9JD4_I_mono","X9JD4_NC_mono","G4YW_I_mono", "G4YW_NC_mono","X9JD4_mDC","G4YW_mDC")
dc<-c("DZQV_pDC","DZQV_mDC", "X925L_pDC","X925L_mDC" , "X9JD4_pDC","G4YW_pDC","G4YW_mDC","G4YW_I_mono", "G4YW_NC_mono")

B<- c("G4YW_B_NSM","G4YW_B_Ex","G4YW_B_SM","G4YW_B_naive","X9JD4_B_Ex","X9JD4_B_SM","X9JD4_B_NSM","X9JD4_B_naive","X925L_B_NSM","X925L_B_Ex","X925L_B_SM","X925L_B_naive","DZQV_B_SM","DZQV_B_naive","DZQV_B_NSM","DZQV_B_Ex")

NK<- c("DZQV_NK","X925L_NK","X9JD4_NK","G4YW_NK")

tcell<-c( "DZQV_CD8_naive","DZQV_CD8_CM","DZQV_CD8_EM","DZQV_CD8_TE","DZQV_TFH","DZQV_Treg","DZQV_Th1","DZQV_Th1.Th17","DZQV_Th17","DZQV_Th2","DZQV_CD4_naive", "X925L_CD8_naive","X925L_CD8_CM","X925L_CD8_EM","X925L_CD8_TE","X925L_TFH","X925L_Treg", "X925L_Th1","X925L_Th1.Th17","X925L_Th17", "X925L_Th2","X925L_CD4_naive","X925L_CD4_TE","X9JD4_CD8_naive", "X9JD4_CD8_CM","X9JD4_CD8_EM","X9JD4_CD8_TE" ,"X9JD4_TFH","X9JD4_Treg", "X9JD4_Th1","X9JD4_Th1.Th17","X9JD4_Th17","X9JD4_Th2","X9JD4_CD4_naive","X9JD4_CD4_TE","G4YW_CD8_naive","G4YW_CD8_CM" ,"G4YW_CD8_EM" ,"G4YW_CD8_TE","G4YW_TFH","G4YW_Treg" ,"G4YW_Th1", "G4YW_Th1.Th17","G4YW_Th17","G4YW_Th2","G4YW_CD4_naive","G4YW_B_naive","G4YW_B_NSM" )

bulk_tpm<-bulk_tpm[,c(Mono,dc,B,NK,tcell)]
genes<-unlist(str_sub(rownames(bulk_tpm),1,15))
bulk_tpm$gene<-genes
bulk_tpm<-bulk_tpm[!duplicated(bulk_tpm$gene),]
rownames(bulk_tpm)<-unique(bulk_tpm$gene)

library(org.Hs.eg.db)
library(dplyr)
a<-data.frame(ensembl_id=bulk_tpm$gene)
g2s=toTable(org.Hs.egSYMBOL)
g2e=toTable(org.Hs.egENSEMBL)
b=merge(a,g2e,by="ensembl_id",all.x=T)
d=merge(b,g2s,by="gene_id",all.x=T)
d<-d[!duplicated(d$symbol),]
d<-d[!is.na(d$symbol),]

bulk_tpm<-bulk_tpm[d$ensembl_id,]
rownames(bulk_tpm)<-d$symbol
#save(bulk_tpm,file="/share/pub/dengcy/GWAS_Multiomics/pbmc_bulkdata/pbmc_bulk_tpm.RData")
bulk_tpm2<-bulk_tpm[,-ncol(bulk_tpm)]
bulk_tpm2<-ceiling(bulk_tpm2)
bulk_tpm2<-as.matrix(bulk_tpm2)
colnames(bulk_tpm2)<-c(rep("monocytes",length(Mono)),
                       rep("DC",length(dc)),
                       rep("B",length(B)),
                       rep("NK",length(NK)),
                       rep("Tcells",length(tcell)))
save(bulk_tpm2,file="/share/pub/dengcy/GWAS_Multiomics/pbmc_bulkdata/pbmc_bulk_tpm3.RData")

library(SeuratObject)
library(Seurat)
 library(scPagwas)
 library(ggplot2)
load("/share/pub/dengcy/GWAS_Multiomics/pbmc_bulkdata/pbmc_bulk_tpm3.RData")
bulk_tpm3<-CreateSeuratObject(counts=bulk_tpm2,assay = "RNA")
bulk_tpm3$celltype<-c(rep("monocytes",length(Mono)),
                      rep("DC",length(dc)),
                      rep("B",length(B)),
                      rep("NK",length(NK)),
                      rep("Tcells",length(tcell)))
bulk_tpm3 <- NormalizeData(bulk_tpm3, normalization.method = "LogNormalize", scale.factor = 10000)
bulk_tpm3 <- ScaleData(bulk_tpm3)
bulk_tpm3<-FindVariableFeatures(object=bulk_tpm3)
Idents(bulk_tpm3)<-bulk_tpm3$celltype
```

### Get the simulated single cell data and normalized

```R
library(scDesign2)
library(copula)    
library(Rtsne)
library(plyr)     
library(reshape2) 
library(gridExtra)
library(ggpubr)
library(cowplot)
library(ggplot2); 
library(stringr)
theme_set(theme_bw());
library(scDesign2)

setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")
load("/share/pub/dengcy/GWAS_Multiomics/pbmc_bulkdata/pbmc_bulk_tpm3.RData")

copula_result <- fit_model_scDesign2(bulk_tpm2,zp_cutoff=0.9,
                                     cell_type_sel=c('monocytes','DC','B','NK','Tcells'), 
                                     sim_method = 'copula',ncores = 20)
saveRDS(copula_result, file = '/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/copula_result_population2.rds')
sim_count_2000 <- simulate_count_scDesign2(copula_result, n_cell_new=2000, sim_method = 'copula',cell_type_prop = c(0.5,0.05,0.2,0.05,0.2))
rownames(sim_count_2000)<-rownames(bulk_tpm2)
saveRDS(sim_count_2000, file = 'sim_count_2000_population.rds')

#
library(SeuratObject)
library(Seurat)
sim_count_2000<-readRDS("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_count_2000_population.rds")

sim_data<-CreateSeuratObject(counts=sim_count_2000,assay = "RNA")
sim_data$celltype<-colnames(sim_count_2000)
sim_data$type<- c(rep(1,1000),rep(0,1000))
sim_data <- NormalizeData(sim_data, normalization.method = "LogNormalize", scale.factor = 10000)
sim_data <- ScaleData(sim_data)
sim_data<-FindVariableFeatures(object=sim_data)

Idents(sim_data)<-sim_data$celltype
#sim_data<-RunTsne(object=sim_data)

sim_data <- FindVariableFeatures(sim_data,nfeatures = 3000)
sim_data <- RunPCA(object = sim_data, assay = "RNA", npcs = 50)

##subfunction:
cluster_pca_umap <- function(obj,assay=NULL, reduction,cluster_res = 0.3){
  #obj2 <- RunPCA(obj, assay = "SCT", reduction = "harmony",verbose = F)
  obj2 <- RunTSNE(object = obj,assay = assay, reduction = reduction, dims = 1:50,check_duplicates = FALSE)
  obj2 <- RunUMAP(object = obj2, assay =assay, reduction = reduction, dims = 1:50,check_duplicates = FALSE)
  obj2 <- FindNeighbors(object=obj2, assay = assay, reduction = reduction, dims = 1:50)
  obj2 <- FindClusters(object=obj2, resolution = cluster_res)
  return(obj2)
}
sim_data<-cluster_pca_umap(obj = sim_data,reduction="pca",cluster_res = 0.3)
saveRDS(sim_data,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds")
Idents(sim_data)<-sim_data$celltype
library(ggpubr)
library(ggplot2)
color_scanpy_patient <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")
pdf("model_groundtruth_monocyte.pdf",width = 6,height = 6)
 DimPlot(sim_data,reduction ="umap",
                         group.by = "celltype",pt.size=0.3,
                         label = F, repel=TRUE)+ 
  umap_theme()+
  scale_colour_manual(name = "celltype", values =color_scanpy_patient[c(2:9,1)]) +
  theme(aspect.ratio=1)
dev.off()

```

### Run scPagwas for monocytecount

```R
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))

setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")
i<-"monocytecount"

Pagwas<-scPagwas_main(Pagwas = NULL,
gwas_data=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"),
                     Single_data =sim_data,
                     output.prefix="modeldata_monocytecount",
                     output.dirs="modeldata_outputv1.10.0",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     seruat_return=T,
                     celltype=F,
                     ncores = 1)
                     
save(Pagwas,file="Pagwas_monocytecount_modeldata_v1.10.0.RData")
```

### Run scDRS

```R
library(SeuratDisk)
library(Seurat)

a<-GetAssayData(Pagwas,assay = "RNA",slot="counts")
Pagwas_modelgroundtruth <- CreateSeuratObject(counts = a,
                           project = "scPagwas",
                           min.cells = 3,
                           min.features = 200)
DefaultAssay(Pagwas_modelgroundtruth) <- "RNA"
SaveH5Seurat(Pagwas_modelgroundtruth, "Pagwas_modelgroundtruth_addata.h5seurat")
Convert("Pagwas_modelgroundtruth_addata.h5seurat", dest="h5ad")

i<-"monocytecount"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_magma_genes.RData"))

scPagwas_genes<-names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation,decreasing=T),])

magmatop1000<-intersect(magma_genes$symbol,rownames(Pagwas))[1:1000]
scPagwastop1000<-scPagwas_genes[1:1000]
 magmatop500<-intersect(magma_genes$symbol,rownames(Pagwas))[1:500]
scPagwastop500<-scPagwas_genes[1:500]
topgene<-list(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500)
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")

a<-data.frame(genes=c(paste(scPagwastop1000,collapse=","),paste(magmatop1000,collapse=","),paste(scPagwastop500,collapse=","),paste(magmatop500,collapse=",")))
rownames(a)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
write.csv(a,file=paste0("scPagwas.magma.Model.topgenes.",i,".csv"))


ab<-lapply(1:100, function(i){
  b<-paste(magma_genes$symbol[1:(10*i)],collapse=",")
  return(b)
})
b<-paste0("MAGMAtop",1:100)
a<-data.frame(genes=unlist(ab))
rownames(a)<-b
write.csv(a,file="magma_scRDS_10_1000.csv")


################################scdrs
import scdrs
from scipy import stats
import pandas as pd
import scanpy as sc
from anndata import AnnData
sc.set_figure_params(dpi=125)
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings

warnings.filterwarnings("ignore")
# load adata
DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/"
H5AD_FILE = os.path.join(DATA_PATH, "Pagwas_modelgroundtruth_addata.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
scdrs.preprocess(adata)


i='monocytecount'
df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_scRDS"+i+".csv", index_col=0)
df_gs = df_gs.loc[["scPagwastop1000","magmatop1000"],:]
df_gs=  df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/"+i+".geneset.gs")

df_score = scdrs.score_cell(
             data=adata,
             gene_list=df_gs['magmatop1000'][0],
             gene_weight=df_gs['magmatop1000'][1],
             ctrl_match_key="mean_var",
             n_ctrl=200,
             weight_opt="vs",
             return_ctrl_raw_score=False,
             return_ctrl_norm_score=True,
             verbose=False)

df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/mono.1000magma.df_modeldata1.csv", sep=",")

df_score = scdrs.score_cell(
             data=adata,
             gene_list=df_gs['scPagwastop1000'][0],
             gene_weight=df_gs['scPagwastop1000'][1],
             ctrl_match_key="mean_var",
             n_ctrl=200,
             weight_opt="vs",
             return_ctrl_raw_score=False,
             return_ctrl_norm_score=True,
             verbose=False)

df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/mono.1000scPagwas.df_modeldata1.csv", sep=",")


method="smultixcan"
gwas="monocytecount"
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/" + method + "/scDRS_result/" + gwas + ".geneset.gs")
n=method + '_top' + str(100)
df_score = scdrs.score_cell(
                data=adata,
                gene_list=df_gs[n][0],
                gene_weight=df_gs[n][1],
                ctrl_match_key="mean_var",
                n_ctrl=200,
                weight_opt="vs",
                return_ctrl_raw_score=False,
                return_ctrl_norm_score=True,
                verbose=False)

df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/scPagwas.magma.PBMC.topgenes.monocytecount.csv", index_col=0)
df_gs = df_gs.loc[["scPagwastop1000","magmatop1000"],:]
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/geneset.gs")


gene_list  = df_gs['scPagwastop1000'][0]
gene_weight  = df_gs['scPagwastop1000'][1]
df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=1000)
df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/mono.pbmc.1000scPagwas.df_model2.csv", sep=",", index=False)


 
gene_list  = df_gs['magmatop1000'][0]
gene_weight  = df_gs['magmatop1000'][1]
df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=1000)
df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Mono1000magma.df_model2.csv", sep=",", index=False)

gene_list  = df_gs['scPagwastop1000'][0]
gene_weight  = df_gs['scPagwastop1000'][1]
df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=1000)
df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Mono1000scPagwas.df_model2.csv", sep=",", index=False)


df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/magma_scRDS_10_1000.csv", index_col=0)
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/geneset.gs")
for i in range(1,101):
 a = 'MAGMAtop' + str(i)
 gene_list  = df_gs[a][0]
 gene_weight  = df_gs[a][1]
 df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=1000)
 df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scDRS/"+a+".scDRS.csv", sep=",", index=False)

```

### Visualization for figure3

```R
library(ggplot2)
library(ggpubr)
library(ggtext)
#################
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0")
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
df_list<-lapply(c("1000scPagwas","1000magma","500scPagwas","500magma"),function(i){
 df<-read.csv(paste0("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Mono",i,".df_model2.csv"))
 return(df$norm_score)
})

df_list<-as.data.frame(df_list)
Pagwas$scDRS.1000scPagwas<-df_list[[1]]
Pagwas$scDRS.1000magma<-df_list[[2]]
Pagwas$scDRS.500scPagwas<-df_list[[3]]
Pagwas$scDRS.500magma<-df_list[[4]]
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


######################
#rank plot
###################

library('ComplexHeatmap')
library(circlize)
####
colors_celltypes=c("#DF7861","#8CC0DE")

all_fortify_can<-all_fortify_can[order(all_fortify_can$scPagwas.TRS.Score1,decreasing=T),]
all_fortify_can$celltype<-as.vector(all_fortify_can$types)
pdf("rankplot_scPagwas.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "scPagwas rank"))
dev.off()

all_fortify_can<-all_fortify_can[order(all_fortify_can$scDRS.1000magma,decreasing=T),]
all_fortify_can$types<-as.vector(all_fortify_can$types)
pdf("rankplot_magma.pdf",height =7,width =2)
print(Heatmap(data.matrix(all_fortify_can$types),
              col =colors_celltypes,
              cluster_columns = F,
              cluster_rows = F,
              border=T,
              row_gap = unit(0.25, "mm"),
              show_parent_dend_line=T,
              name = "MAGMA rank"))
dev.off()

############
#proportion percent plot
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
a <- data.frame(
  group = names(percent1),
  value = percent1)

p1<-ggdonutchart(a, "value", label = rev(a$value),
                 fill = "group", color = "white",
                 palette = c("#DF7861","#8CC0DE") )

####################magma
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
p3<-ggdonutchart(a, "value", label = rev(a$value),
                 fill = "group", color = "white",
                 palette = c("#DF7861","#8CC0DE") )

pdf("model_percent_monocyte.pdf",width =8,height = 8)
ggpubr::ggarrange(p1,p3,nrow = 2,ncol =1)
dev.off()


```

## Rolypoly for single cell

```R
lapply(list.files("/share/pub/dengcy/GWAS_Multiomics/Pkg/rolypoly-master/rolypoly-master/R/"),function(x){
source(paste0("/share/pub/dengcy/GWAS_Multiomics/Pkg/rolypoly-master/rolypoly-master/R/",x))
})

library(Seurat)
library(SeuratObject)
library("dplyr")
library("foreach")
library("ggplot2")
library("glmnet")
library("MASS")
library("Matrix")
library("matrixcalc")
library("doParallel")
library("foreach")
library(plyr)

cl.cores=detectCores(logical=F)
cl=makeCluster(cl.cores-2)
cl=makeCluster(5)
registerDoParallel(cl)

load("/share/pub/dengcy/Singlecell/COVID19/data/covid_ld.RData")
library("data.table")
load("/share/pub/dengcy/Singlecell/COVID19/1.rolypoly_result/geneid_df1.RData")
i<-"monocytecount"
suppressMessages(gwas_data <- bigreadr::fread2(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt")))
 Single_data =readRDS("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds")
 Single_mat<-GetAssayData(Single_data,assay="RNA")

 geneid_df1<-geneid_df1[!duplicated(geneid_df1$label),]
ld_path <- "/share/pub/dengcy/Singlecell/COVID19/data/LD"
rolypoly_result <- rolypoly_roll(
  gwas_data =as.data.frame(gwas_data),
  block_annotation = geneid_df1,
  block_data =as.data.frame(Single_mat[,1:3]) ,
  ld_folder =ld_path,
  bootstrap_iters = 100
)
#a<-rolypoly_result$block_heritability_contribution
save(rolypoly_result,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/rolypoly/rolypoly_monocytecount_modeldata.RData")
```

Integrate rolypoly：

```
rolypoly_score<-unlist(lapply(1:20,function(i){
load(paste0("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/rolypoly_modeldata",i,".RData"))
return(rolypoly_result$block_heritability_contribution)
}))
```

#### Put the rolypoly  result to scPagwas: 

```R
 load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_nk_v1.9.1.RData") 
rolypoly_score<-unlist(lapply(1:20,function(i){
load(paste0("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/rolypoly_modeldata",i,".RData"))
return(rolypoly_result$block_heritability_contribution)
}))

rolypoly_p<-unlist(lapply(1:20,function(i){
load(paste0("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/rolypoly_modeldata",i,".RData"))
return(rolypoly_result$bootstrap_results$bp_value[-1])
}))
Pagwas$rolypoly_score<-rolypoly_score
Pagwas$rolypoly_p<-rolypoly_p
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_nk_v1.9.1.RData")
```

## model1:DCtest, Monocytcount: Compare with other methods

### 

### 2.Integrate result to seruat

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth")
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")

for(i in c(1:10)){
    j<-100*i
scPagwas_topgenes <- names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation, decreasing = T), ])[1:j]
Pagwas <- Seurat::AddModuleScore(Pagwas, assay = "RNA", list(scPagwas_topgenes), name = paste0("scPagwas.",j,"topgenes.Score"))
}

i<-"monocytecount"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_magma_genes.RData"))

for(i in c(1:10)){
    j<-100*i
magmatop1000<-intersect(magma_genes$symbol[1:j],rownames(Pagwas))
Pagwas <- Seurat::AddModuleScore(Pagwas, assay = "RNA", list(magmatop1000), name = paste0("magma.",j,"seruat.Score"))
}
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")



library(scPagwas)
suppressMessages(library(Seurat))
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggtext)

df<-read.csv('/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/mono.1000scPagwas.df_modeldata1.csv')
Pagwas@meta.data$scPagwas_monocytecount_scDRS<-df$norm_score
"raw_score"
df<-read.csv('/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/mono.1000magma.df_modeldata1.csv')
Pagwas@meta.data$magma_monocytecount_scDRS<-df$raw_score


scDRS_re<-read.csv('/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0/ThreeMethods_scDRSresult_model1.csv')
Pagwas@meta.data$smultixcan_scDRS<-scDRS_re$smultixcan_top100_monocytecount
Pagwas@meta.data$spredixcan_scDRS<-scDRS_re$spredixcan_top100_monocytecount
Pagwas@meta.data$TWAS_scDRS<-scDRS_re$TWAS_top100_monocytecount
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")
```

### 3.UMAP plot

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
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")

all_fortify_can <- fortify.Seurat.umap(Pagwas)

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas100), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("scPagwas.topgenes.Score1")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0/Umap.modelmonocytecount_model1.scPagwas.topgenes.Score1.pdf",width = 6, height = 6)
print(p1)
dev.off()

#all_fortify_can$scDRS.1000scPagwas[1]<-9

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas_monocytecount_scDRS), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("scPagwas.1000_scDRS")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0/Umap.modelmonocytecount_model1.scPagwas.1000_scDRS_score.pdf",width = 6, height = 6)
print(p1)
dev.off()


    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scDRS.1000magma), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("magma.1000_scDRS")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0/Umap.modelmonocytecount_model1.magma.1000_scDRS_score.pdf",width = 6, height = 6)
print(p1)
dev.off()

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

save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")

save(all_fortify_can,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0/all_fortify_can.RData")
```

### 4.Rank plot

```R
######################
#rank plot
###################

library('ComplexHeatmap')
library(circlize)
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0")
####
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")
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
```

### 5.Percent plot

```R

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

        p1<-ggdonutchart(df, "value", label = "group",
                         fill = "group", color = "white",
                         palette = c("#DF7861", "#D4E2D4") )
        return(p1)
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

pdf("percent.methods_scdrs.pdf",width = 10,height = 8)
ggpubr::ggarrange(p1,p2,p3,p4,p5,nrow = 2)
dev.off()
```

![image-20221025100348175](D:\OneDrive\GWAS_Multiomics\scPagwas_reproduce\scPagwas_reproduce\image-20221025100348175.png)

### 4.Compare with other score methods

```R
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
library(testSctpa)
library(VISION)
library(AUCell)
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0")
#分别对以下数据进行计算
#load("Pagwas_monocytecount_nk_v1.9.1.RData")
load("Pagwas_monocytecount_modeldata_v1.10.0.RData")
##输出mgt文件
#gwass<-c('eosinophilcount', 'basophilcount','LymphocytePercent','monocytecount' ,'neutrophilcount', 'WhiteBloodCellcount', 'Hemoglobinconcen', 'MeanCorpuscularHemoglobin', 'MeanCorpusVolume')
gwass<-c('monocytecount')
methods<-c('smultixcan', 'spredixcan' ,'TWAS')
for (gwas in gwass) {
	for (method in methods) {
		files<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.csv"))

genes<-read.csv(files)
rn= glue::glue("{method}_top")
b<-paste0(rn,1000)

a<-c(b,"NA",genes$gene_name)
a<-matrix(a,nrow=1)
out<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
write.table(a,file=out,append = T,col.names = F,row.names = F,quote = F,sep = "\t")
	}
}

##计算得分

counts = load_counts()
se_oj = CreateSeuratObject(counts)
se_oj = cal_PAS(seurat_object = se_oj,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = 'log',
              species = 'mouse', 
              pathway='kegg')

DefaultAssay(object = Pagwas) <- "RNA"

methods<-c('smultixcan', 'spredixcan' ,'TWAS')
gwas='monocytecount'
m_l<-lapply(methods,function(method){
  gmt<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
  DefaultAssay(object = Pagwas) <- "RNA"
  auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)

  rn= glue::glue("{method}-top")
  b<-paste0(rn,1000)
  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
  return(auccell_df[,b])
	})
auc_df<-as.data.frame(m_l)
colnames(auc_df)<-paste0(methods,"_monocytecount_AUCell")
Pagwas@meta.data<-cbind(Pagwas@meta.data,auc_df)

###scPagwas magma
auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
gmt_file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scPagwas.magma.Model.monocytecount.gmt")

auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
Pagwas$magma_monocytecount_AUCell<-auccell_df$magmatop1000 
Pagwas$scPagwas_monocytecount_AUCell<-auccell_df$scPagwastop1000
####################
######Vision
##################
m_l<-lapply(methods,function(method){
  gmt<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
  DefaultAssay(object = Pagwas) <- "RNA"

  auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'Vision',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)

  rn= glue::glue("{method}-top")
  b<-paste0(rn,1000)
  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
  return(auccell_df[,b])
	})
df<-as.data.frame(m_l)
colnames(df)<-paste0(methods,"_monocytecount_Vision")
Pagwas@meta.data<-cbind(Pagwas@meta.data,df)

###scPagwas magma
auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'Vision',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
gmt_file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scPagwas.magma.Model.monocytecount.gmt")

auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
Pagwas$magma_monocytecount_Vision<-auccell_df$magmatop1000 
Pagwas$scPagwas_monocytecount_Vision<-auccell_df$scPagwastop1000
save(Pagwas,file="Pagwas_monocytecount_modeldata_v1.10.0.RData")

#1.umap plot
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0")
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


library('ComplexHeatmap')
library(circlize)
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0")
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

        p1<-ggdonutchart(df, "value", label = "group",
                         fill = "group", color = "white",
                         palette = c("#DF7861", "#D4E2D4") )
        return(p1)
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


pdf("percent.score.methods_scPagwas.pdf",width = 10,height = 8)
ggpubr::ggarrange(p1,p2,p3,p4,nrow = 2)
dev.off()

save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")
```

## model1: Lymphocyte

```R
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

setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount")
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Pagwas_Lymphocytecount_modeldata_v1.10.0.RData")


scDRS_re<-read.csv('/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_outputv1.10.0/ThreeMethods_scDRSresult_model1.csv')
Pagwas@meta.data$smultixcan_scDRS_Lymphocytecount3<-scDRS_re$smultixcan_top100_Lymphocytecount3
Pagwas@meta.data$spredixcan_scDRS_Lymphocytecount3<-scDRS_re$spredixcan_top100_Lymphocytecount3
Pagwas@meta.data$TWAS_scDRS_Lymphocytecount3<-scDRS_re$TWAS_top100_Lymphocytecount3
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Pagwas_Lymphocytecount_modeldata_v1.10.0.RData")
```

### scDRS

```R

import scdrs
from scipy import stats
import pandas as pd
import scanpy as sc
from anndata import AnnData
sc.set_figure_params(dpi=125)
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings

warnings.filterwarnings("ignore")
# load adata
DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/"
H5AD_FILE = os.path.join(DATA_PATH, "Pagwas_modelgroundtruth_addata.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
#sc.tl.pca(adata, svd_solver="arpack")
#sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
scdrs.preprocess(adata)

df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scPagwas.Model.magmatopgenes_scRDSLymphocytecount3.csv", index_col=0)
df_gs = df_gs.loc[["scPagwastop1000","magmatop1000"],:]
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/geneset.gs")

gene_list  = df_gs['magmatop1000'][0]
gene_weight  = df_gs['magmatop1000'][1]
df_score = scdrs.score_cell(
                data=adata,
                gene_list=gene_list,
                gene_weight=gene_weight,
                ctrl_match_key="mean_var",
                n_ctrl=200,
                weight_opt="vs",
                return_ctrl_raw_score=False,
                return_ctrl_norm_score=True,
                verbose=False)

df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount31000magma.df_model.csv", sep=",")

gene_list  = df_gs['scPagwastop1000'][0]
gene_weight  = df_gs['scPagwastop1000'][1]
df_score = scdrs.score_cell(
                data=adata,
                gene_list=gene_list,
                gene_weight=gene_weight,
                ctrl_match_key="mean_var",
                n_ctrl=200,
                weight_opt="vs",
                return_ctrl_raw_score=False,
                return_ctrl_norm_score=True,
                verbose=False)
df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount31000scPagwas.df_model.csv", sep=",")
```

### Umap

```R
scDRS_re<-read.csv('/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount31000magma.df_model.csv')
Pagwas$magma_scDRS_Lymphocytecount3.score <- scDRS_re$norm_score

scDRS_re<-read.csv('/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount31000scPagwas.df_model.csv')
Pagwas$scPagwas_scDRS_Lymphocytecount3.score <- scDRS_re$norm_score

load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/Lymphocytecount3_Hema_bmmc_scPagwas_v1.10.0.RData")

genes<-names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation, decreasing = T), ])[1:1000]

load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Pagwas_Lymphocytecount_modeldata_v1.10.0.RData")

Pagwas <- Seurat::AddModuleScore(Pagwas, assay = "RNA", list(genes), name = "scPagwas_seurat_score")

save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Pagwas_Lymphocytecount_modeldata_v1.10.0.RData")

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

```



### score methods

```R
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
library(testSctpa)
library(VISION)
library(AUCell)


counts = load_counts()
se_oj = CreateSeuratObject(counts)
se_oj = cal_PAS(seurat_object = se_oj,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = 'log',
              species = 'mouse', 
              pathway='kegg')

DefaultAssay(object = Pagwas) <- "RNA"

#methods<-c('smultixcan', 'spredixcan' ,'TWAS')
#gwass<-c('Lymphocytecount3')

  gmt<-'/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scPagwas.magma.Model.Lymphocytecount3.gmt'
  DefaultAssay(object = Pagwas) <- "RNA"
  auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)
  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))

Pagwas@meta.data$scPagwas_Lymphocytecount3_AUCell<-auccell_df$scPagwastop1000

  auccell = cal_PAS(seurat_object = Pagwas,
                       tool = 'Vision',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)
  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))

Pagwas@meta.data$scPagwas_Lymphocytecount3_Vision<-auccell_df$scPagwastop1000


all_fortify_can <- fortify.Seurat.umap(Pagwas)

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas_Lymphocytecount3_Vision), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("scPagwas_Lymphocytecount3_Vision")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Umap.scPagwas_Lymphocytecount3_Vision_score.pdf",width = 6, height = 6)
print(p1)
dev.off()

    p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas_Lymphocytecount3_AUCell), size = 0.3, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#F7F5F2",high="#FD5D5D")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("scPagwas_Lymphocytecount3_AUCell")

pdf(file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Umap.scPagwas_Lymphocytecount3_AUCell_score.pdf",width = 6, height = 6)
print(p1)
dev.off()

save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Pagwas_Lymphocytecount_modeldata_v1.10.0.RData")
```

### Percent

```R

all_fortify_can$type<-"Lymphocyte"
all_fortify_can$type[all_fortify_can$celltype=="monocytes"]<-"non_Lymphocyte"
all_fortify_can$type[all_fortify_can$celltype=="DC"]<-"non_Lymphocyte"

percent_plot<-function(df1){
        n<-nrow(df1)
        b<-rep("top",n)
        b[(0.55*n+1):n]<-"bottom"

        df1$rg<-b
        print(table(data.frame(df1$type,df1$rg)))
        t1<-table(data.frame(df1$type,df1$rg))
        percent1<-t1[,2]/sum(t1[,2])
print(percent1)
        return(p1)
}

df1<-all_fortify_can[order(all_fortify_can$scPagwas_scDRS_Lymphocytecount3.score,decreasing=T),]
p1<-percent_plot(df1)
#Lymphocyte non_Lymphocyte
#     0.8181818      0.1818182

df1<-all_fortify_can[order(all_fortify_can$scPagwas_Lymphocytecount3_AUCell,decreasing=T),]
p1<-percent_plot(df1)
# Lymphocyte non_Lymphocyte
#     0.8181818      0.1818182

df1<-all_fortify_can[order(all_fortify_can$scPagwas_Lymphocytecount3_Vision,decreasing=T),]
p1<-percent_plot(df1)
#   Lymphocyte non_Lymphocyte
#     0.8181818      0.1818182

df1<-all_fortify_can[order(all_fortify_can$scPagwas_seurat_score1,decreasing=T),]
p1<-percent_plot(df1)
# Lymphocyte non_Lymphocyte
#     0.8181818      0.1818182


df1<-all_fortify_can[order(all_fortify_can$magma_scDRS_Lymphocytecount3.score,decreasing=T),]
p2<-percent_plot(df1)
# Lymphocyte non_Lymphocyte
#     0.4927273      0.5072727

df1<-all_fortify_can[order(all_fortify_can$smultixcan_scDRS,decreasing=T),]
p3<-percent_plot(df1)
# Lymphocyte non_Lymphocyte
#     0.2809091      0.7190909



df1<-all_fortify_can[order(all_fortify_can$spredixcan_scDRS,decreasing=T),]
p4<-percent_plot(df1)
# Lymphocyte non_Lymphocyte
#     0.1672727      0.8327273

df1<-all_fortify_can[order(all_fortify_can$TWAS_scDRS,decreasing=T),]
p5<-percent_plot(df1)
#  Lymphocyte non_Lymphocyte
#     0.1072727      0.8927273


pdf("percent.methods_scdrs.pdf",width = 10,height = 8)
ggpubr::ggarrange(p1,p2,p3,p4,p5,nrow = 2)
dev.off()
```



## model2:DC

```
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
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_nk_v1.9.1.RData")
all_fortify_can <- fortify.Seurat.umap(Pagwas)

percent_plot<-function(df1){
        n<-nrow(df1)
        b<-rep("top",n)
        b[(0.5*n+1):n]<-"bottom"

        df1$rg<-b
        print(table(data.frame(df1$type,df1$rg)))
        t1<-table(data.frame(df1$type,df1$rg))
        percent1<-t1[,2]/sum(t1[,2])
print(percent1)
        #return(p1)
}

df1<-all_fortify_can[order(all_fortify_can$scPagwas.topgenes.Score1,decreasing=T),]
p1<-percent_plot(df1)

all_fortify_can$rolypoly_p<- -log2(all_fortify_can$rolypoly_p)

df1<-all_fortify_can[order(all_fortify_can$rolypoly_p,decreasing=T),]
p1<-percent_plot(df1)
# 0     1
#0.504 0.496

```



## sub-functions

```
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

```

