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
                     output.dirs="modeldata_output",
                     block_annotation = block_annotation,
                     assay="RNA",
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     seruat_return=T,
                     celltype=F,
                     ncores = 1)
save(Pagwas,file="Pagwas_monocytecount_modeldata.RData")
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
```

```python
################################scdrs
import scdrs
from scipy import stats
import pandas as pd
import seruat as sc
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

### 2.Integrate result to seruat
### 4.Compare with other score methods

```R
########### S-PrediXcan, S-MutiXcan
#/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits

#library(org.Hs.eg.db)
#library(dplyr)
#a<-data.frame(ensembl_id=rownames(exprSet5))
#g2s=toTable(org.Hs.egSYMBOL)
#g2e=toTable(org.Hs.egENSEMBL)
#s2e=merge(g2s,g2e,by="gene_id",all.x=T)
#save(s2e,file="/share/pub/dengcy/refGenome/symbol2emsembl.RData")
#install.packages("WebGestaltR")
Args <- commandArgs(T)

gwas = "Lymphocytecount3"#print(Args[1])
method ='spredixcan' #print(Args[2])
os.file= file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/"))
scDRS.file= file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.scDRS.multiMethods.genes.csv"))
enrich.file= file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.csv"))

setwd(os.file)

if(method =="smultixcan"){

resultfile ='/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/Lymphocytecount3/test/smultixcan/eqtl/COVID_GWAS_round_4_mashr_smultixcan_eqtl.txt'

  resultfile = file.path(glue::glue("/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/{gwas}/{method}/eqtl/COVID_GWAS_round_4_mashr_smultixcan_eqtl.txt"))
  rdf<-read.table(resultfile,header=T)
  rdf<-rdf[,c("gene_name","pvalue")]

}

if(method =="spredixcan"){
resultfile ='/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/Lymphocytecount3/test/spredixcan/eqtl/mashr/COVID_round_4th_GWAS_mashr_Whole_Blood.db.csv'

  resultfile = file.path(glue::glue("/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/{gwas}/{method}/eqtl/mashr/COVID_round_4th_GWAS_mashr_Whole_Blood.db.csv"))
  rdf<-read.csv(resultfile)
  rdf<-rdf[,c("gene_name","pvalue")]

}

if(method =="TWAS"){

	library(stringr) 
	load("/share/pub/dengcy/refGenome/symbol2emsembl.RData")
	colnames(s2e)<-c("gene_id","gene_name","ID")
  resultfile = file.path(glue::glue("/share2/pub/chenchg/chenchg/TWAS/result/{gwas}_result/{gwas}_merge.dat"))
  rdf<-read.table(resultfile,header=T)
  rdf<-rdf[,c("ID","TWAS.P")]
  rdf$ID<-unlist(str_sub(rdf$ID,1,15))
  rdf=merge(rdf,s2e,by="ID",all.x=T)
  rdf=rdf[,c("gene_name","TWAS.P")]
  colnames(rdf)<-c("gene_name","pvalue")
}

rdf<-rdf[order(rdf$pvalue,decreasing=F),]

ab<-lapply(1:100, function(i){
  b<-paste(rdf$gene_name[1:(10*i)],collapse=",")
  return(b)
})

rn= glue::glue("{method}_top")
b<-paste0(rn,1:100)
a<-data.frame(genes=unlist(ab))
rownames(a)<-b
write.csv(a,file=scDRS.file)

rdf<-rdf[1:1000,]
write.csv(rdf,file=enrich.file)
```

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
setwd("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_output")

load("Pagwas_monocytecount_modeldata.RData")

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
save(Pagwas,file="Pagwas_monocytecount_modeldata.RData")
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
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Pagwas_Lymphocytecount_modeldata.RData")


scDRS_re<-read.csv('/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/modeldata_output/ThreeMethods_scDRSresult_model1.csv')
Pagwas@meta.data$smultixcan_scDRS_Lymphocytecount3<-scDRS_re$smultixcan_top100_Lymphocytecount3
Pagwas@meta.data$spredixcan_scDRS_Lymphocytecount3<-scDRS_re$spredixcan_top100_Lymphocytecount3
Pagwas@meta.data$TWAS_scDRS_Lymphocytecount3<-scDRS_re$TWAS_top100_Lymphocytecount3
save(Pagwas,file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Lymphocytecount/Pagwas_Lymphocytecount_modeldata.RData")
```

### scDRS

```R
import scdrs
from scipy import stats
import pandas as pd
import seruat as sc
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

```