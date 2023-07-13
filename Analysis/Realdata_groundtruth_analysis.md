# Monocytecount Realdata_groundtruth
the real data groudtruth data are produce from BMMC and PBMC data, we mainly used the BMMC groundtruth data.
## 一， BMMC

### 1.Select the example celltypes for compare

```R
library(ggpubr)
library(ggplot2)
library(ggtext)

suppressMessages(library(Seurat))

##################
#GET EXAMPLE
###############
setwd("D:/OneDrive/GWAS_Multiomics/realgroundtruth")
load("D:/OneDrive/GWAS_Multiomics/Compare/Hema_test2/monocytecount_bmmc_Pagwas1.10.0.RData")
table(Pagwas$celltypes)
meta.data<-Pagwas@meta.data
celltypes<-c("11_CD14.Mono.1","12_CD14.Mono.2","13_CD16.Mono","09_pDC","10_cDC",
             "17_B","19_CD8.N","20_CD4.N1","25_NK")

#colors_celltypes=c("#F38BA0","#FFBCBC","#EDF6E5","#B5EAEA")
df<-meta.data[meta.data$celltypes %in% celltypes,]
df$types<-"monocyte"
df$types[df$celltypes=="09_pDC"]<-"non_monocyte"
df$types[df$celltypes=="10_cDC"]<-"non_monocyte"
df$types[df$celltypes=="19_CD8.N"]<-"non_monocyte"
df$types[df$celltypes=="20_CD4.N1"]<-"non_monocyte"
df$types[df$celltypes=="25_NK"]<-"non_monocyte"
df$types[df$celltypes=="17_B"]<-"non_monocyte"

groundtruth_samples<-c(rownames(df)[sample(which(df$celltypes=="09_pDC"),100)],
                       rownames(df)[sample(which(df$celltypes=="10_cDC"),100)],
                       rownames(df)[sample(which(df$celltypes=="19_CD8.N"),1500)],
                       rownames(df)[sample(which(df$celltypes=="17_B"),1000)],
                       rownames(df)[sample(which(df$celltypes=="20_CD4.N1"),1500)],
                       rownames(df)[sample(which(df$celltypes=="25_NK"),800)],
  rownames(df)[sample(which(df$types=="monocyte"),5000)])

df<-df[groundtruth_samples,]
df$sig<-1
df$sig[df$CellScaleqValue>0.01]<-0
table(df[,c("types","sig")])

Pagwas_groundtruth<-Pagwas[,groundtruth_samples]

Pagwas_groundtruth <- FindVariableFeatures(Pagwas_groundtruth,nfeatures = 3000)
Pagwas_groundtruth <- RunPCA(object = Pagwas_groundtruth, assay = "RNA", npcs = 50)

##subfunction:
cluster_pca_umap <- function(obj,assay=NULL, reduction,cluster_res = 0.3){
  #obj2 <- RunPCA(obj, assay = "SCT", reduction = "harmony",verbose = F)
  obj2 <- RunTSNE(object = obj,assay = assay, reduction = reduction, dims = 1:50,check_duplicates = FALSE)
  obj2 <- RunUMAP(object = obj2, assay =assay, reduction = reduction, dims = 1:50,check_duplicates = FALSE)
  obj2 <- FindNeighbors(object=obj2, assay = assay, reduction = reduction, dims = 1:50)
  obj2 <- FindClusters(object=obj2, resolution = cluster_res)
  return(obj2)
}
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

load("Pagwas_groundtruth.monocyte.RData")
Pagwas_groundtruth<-cluster_pca_umap(obj = Pagwas_groundtruth,reduction="pca",cluster_res = 0.3)
save(Pagwas_groundtruth,file = "Pagwas_groundtruth.monocyte.RData")


color_scanpy_patient <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")
pdf("groundtruth_monocyte_test8.13.pdf",width = 6,height = 6)
 DimPlot(Pagwas_groundtruth,reduction ="umap",
                         group.by = "celltypes",pt.size=0.2,
                         label = F, repel=TRUE)+ 
  umap_theme()+labs(x="TSNE",y="")+
  scale_colour_manual(name = "celltypes", values =color_scanpy_patient[c(2:9,1)]) +
  theme(aspect.ratio=1)
dev.off()

```



### 4.Trait associated genes for 10 to 1000

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/monocytecount_bmmc_Pagwas1.10.0.RData")
load("Pagwas_groundtruth.monocyte.RData")

#table(Pagwas$celltypes)
meta.data<-Pagwas@meta.data

oa<-order(Pagwas@misc$gene_heritability_correlation,decreasing = T)

topgene<-lapply(1:100, function(i){
  g1<-rownames(Pagwas@misc$gene_heritability_correlation)[oa[1:(10*i)]]
})
names(topgene)<-paste0(1:100,"times_topgene")
Pagwas_groundtruth <- Seurat::AddModuleScore(Pagwas_groundtruth, 
                                 assay = "RNA", 
                                 topgene,name="scPagwas")

i<-"monocytecount"

load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_magma_genes.RData"))

#load("Pagwas_groundtruth.monocyte.RData")

#Pagwas_groundtruth<-Pagwas[,colnames(Pagwas_groundtruth)]
Pagwas_groundtruth$types<-"monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="09_pDC"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="10_cDC"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="19_CD8.N"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="20_CD4.N1"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="25_NK"]<-"non_monocyte"
Pagwas_groundtruth$types[Pagwas_groundtruth$celltypes=="17_B"]<-"non_monocyte"

percent1<-lapply(paste0("scPagwas",1:100), function(i){
  a<-Pagwas_groundtruth@meta.data$types[order(Pagwas_groundtruth@meta.data[,i],decreasing = T)[1:5000]]
  return(sum(a=="monocyte")/5000)
})

#####scDRS

ab<-lapply(1:100, function(i){
  b<-paste(magma_genes$symbol[1:(10*i)],collapse=",")
  return(b)
})

b<-paste0("MAGMAtop",1:100)
a<-data.frame(genes=unlist(ab))
rownames(a)<-b
write.csv(a,file="magma_scRDS_10_1000.csv")

########——————————————————————————————————————————————
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("Pagwas_groundtruth.monocyte.RData")
library(SeuratDisk)
library(Seurat)
a<-GetAssayData(Pagwas_groundtruth,assay = "RNA",slot="counts")
Pagwas_readgroundtruth <- CreateSeuratObject(counts = a,
                           project = "scPagwas",
                           min.cells = 3,
                           min.features = 200)


DefaultAssay(Pagwas_readgroundtruth) <- "RNA"
SaveH5Seurat(Pagwas_readgroundtruth, "Pagwas_groundtruth_addata.h5seurat")
Convert("Pagwas_groundtruth_addata.h5seurat", dest="h5ad")
```

### scDRS for magma and scPagwas

```python

DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/"
H5AD_FILE = os.path.join(DATA_PATH, "Pagwas_groundtruth_addata.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
#sc.tl.pca(adata, svd_solver="arpack")
#sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
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

df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/mono.bmmc.1000magma.df_model.csv", sep=",")

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

df_score.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/mono.bmmc.1000scPagwas.df_model.csv", sep=",")
```


### 8.Compare with other score methods

```R
files<-"/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/mono.1000magma.df_modeldata1.csv"

magma_genes<-read.table("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/monocytecountmagmatop1000.txt",header=T)

a<-c("MAGMA1000","NA",magma_genes$x)
a<-matrix(a,nrow=1)
out<-"/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/monocytecountmagmatop1000.genes.mgt"
write.table(a,file=out,append = T,col.names = F,row.names = F,quote = F,sep = "\t")


  auccell = cal_PAS(seurat_object = Pagwas_groundtruth,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/monocytecountmagmatop1000.genes.mgt")

  auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
  Pagwas_groundtruth$magma_monocytecount_AUCell<-auccell_df$MAGMA1000 

df1<-all_fortify_can[order(all_fortify_can$magma_monocytecount_AUCell,decreasing=T),]
p<-percent_plot(df1)

save(Pagwas_groundtruth,file="Pagwas_groundtruth_bmmc_monocyte_v10.RData")
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
setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")

load('Pagwas_groundtruth_bmmc_monocyte_v10.RData')

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
  DefaultAssay(object = Pagwas_groundtruth) <- "RNA"
  auccell = cal_PAS(seurat_object = Pagwas_groundtruth,
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
Pagwas_groundtruth@meta.data<-cbind(Pagwas_groundtruth@meta.data,auc_df)

###scPagwas magma
auccell = cal_PAS(seurat_object = Pagwas_groundtruth,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
gmt_file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scPagwas.magma.Model.monocytecount.gmt")

auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
Pagwas_groundtruth$magma_monocytecount_AUCell<-auccell_df$magmatop1000 
Pagwas_groundtruth$scPagwas_monocytecount_AUCell<-auccell_df$scPagwastop1000

####################
######Vision
##################
m_l<-lapply(methods,function(method){
  gmt<-file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.mgt"))
  DefaultAssay(object = Pagwas_groundtruth) <- "RNA"

  auccell = cal_PAS(seurat_object = Pagwas_groundtruth,
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
Pagwas_groundtruth@meta.data<-cbind(Pagwas_groundtruth@meta.data,df)

###scPagwas magma
auccell = cal_PAS(seurat_object = Pagwas_groundtruth,
                       tool = 'Vision',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
gmt_file="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/scPagwas.magma.Model.monocytecount.gmt")

auccell_df<- as.data.frame(t(GetAssayData(auccell, slot="data", assay="PAS")))
Pagwas_groundtruth$magma_monocytecount_Vision<-auccell_df$magmatop1000 
Pagwas_groundtruth$scPagwas_monocytecount_Vision<-auccell_df$scPagwastop1000
save(Pagwas_groundtruth,file="Pagwas_groundtruth_bmmc_monocyte_v10.RData")
```



## 二、PBMC

### 1.Run Pagwas and get top1000 genes

```R
library(scPagwas)
library(ggplot2)
suppressMessages(library(Seurat))
suppressMessages(library("dplyr"))
i<-"monocytecount"
Pagwas_pbmc_groundtruth<-scPagwas_main(Pagwas = NULL,
                                       gwas_data=paste0("/share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/",i,"_prune_gwas_data.txt"),
                                       Single_data ="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_pbmc_realgroundtruth.rds",
                                       output.prefix="real_pbmc_monocytecount",
                                       output.dirs="real_pbmc_monocytecount",
                                       block_annotation = block_annotation,
                                       assay="RNA",
                                       Pathway_list=Genes_by_pathway_kegg,
                                       chrom_ld = chrom_ld,
                                       singlecell=T,
                                       seruat_return=T,
                                       celltype=F,
                                       ncores = 1)
save(Pagwas_pbmc_groundtruth,file = "/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_pbmc_groundtruth_monocytecount.RData")

load("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_pbmc_groundtruth_monocytecount.RData")
i<-"monocytecount"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_magma_genes.RData"))

scPagwas_genes<-names(Pagwas_pbmc_groundtruth@misc$gene_heritability_correlation[order(Pagwas_pbmc_groundtruth@misc$gene_heritability_correlation,decreasing=T),])

magmatop1000<-intersect(magma_genes$symbol,rownames(Pagwas_pbmc_groundtruth))[1:1000]
scPagwastop1000<-scPagwas_genes[1:1000]
magmatop500<-intersect(magma_genes$symbol,rownames(Pagwas_pbmc_groundtruth))[1:500]
scPagwastop500<-scPagwas_genes[1:500]
topgene<-list(scPagwastop1000,magmatop1000,scPagwastop500,magmatop500)
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")

a<-data.frame(genes=c(paste(scPagwastop1000,collapse=","),paste(magmatop1000,collapse=","),paste(scPagwastop500,collapse=","),paste(magmatop500,collapse=",")))
rownames(a)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
write.csv(a,file=paste0("scPagwas.magma.PBMC.topgenes.",i,".csv"))

```



### 3.Get top genes for 10 to 1000

```R
i<-"monocytecount"
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/",i,"_magma_genes.RData"))

ab<-lapply(1:100, function(i){
  b<-paste(magma_genes$symbol[1:(10*i)],collapse=",")
  return(b)
})
b<-paste0("MAGMAtop",1:100)
a<-data.frame(genes=unlist(ab))
rownames(a)<-b
write.csv(a,file="/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/PBMC_monocyte_magma_scRDS_10_1000.csv")


setwd("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth")
load("Pagwas_pbmc_groundtruth_monocytecount.RData")
library(SeuratDisk)
library(Seurat)
a<-GetAssayData(Pagwas_pbmc_groundtruth,assay = "RNA",slot="counts")
Pagwas_pbmc_groundtruth <- CreateSeuratObject(counts = a,
                                             project = "scPagwas",
                                             min.cells = 3,
                                             min.features = 200)
DefaultAssay(Pagwas_pbmc_groundtruth) <- "RNA"
SaveH5Seurat(Pagwas_pbmc_groundtruth, "Pagwas_pbmc_groundtruth.h5seurat")
Convert("Pagwas_pbmc_groundtruth.h5seurat", dest="h5ad") 

```

### 4.Run scDRS for top 10-1000 genes

```python
##scDRS
# load adata
DATA_PATH = "/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/"
H5AD_FILE = os.path.join(DATA_PATH, "Pagwas_pbmc_groundtruth.h5ad")
adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
scdrs.preprocess(adata)

df_gs=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/PBMC_monocyte_magma_scRDS_10_1000.csv", index_col=0)
df_gs= df_gs.reset_index().to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/geneset.gs", sep="\t", index=False)
df_gs = scdrs.util.load_gs("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/geneset.gs")

for i in range(1,101):
 a = 'MAGMAtop' + str(i)
 gene_list  = df_gs[a][0]
 gene_weight  = df_gs[a][1]
 df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=1000)
 df_res.to_csv("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/PBMC/"+a+".scDRS.csv", sep=",", index=False)
#######r
```
