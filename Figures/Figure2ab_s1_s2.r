### 3.Gene ranked plot

###Figure2ab

library(ggplot2)
library(reshape2)
library(ggpubr)
traits<-c("eosinophilcount","basophilcount","LymphocytePercent","Lymphocytecount3",
             "monocytecount","neutrophilcount","WhiteBloodCellcount","Hemoglobinconcen",
             "MeanCorpuscularHemoglobin","MeanCorpusVolume")

lapply(traits,function(i){
  print(i)
  n_l<-lapply(c("scPagwas","magma","TWAS","spredixcan","smultixcan"),function(j){
    
    if(j=="magma"){
      if(i=="MeanCorpuscularHemoglobin"){
        goresult<-NULL
        a1<-a2<-a3<-a4<-0
      }else{
        goresult<-read.delim(paste0("D:/OneDrive/GWAS_Multiomics/Compare/goanalysis/",i,"magmatop1000/enrichment_results_wg_result.txt"),header = T)
        a1<-sum(goresult$FDR<0.01)
        a2<-sum(goresult$FDR<0.001)
        a3<-sum(goresult$FDR<0.0001)
        a4<-sum(goresult$FDR<0.00001)
      }
      
    }else if(j=="scPagwas"){
      goresult<-read.delim(paste0("D:/OneDrive/GWAS_Multiomics/Compare/goanalysis/",i,"scPagwastop1000/enrichment_results_wg_result.txt"),header = T)
      a1<-sum(goresult$FDR<0.01)
      a2<-sum(goresult$FDR<0.001)
      a3<-sum(goresult$FDR<0.0001)
      a4<-sum(goresult$FDR<0.00001)
    }else{
      file= file.path(glue::glue("D:/OneDrive/GWAS_Multiomics/Compare/goanalysis/{j}/Project_{i}/enrichment_results_{i}.txt"))
      
      if( file.exists(file)){
        goresult<-read.delim(file,header = T)
        a1<-sum(goresult$FDR<0.01)
        a2<-sum(goresult$FDR<0.001)
        a3<-sum(goresult$FDR<0.0001)
        a4<-sum(goresult$FDR<0.00001) 
      }else{
        a1<-a2<-a3<-a4<-0
      }
    }
    return(c(a1,a2,a3,a4))
  })
  

  go_df <- data.frame(scPagwas=n_l[[1]],
                      magma=n_l[[2]],
                      TWAS=n_l[[3]],
                      spredixcan=n_l[[4]],
                      smultixcan=n_l[[5]],
                      FDR=c("1e-02","1e-03","1e-04","1e-05"))
  ##plot
  gg_go_df<-melt(go_df,id.vars = "FDR")
  gg_go_df$FDR<-factor(gg_go_df$FDR,levels = c("1e-05","1e-04","1e-03","1e-02"))
  
  setwd("D:/OneDrive/GWAS_Multiomics/Compare/goanalysis/")
  pdf(paste0("dorplot_",i,"_goranalysis_allgenes.pdf"),height =3,width=5)
  print(ggdotchart(gg_go_df, x="FDR", y="value", color = "variable",          
             palette = c('#e8505b','#f9d56e','#f3ecc2','#14b1ab','#d4b5b0'), #
             sorting = "none",    # 
             size =1,dot.size=4,
             rotate = T,label="value",
             add = "segments", #
             main=i, ylab="Numbers of significant GO terms",
             ggtheme = theme_pubr()) )
  
  dev.off()
  })
  

lapply(traits,function(i){
 if(i=="MeanCorpuscularHemoglobin"){
    magma_goresult<-NULL
    scPagwas_goresult<-read.delim(paste0("E:/OneDrive/GWAS_Multiomics/Compare/goanalysis/",i,"scPagwastop1000/enrichment_results_wg_result.txt"),header = T)
    write.csv(scPagwas_goresult,file=paste0("E:/OneDrive/GWAS_Multiomics/Manuscripts/Supplementfiles/top1000goanalysis/",i,"scPagwas_top1000genes_goresult.csv"))

  }else{
    magma_goresult<-read.delim(paste0("E:/OneDrive/GWAS_Multiomics/Compare/goanalysis/",i,"magmatop1000/enrichment_results_wg_result.txt"),header = T)
    scPagwas_goresult<-read.delim(paste0("E:/OneDrive/GWAS_Multiomics/Compare/goanalysis/",i,"scPagwastop1000/enrichment_results_wg_result.txt"),header = T)
    write.csv(magma_goresult,file=paste0("E:/OneDrive/GWAS_Multiomics/Manuscripts/Supplementfiles/top1000goanalysis/",i,"magma_top1000genes_goresult.csv"))
    write.csv(scPagwas_goresult,file=paste0("E:/OneDrive/GWAS_Multiomics/Manuscripts/Supplementfiles/top1000goanalysis/",i,"scPagwas_top1000genes_goresult.csv"))

}
})



#### Enrichment the ranked gene


library(scPagwas)
suppressMessages(library(Seurat))
traits<-c("eosinophilcount","basophilcount","LymphocytePercent","Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

setwd("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/")
for(i in traits){
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/scDRS/scPagwas.magmatopgenes_",i,".RData"))
names(topgene)<-c("scPagwastop1000","magmatop1000","scPagwastop500","magmatop500")
}

load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
geneList<-Pagwas@misc$gene_heritability_correlation[,1]
geneList=sort(geneList,decreasing = T) 
select_cells<-c("01_HSC","02_Early.Eryth","05_CMP.LMPP","09_pDC","12_CD14.Mono.2","19_CD8.N")
set.seed(1234)
a_list<-lapply(select_cells,function(x){
  a<-subset(Pagwas,idents =x) 
  a<-a[,sample(1:ncol(a),500)]
    return(a)
})
Pagwas_rank<-merge(x = a_list[[1]],y =  a_list[2:length(select_cells)])
save(Pagwas_rank,file="Pagwas_rank.RData")
data_mat<- GetAssayData(Pagwas_rank,assay ="RNA")
##################
a<-apply(data_mat,1,mean)
a<-sort(a,decreasing=T)                
rank_gene<-names(a)
save(rank_gene,file="expr_rank_gene.RData")
g1<-intersect(rank_gene[1:(1/2*length(a))],topgene$scPagwastop1000)
g2<-intersect(rank_gene[(1/2*length(a)):length(a)],topgene$scPagwastop1000)

m1<-intersect(rank_gene[1:(1/2*length(a))],topgene$magmatop1000)
m2<-intersect(rank_gene[(1/2*length(a)):length(a)],topgene$magmatop1000)

gene_list<-list(g1,g2,m1,m2)
names(gene_list)<-c("g1","g2","m1","m2")
lapply(names(gene_list), function(x){
  write.csv(gene_list[[x]],file=paste0(x,"split2_gene.csv"),row.names = F)
})
write.csv(degene,file="de_gene_monocyte.csv")
```

#### Ranked gene visualize Figure2b,FigureS1

setwd("E:/OneDrive/GWAS_Multiomics/Compare/5.16compareresult")
load("D:/OneDrive/GWAS_Multiomics/Compare/goanalysis/expr_rank_gene.RData")
load("D:/OneDrive/GWAS_Multiomics/Compare/goanalysis/Lymphocytecount3_genes_split2.RData")

percent1<-c(length(c(gene_list$g1))/length(unlist(gene_list[1:2])),
            length(c(gene_list$g2))/length(unlist(gene_list[1:2]))
)
percent2<-c(length(c(gene_list$m1))/length(unlist(gene_list[3:4])),
            length(c(gene_list$m2))/length(unlist(gene_list[3:4]))
)

b1<-rep(0,length(rank_gene))
b1[rank_gene %in% unlist(gene_list[1:2])]<-1
b2<-rep(0,length(rank_gene))
b2[rank_gene %in% unlist(gene_list[3:4])]<-1

a<-c(rep(1,length(rank_gene)*0.5),rep(2,length(rank_gene)*0.5))
rank_df<-data.frame(rank_gene,a,b1,b2)

library(ggpubr)
df <- data.frame(
  group = c("scPagwas","others"),
  value = c(percent1[1],1-percent1[1]))

p1<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#FFE194", "#4C4C6D") )
df <- data.frame(
  group = c("scPagwas","others"),
  value = c(percent1[2],1-percent1[2]))

p2<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#FFE194", "#4C4C6D")  )  
df <- data.frame(
  group = c("scPagwas","others"),
  value = c(percent2[1],1-percent2[1]))
p3<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#FFE194", "#4C4C6D")  )  

df <- data.frame(
  group = c("scPagwas","others"),
  value = c(percent2[2],1-percent2[2]))
p4<-ggdonutchart(df, "value", label = "group",
                 fill = "group", color = "white",
                 palette = c("#FFE194", "#4C4C6D")  )  
setwd("D:/OneDrive/GWAS_Multiomics/Compare")

pdf("percent.for.toprankgenes.lymphocyte.pdf",width = 10,height = 10)
ggpubr::ggarrange(p1,p2,p3,p4,nrow = 2,ncol = 4)
dev.off()

setwd("D:/OneDrive/GWAS_Multiomics/Compare")
library('ComplexHeatmap')
library(circlize)
col_fun = colorRamp2(c(0, 1), c("#F6F6F6", "#161D6F"))
col_fun(seq(0, 1))

pdf("lymphocyte.heatmap_rank_PAGWASgene.pdf",width = 2)
Heatmap(data.matrix(rank_df$b1),
        col = col_fun,
        #left_annotation = ha, 
        cluster_columns = F,
        cluster_rows = F,
        color_space="HLS",
        border=T,
        row_gap = unit(0.25, "mm"),
        show_parent_dend_line=T,
        #name = "-log2(p)",
        #row_order =order(rdf$types),
        show_row_names=T,
        show_column_names=T
        #row_split=rdf$phenotypes
)
dev.off()

pdf("lymphocyte.heatmap_rank_MAGMAgene.pdf",width = 2)
Heatmap(data.matrix(rank_df$b2),
        col = col_fun,
        #left_annotation = ha, 
        cluster_columns = F,
        cluster_rows = F,
        color_space="HLS",
        border=T,
        row_gap = unit(0.25, "mm"),
        show_parent_dend_line=T,
        #name = "-log2(p)",
        #row_order =order(rdf$types),
        show_row_names=T,
        show_column_names=T
        #row_split=rdf$phenotypes
)
dev.off()

###FigureS2
library(ggplot2)
library(ggpubr)
library(patchwork)
library(readr)
library(dplyr)
library(tidyverse)

GO_overlap <- readr::read_csv("./gather_overlap.csv")
pp <- ggdotchart(
  data = GO_overlap,
  x='Phenotype',
  y='Number of significant GO-BP terms',
  add = "segments",
  color = 'Methods',
  palette =c('#F6D569','#E64C60','#D2B4AF','#33ACAA','#F2EAC4'),
  group = 'Methods',
  label = 'Number of significant GO-BP terms',
  dot.size = 3,
  sorting = "ascending"
) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0))
ggsave(filename = './overlap.pdf',plot = pp,width = 14,height = 7)
```