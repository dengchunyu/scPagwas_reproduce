#################
#PC test for scPagwas
library(scPagwas)
 library(reshape2)
 library(ggplot2)
library(ggpubr)    # as_ggplot
library(cowplot)   # draw_plot_label
library(ggtext)
library(psych)
#install.packages("/share/pub/dengcy/software/scPagwas_1.10.4.tar.gz",repos=NULL,type="source")

setwd("/share/pub/dengcy/GWAS_Multiomics/compare")
library(scPagwas)
suppressMessages(library(Seurat))
 #Input pathway gene list, you can construct with youself.

Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")
Pagwas <- Single_data_input(Pagwas=NULL,
                                assay="RNA",
                                Single_data=Single_data,
                                Pathway_list=Genes_by_pathway_kegg)
Pathway_list=Genes_by_pathway_kegg
Pa.len <- unlist(lapply(Pathway_list, function(Pa) length(Pa)))

  # Valid the total genes of pathway.
  pa_gene <- unique(unlist(Pathway_list))
  pana <- names(Pathway_list)[which(
    unlist(lapply(
      Pathway_list,
      function(Pa) {
        length(intersect(
          Pa,
          Pagwas$VariableFeatures
        ))
      }
    )) > 2
  )]
  Pathway_list <- Pathway_list[pana]
  # keep the raw pathway
  Pagwas$rawPathway_list <- Pathway_list
  celltypes <- as.vector(unique(Pagwas$Celltype_anno$annotation))
  pana_list <- lapply(celltypes, function(celltype) {
    scCounts <- Pagwas$data_mat[
      ,
      Pagwas$Celltype_anno$cellnames[
        Pagwas$Celltype_anno$annotation == celltype
      ]
    ]
    scCounts <- as_matrix(scCounts)
    scCounts <- scCounts[rowSums(scCounts) != 0, ]
    proper.gene.names <- rownames(scCounts)
    pana <- names(Pathway_list)[which(
      unlist(lapply(
        Pathway_list,
        function(Pa) {
          length(intersect(
            Pa, proper.gene.names
          ))
        }
      )) > 2
    )]
    return(pana)
  })

  Pagwas$Pathway_list <- Pathway_list[Reduce(intersect, pana_list)]
  Pagwas$rawPathway_list <- Pathway_list[Reduce(intersect, pana_list)]

 scPCAscore_list <- lapply(celltypes, function(celltype) {
    scCounts <- Pagwas$data_mat[
      ,
      Pagwas$Celltype_anno$cellnames[
        Pagwas$Celltype_anno$annotation == celltype
      ]
    ]
    scCounts <- as_matrix(scCounts)
    #

      Pathway_list = Pagwas$Pathway_list
        nPcs <- 1
  scCounts <- t(scCounts)

  cm <- Matrix::colMeans(scCounts)
  proper.gene.names <- colnames(scCounts)

  ###### calculate the pca for each pathway terms.
  #pdf(paste0(celltype,".fa.parallel.plot.pdf"))
  prl <- lapply(Pathway_list, function(Pa_id) {
    lab <- proper.gene.names %in% intersect(proper.gene.names, Pa_id)
       result <- tryCatch(
      {
         p1 <- irlba::prcomp_irlba( scCounts[, lab], n=10)
         ps<-summary(p1)
         pr<-ps[7][[1]][2,]/ps[7][[1]][3,10]
         return(pr)
      },
      error = function(e) {
        return(NULL)
      }
    )
   
    #mat<-na.omit(scCounts[, lab])
    #mat<-mat[,colSums(mat)!=0]
    #print(
    #fa.parallel(mat,fa='pc',n.iter = 100,show.legend = F,ncomp=10,main = 'Scree plot with parallel analysis'))
   #pc <- principal(mat,nfactors = 1)
   
    #pc <- principal(mat,nfactors = 1,scores = T)

    
  })
 prl<-prl[!unlist(lapply(prl,is.null))]
 prd<-as.data.frame(t(as.data.frame(prl)))
 prd2<-apply(prd,2,function(x){
 	a<-summary(x)
    x[x>a[5]]<-a[5]
    return(x)
 	})
prd2<-as.data.frame(prd2)
 prd2$pathway<-rownames(prd)

 gg_prd<-melt(prd2,id.var="pathway")
p<-ggplot(gg_prd, aes(x = variable, y =value,fill=variable)) + 
 geom_boxplot(outlier.size = 0.5,alpha=0.5)+theme_classic() +
 geom_jitter(color="black",size=0.5,alpha=0.15)+
 theme(legend.position="none")+
 labs(title = celltype)
pdf(paste0(celltype,".varProportion.plot.pdf"))
print(p)
dev.off()

})

##################
#
 celltype<-"03_Late.Eryth"
 scCounts <- Pagwas$data_mat[
      ,
      Pagwas$Celltype_anno$cellnames[
        Pagwas$Celltype_anno$annotation == celltype
      ]
    ]
 scCounts <- as_matrix(scCounts)
    #

  Pathway_list = Pagwas$Pathway_list
        nPcs <- 1
  scCounts <- t(scCounts)

  cm <- Matrix::colMeans(scCounts)
  proper.gene.names <- colnames(scCounts)

  ###### calculate the pca for each pathway terms.
  pdf(paste0(celltype,".fa.parallel.plot.pdf"))
  lapply(Pathway_list[1:50], function(Pa_id) {
    lab <- proper.gene.names %in% intersect(proper.gene.names, Pa_id)
       result <- tryCatch(
      {
         #p1 <- irlba::prcomp_irlba( scCounts[, lab], n=10)
         mat<-scCounts[, lab]
         mat<-mat[,colSums(mat)!=0]
    print(
    fa.parallel(mat,fa='pc',n.iter = 100,show.legend = F,main = 'Scree plot with parallel analysis'))

         #ps<-summary(p1)
         #pr<-ps[7][[1]][2,]/ps[7][[1]][3,10]
    
      },
      error = function(e) {
        return(NULL)
      })
       })
    
dev.off()