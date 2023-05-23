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

