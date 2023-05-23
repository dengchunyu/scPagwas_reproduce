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
