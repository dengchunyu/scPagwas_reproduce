
#=====================================================================================
#
#  Code chunk 3   (2022-05-16-for AD bulk) GSE109887
#
#=====================================================================================

#boxplot_data_TOP50.txt
data <- read.table("boxplot_data_TOP50.txt", header = TRUE)

p<-ggplot(data,aes(x= Pheno,y=GSK3B,fill= Pheno))+geom_violin(trim=FALSE) 
p +geom_boxplot(width=0.1,fill="white")+ theme_classic()+scale_fill_manual(values=c("grey","#f3715c"))

p<-ggplot(data,aes(x= Pheno,y=CREB1,fill= Pheno))+geom_violin(trim=FALSE) 
p +geom_boxplot(width=0.1,fill="white")+ theme_classic()+scale_fill_manual(values=c("grey","#f3715c"))


#=====================================================================================
#
#  Code chunk 4   (2022-05-16-for AD bulk) GSE15222
#
#=====================================================================================

#boxplot_data_TOP50.txt
data <- read.table("boxplot_data_TOP50.txt", header = TRUE)

p<-ggplot(data,aes(x= Pheno,y=GSK3B,fill= Pheno))+geom_violin(trim=FALSE) 
p +geom_boxplot(width=0.1,fill="white")+ theme_classic()+scale_fill_manual(values=c("grey","#f3715c"))

p<-ggplot(data,aes(x= Pheno,y=CREB1,fill= Pheno))+geom_violin(trim=FALSE) 
p +geom_boxplot(width=0.1,fill="white")+ theme_classic()+scale_fill_manual(values=c("grey","#f3715c"))




