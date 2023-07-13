
library(reshape2)
library(ggplot2)
library(paletteer)
#pal <- paletteer_d("ggsci::nrc_npg")[c(1,3,4,9,5,2,6,8,10)]
pal <- paletteer_d("ggsci::nrc_npg")[c(1,2,4,5,6)]
#pal <- paletteer_d("ggsci::nrc_npg")[c(1,5,10,2,6)]

setwd("F:\\02-温医大-科研数据\\0000-scRNA\\000-单细胞+GWAS-软件开发\\001-scRNA-GWAS-软件开发\\00-Data-analysis\\01-比较GWAS-COVID-19数据集")

#data <- read.table("severe.txt", sep="\t", header=T)
data <- read.table("severe_non_scDRS.txt", sep="\t", header=T)
data2 <- melt(data)
#利用levels对数据进行排序
data2$annotation <- factor(data2$annotation,levels = data$annotation)


#data3 <- read.table("./moderate.txt", sep="\t", header=T)
data3 <- read.table("./moderate_non_scDRS.txt", sep="\t", header=T)
data4 <- melt(data3)
#利用levels对数据进行排序
data4$annotation <- factor(data4$annotation,levels = data3$annotation)


#data5 <- read.table("./mild.txt", sep="\t", header=T)
data5 <- read.table("./mild_non_scDRS.txt", sep="\t", header=T)
data6 <- melt(data5)
#利用levels对数据进行排序
data6$annotation <- factor(data6$annotation,levels = data5$annotation)


#data7 <- read.table("./normal.txt", sep="\t", header=T)
data7 <- read.table("./normal_non_scDRS.txt", sep="\t", header=T)

data8 <- melt(data7)
#利用levels对数据进行排序
data8$annotation <- factor(data8$annotation,levels = data7$annotation)

opar1 <- par(no.readonly=TRUE)
par(mfrow=c(2,2))


ggplot(data2, aes(x=annotation, y=value))+
  geom_point(aes(color=variable), size=3)+
  scale_color_manual(values=pal)+
  scale_y_sqrt("-Log10(P-value)", breaks=c(0.5,2,5,10,15))+
  theme_classic()+
  theme(axis.text.x =  element_text(angle = 90, hjust = 1, vjust = .5))+
  #geom_hline(yintercept =-log10(0.05), color="red", lty=12)
  geom_hline(yintercept =-log10(0.003846154), color="red", lty=11)



ggplot(data4, aes(x=annotation, y=value))+
  geom_point(aes(color=variable), size=3)+
  scale_color_manual(values=pal)+
  scale_y_sqrt("-Log10(P-value)", breaks=c(0.5,2,5,10,15))+
  theme_classic()+
  theme(axis.text.x =  element_text(angle = 90, hjust = 1, vjust = .5))+
  #geom_hline(yintercept =-log10(0.05), color="red", lty=12)
  geom_hline(yintercept =-log10(0.003846154), color="red", lty=11)


ggplot(data6, aes(x=annotation, y=value))+
  geom_point(aes(color=variable), size=3)+
  scale_color_manual(values=pal)+
  scale_y_sqrt("-Log10(P-value)", breaks=c(0.5,2,5,10,15))+
  theme_classic()+
  theme(axis.text.x =  element_text(angle = 90, hjust = 1, vjust = .5))+
  #geom_hline(yintercept =-log10(0.05), color="red", lty=12)
  geom_hline(yintercept =-log10(0.003846154), color="red", lty=11)


ggplot(data8, aes(x=annotation, y=value))+
  geom_point(aes(color=variable), size=3)+
  scale_color_manual(values=pal)+
  scale_y_sqrt("-Log10(P-value)", breaks=c(0.5,2,5,10,15))+
  theme_classic()+
  theme(axis.text.x =  element_text(angle = 90, hjust = 1, vjust = .5))+
  #geom_hline(yintercept =-log10(0.05), color="red", lty=12)
  geom_hline(yintercept =-log10(0.003846154), color="red", lty=11)

par(opar1)
