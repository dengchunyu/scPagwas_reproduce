
library(reshape2)
library(ggplot2)
library(paletteer)
#pal <- paletteer_d("ggsci::nrc_npg")[c(1,3,4,9,5,2,6,8,10)]
pal <- paletteer_d("ggsci::nrc_npg")[c(1,2,4,5,6)]
#pal <- paletteer_d("ggsci::nrc_npg")[c(1,5,10,2,6)]


#########-----AD-scPagwas-benchmarking analysis
data22 <- read.table("AD_dotplot.txt", sep="\t", header=T)

data23 <- melt(data22)
#利用levels对数据进行排序
data23$annotation <- factor(data23$annotation,levels = data22$annotation)


ggplot(data23, aes(x=annotation, y=value))+
  geom_point(aes(color=variable), size=3)+
  scale_color_manual(values=pal)+
  scale_y_sqrt("-Log10(P-value)", breaks=c(0.5,2,5,10,15))+
  theme_classic()+
  theme(axis.text.x =  element_text(angle = 90, hjust = 1, vjust = .5))+
  #geom_hline(yintercept =-log10(0.05), color="red", lty=12)
  geom_hline(yintercept =-log10(0.05), color="red", lty=11)

