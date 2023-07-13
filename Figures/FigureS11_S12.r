####################1000 test 
#Figure S12
setwd("/share/pub/dengcy/GWAS_Multiomics/test/1000test")
for(i in c("p0.001_a0.72","p0.001_a_0","p0.01_a0.72","p0.01_a_0","p0.1_a0.72","p0.1_a_0")){
  #导入raw result data
  re_ll<-list()
  for(i in c(2000,4000,6000,8000,10000)){
  load("liver_raw_results_",i,".RData")
  re_l<-list()
  for(j in c(50,100,150,200,250,300)){
    load("liver_perm_results_",i,"_",j,".RData")
    cor_reult<-lapply(1:length(perm_results),function(x){
      cor(bootstrap_results$bt_value,perm_results$bt_value,method="spearman")
    })
    cor_reult<-as.data.frame(cor_reult)
    cor_reult$pathway_number<-j
    re_l[[j]]<-cor_reult
  }
  re_l<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),re_l)
re_l$nfeatures<-i
re_ll[[i]]<-re_l
  }
  re_ll<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),re_ll)
p<- ggplot(re_l,aes(x=pathway_number,y=cor_reult,fill=nfeatures))+
  geom_boxplot(outlier.size=0.1,width= 0.2,position=position_dodge(0.7),alpha=0.6)+ #箱式图异常值大小调整
  geom_jitter(color="black" ,size=0.5,alpha=0.3)+theme(legend.position="none")+ #不需要图例
  theme_bw()+scale_fill_manual(values=c("#D1E8E4","#C37B89","#BCCC9A","#FFC898","#CFB784",
  "#C6D57E","#79B4B7","#F6A9A9","#3E7C17","#F4A442","#1597E5","#9CC094")) #填充颜色调整
#save plot
pdf(paste0("liver_",i,".pdf"))
print(p)
dev.off()
}
###################
#barplot figure s11
######################
for(i in c("p0.001_a0.72","p0.001_a_0","p0.01_a0.72","p0.01_a_0","p0.1_a0.72","p0.1_a_0")){
  load(paste0("liver_perm_results_",i,".RData"))

perm_results$adjp <- rep(0,nrow(perm_results))
perm_results$adjp[perm_results$bp_value<0.01] <- 1

gg<-perm_results %>% 
  group_by(annotation) %>% 
  summarise(significant=sum(adjp)) # 计算列指定列名
gg<-as.data.frame(gg)
##柱状图
library(ggplot2) 
library(ggthemes) 

p1<-ggplot(gg, aes(x=annotation, y=significant,fill=annotation))+ 
ylim(0, 1000)+  
geom_bar(position = "dodge",stat = "identity")+ theme_bw()+   
labs(x = "", y = "Number of significant result", title = "1000 test")+   
geom_hline(aes(yintercept=50), colour="#990000", linetype="dashed")+
coord_flip()
pdf(paste0("liver_",i,".pdf"))
print(p1)
dev.off()
}
