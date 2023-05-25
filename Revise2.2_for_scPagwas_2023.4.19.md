# review2 的方法比较问题以及通路问题的解决方案

## 2.1 如何证明我们的通路结果是更具有功能性的

### 2.1.1基于kegg的数据从所有基因表达谱中随机构建通路数据，
进行计算


```R
library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add")
#Genes_by_pathway_kegg
modeldata<-readRDS("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds")
#获得整体基因
genes_all<-rownames(modeldata)
set.seed(1234)
#根据Genes_by_pathway_kegg的通路中基因的数量，随机抽取基因，获得100个随机通路数据
Genes_by_pathway_kegg_random<-lapply(1:100,function(i){
    pl<-lapply(Genes_by_pathway_kegg,function(j){
        sample(genes_all,length(j))
    })
return(pl)
})
save(Genes_by_pathway_kegg_random,file="Genes_by_pathway_kegg_random.RData")
```

运行模拟数据计算scPagwas结果

```R
library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result")
#读取随机抽取的通路数据
load("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Genes_by_pathway_kegg_random.RData")
for(i in 1:100){
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_random_kegg", 
                     output.prefix=paste0("kegg_",i),
                     iters_singlecell = 10,
                     block_annotation = block_annotation,
                     Pathway_list=Genes_by_pathway_kegg_random[[i]],
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file=paste0("model_random_kegg_",i,".RData"))
}

```

结果的准确性计算，计算top50%的基因中单核细胞的比例，以及计算top1000基因显著富集的通路的数量

```R
library(scPagwas)
library(Seurat)
library(WebGestaltR)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result")
percent_list<-c()
refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
pa_num<-list()
outputDirectory <- getwd()
for(i in 1:100){
    load(paste0("model_random_kegg_",i,".RData"))
    df<-Pagwas@meta.data
    #将df按照scPagwas.TRS.Score1的列从大到小排序，获得celltype列中前50%的数据中monocytes元素的占比
    df<-df[order(df$scPagwas.TRS.Score1,decreasing = T),]
    df<-df[1:round(nrow(df)/2),]
    #计算monocytes在celltypes列中的占比
    df<-table(df$celltype)
    percent_list[i]<-df["monocytes"]/sum(df)
    scPagwas_topgenes <- names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation, decreasing = T), ])[1:1000]
    result<-WebGestaltR(enrichMethod="ORA", organism="hsapiens",
        enrichDatabase="geneontology_Biological_Process_noRedundant", 
        interestGene=scPagwas_topgenes,
        interestGeneType="genesymbol", 
        referenceGeneFile=refFile,
        referenceGeneType="genesymbol", 
        isOutput=FALSE,
        projectName=paste0("random_",i))
        a1<-sum(result$FDR<0.01)
        a2<-sum(result$FDR<0.001)
        a3<-sum(result$FDR<0.0001)
        a4<-sum(result$FDR<0.00001)
        pa_num[[i]]<-c(a1,a2,a3,a4)
}
save(percent_list,file="random_kegg_percent_list.RData")
save(pa_num,file="random_kegg_pa_num.RData")

```
### 2.1.2 结果可视化

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result")
load("random_kegg_percent_list.RData")
load("random_kegg_pa_num.RData")
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(ggstatsplot)
##画个相关性点图，percent结果越准确，通路结果越多，相关性比较显著
num_l<-unlist(lapply(pa_num,function(i) i[1]))
df<-data.frame(percent=percent_list,pa_num=num_l)
#画相关性点图，图中标注相关性系数和显著性水平
#在xy轴旁边加入分布图
library(ggstatsplot)
p<-ggscatterstats(df, 
                x = percent, y = pa_num,
               type = "pearson",
               centrality.para = "mean",    
               margins = "both",
               xfill = "#009E73",
               yfill = "#D55E00",
               xlab = "Percent of monocytes in top 50% cells (Sensitivity)",
          ylab = "Number of significant pathways for top 1000 genes")

pdf("random_kegg_percent.pdf",width=8,height=8)
p
dev.off()

png("random_kegg_percent.png",width=1700,height=1700,units = "px",res = 300)
p
dev.off()

#计算random_kegg中每个大的基因list中包含的所有基因list显著富集的通路结果
load("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Genes_by_pathway_kegg_random.RData")
library(WebGestaltR)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result")
#pa_num2<-list()
refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")

for(i in 3:30){
  x<-Genes_by_pathway_kegg_random[[i]]
  a1_list<- lapply(x,function(y){
    #获取报错信息
    tryCatch({
    result<-WebGestaltR(enrichMethod="ORA", organism="hsapiens",
        enrichDatabase="geneontology_Biological_Process_noRedundant", 
        interestGene=y,
        interestGeneType="genesymbol", 
        referenceGeneFile=refFile,
        referenceGeneType="genesymbol", 
        isOutput=FALSE,
        projectName=paste0("random_",i))
        a1<-sum(result$FDR<0.01)
        return(a1)
    },error=function(e) return(0))
        
  })
  panum<-unlist(a1_list)
  save(panum,file=paste0("random_kegg",i,"_enrichment.RData"))
}
```

```shell
#!/bin/bash
#SBATCH -e enrich.err
#SBATCH -o enrich.out
#SBATCH -J enrich
#SBATCH -w in007
#SBATCH --mem=10000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate R4.2
cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result
Rscript 2.r 3

```

上面的随机方式是基于KEGG中所有基因作为基因池，进行分析的，这样不涉及其他基因。但是结果仍然有较大差异，是什么原因导致的呢？这里我们查看所有结果中过滤后的kegg通路的数量比较其差别。


### 2.1.2 将通路中随机比例从0.1到1，每次增加0.1，计算结果

上面的随机方式结果不够平滑，因此我们将通路中加入的随机基因比例从0.05到0.1,0.2,计算结果
这次随机和之前随机不同的是，加入的新基因是去掉KEGG基因之后的其他基因，而不是所有基因
计算结果，观察到底是功能性基因更好还是随机的基因更好?

```R
modeldata <- readRDS("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds")
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result2")
library(scPagwas)
#获得0.1到1的数据list，每次增加0.1
percent_list<-c(0.05,0.1,0.15,0.2,0.25)
genes_all<-rownames(modeldata)
genes_other <- genes_all[!(genes_all %in% unique(unlist(Genes_by_pathway_kegg)))]
percent_list<-c(0.15,0.25)
for(i in percent_list){
Genes_by_pathway_kegg_random<-lapply(1:20,function(i){
    pl<-lapply(Genes_by_pathway_kegg,function(j){
        #从剩余的基因中随机抽取i*length(j)个基因
        num <- floor(i*length(j))
        other <- sample(genes_other,num)
        j<-c(j,other)
    })
return(pl)
})
save(Genes_by_pathway_kegg_random,file=paste0("Genes_by_pathway_kegg_random_",i,".RData"))
}

```

### 2.1.3 运行模拟数据计算上面随机kegg的scPagwas结果

```R
library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result2")
Args <- commandArgs(T)
#gwas='finngen_r7_c3_breast_exallc'
i = print(Args[1])
percent_list<-c(0.05,0.1,0.15,0.2,0.25)
#for(i in percent_list){
i<-0.15
    load(paste0("Genes_by_pathway_kegg_random_",i,".RData"))
    TRS_score<-list()
    gPas_score<-list()
    for(j in 1:20){
        Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs=paste0("random_kegg_",i,"_",j),
output.prefix="",
            block_annotation = block_annotation,
            Pathway_list=Genes_by_pathway_kegg_random[[j]],
            chrom_ld = chrom_ld,
            iters_singlecell = 10,
            singlecell=T,
            celltype=F
        )
        TRS_score[j]<-Pagwas@meta.data$scPagwas.TRS.Score1
        gPas_score[j]<-Pagwas@meta.data$scPagwas.gPAS.score
    }
    save(TRS_score,gPas_score,file=paste0("model_random_kegg_TRS_",i,".RData"))
#}
```


```shell
#!/bin/bash
#SBATCH -e model9.err
#SBATCH -o model9.out
#SBATCH -J model9
#SBATCH -w in007
#SBATCH --mem=50000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate R4.2
cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result2
Rscript 1.r 0.1

```

### 2.1.4 结果的可视化

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result2")
library(scPagwas)
library(Seurat)
library(WebGestaltR)
refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
outputDirectory <- getwd()
#percent_list<-seq(0.1,1,0.1)
Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds"
scdata<-readRDS(Single_data)
anno<-scdata$celltype
percent_list<-c(0.05,0.1,0.15,0.2)
random_kegg_sensitivity <- lapply(percent_list,function(i){
    load(paste0("model_random_kegg_TRS_",i,".RData"))
    #将df按照scPagwas.TRS.Score1的列从大到小排序，获得celltype列中前50%的数据中monocytes元素的占比
    pcl<-lapply(TRS_score,function(x){
      df<-data.frame(TRS=x,celltype=anno)
    df<-df[order(df$TRS,decreasing = T),]
    df<-df[1:round(nrow(df)/2),]
    #计算monocytes在celltypes列中的占比
    df<-table(df$celltype)
    pc<-df["monocytes"]/sum(df)
    return(pc)
    })
return(unlist(pcl))
})
random_kegg_sensitivity<-as.data.frame(do.call(rbind,random_kegg_sensitivity))
random_kegg_sensitivity<-t(as.matrix(random_kegg_sensitivity))
random_kegg_sensitivity<-as.data.frame(random_kegg_sensitivity)
colnames(random_kegg_sensitivity)<- paste0("addnoise",percent_list)
random_kegg_sensitivity$RandomTimes<-paste0("times",1:20)
#melt数据
library(reshape2)
random_kegg_sensitivity<-melt(random_kegg_sensitivity,id.vars = "RandomTimes")
random_kegg_sensitivity$noise_percent<-rep(percent_list,20)
save(random_kegg_sensitivity,file="random_kegg_sensitivity.RData")

setwd("D:/OneDrive/GWAS_Multiomics/Manuscripts/Revise_comments/pathway_add/Random_kegg_result2")
load("random_kegg_sensitivity.RData")
#画图，横坐标为加入噪音的比例，纵坐标为敏感性
library(ggplot2)
library(ggpubr)
p<-ggboxplot(random_kegg_sensitivity, x = "variable", y = "value",color="#60281E",add = "jitter",size = 0.5,palette = "jco",ylab = "Sensitivity",xlab = "Noise percent",ggtheme = theme_bw())

pdf("random_kegg_noise_sensitivity.pdf")
p
dev.off()
png("random_kegg_noise_sensitivity.png")
p
dev.off()
```

上面的结果得出，加入噪音不是影响结果的真正原因，我们猜测真正的原因是和这些基因是否能够覆盖大部分的单核细胞的差异表达基因有关系。
经过测试，和高变基因覆盖率没关系，和单核细胞的差异表达基因覆盖率也没有关系。

```R
#计算这些通路的基因对于单细胞高变基因的覆盖率
library(Seurat)
Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds"
scdata<-readRDS(Single_data)
#计算单细胞中monocytes细胞相对于其他细胞的差异表达基因
allmarker<-FindAllMarkers(scdata,logfc.threshold = 0.25,only.pos = T, min.pct = 0.25,thresh.use = 0.25)
mono_gene<-rownames(allmarker)[allmarker$cluster=="monocytes"]

setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/pathway_add/Random_kegg_result2")
percent_list<-seq(0.1,1,0.1)
gene_cover<-lapply(percent_list,function(i){
    load(paste0("Genes_by_pathway_kegg_random_",i,".RData"))
    pathway_genes_cover<-c()
  for(j in 1:20){
    #获得当前随机下的通路基因
    current_pathway_genes <- unique(unlist(Genes_by_pathway_kegg_random[[j]]))
    #计算当前随机下的通路基因对于单细胞高变基因的覆盖率
    pathway_genes_cover[j]<-length(intersect(current_pathway_genes,mono_gene))/length(mono_gene)
  }
return(pathway_genes_cover)
})
#转化为dataframe
gene_cover<-as.data.frame(do.call(rbind,gene_cover))
gene_cover<-t(as.matrix(gene_cover))
gene_cover<-as.data.frame(gene_cover)
colnames(gene_cover)<- paste0("addnoise",percent_list)
gene_cover$RandomTimes<-paste0("times",1:20)
#melt数据
library(reshape2)
gene_cover<-melt(gene_cover,id.vars = "RandomTimes")
random_kegg_sensitivity$gene_cover<-gene_cover$value
```

## 3. 随机抽取top1000基因进行测试

```R
library(scPagwas)
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")
percent_list<-seq(0.1,1,0.1)
gene_list <-names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation, decreasing = T), ])[1:1000]
#单细胞表达谱
data<-Pagwas@assays$RNA@data
```



## 2.2 模拟构建gwas数据，并计算scPagwas结果

模拟构建gwas数据是为了验证结果是否是假阳性，因为我们的数据是真实的GWAS数据，因此我们需要模拟构建gwas数据，然后计算scPagwas结果，观察结果是否是假阳性。
另外，也需要观察，是否最后产生的基因最后是否也是功能性基因。如果全部是功能基因，说明我们的方法最后确实是基于功能基因的，那么真实数据中找到的功能基因就不具有真实意义，如果不是功能基因，说明在使用真实数据中找到的功能基因具有真实意义。

### 2.2.1 模拟构建gwas数据

```R
# 读入现有的GWAS summary statistics文件，获取有效的SNP列表
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/simulation_gwas/simudata")
gwas_res <- read.table("/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt", header=T)
snp_list <- gwas_res$rsid

# 定义一个函数，用于产生单个随机模拟数据集的GWAS结果
simulate_gwas <- function(snp_list) {
  n_snps <- length(snp_list)

  # 随机生成每个SNP的beta和se估计值
  beta_estimates <- rnorm(n_snps, mean=0, sd=0.1)
  se_estimates <- rnorm(n_snps, mean=0, sd=0.05)

  # 根据beta和se估计值计算每个SNP的t值和p值
  t_values <- beta_estimates / se_estimates
  p_values <- 2 * pt(-abs(t_values), df=nrow(gwas_res)-1)

  # 保存最终结果为GWAS结果格式
  gwas_simulated <- gwas_res
  gwas_simulated$beta <- beta_estimates
  gwas_simulated$se <- se_estimates
    gwas_simulated$p <- p_values
  return(gwas_simulated)
}

# 使用for循环的方式产生100个随机模拟数据集的GWAS结果，并保存为文件
for (i in 1:100) {
  gwas_simulated <- simulate_gwas(snp_list)
  write.table(gwas_simulated, file=paste0("gwas_simulated_",i,".txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}
```

### 2.2.2 运行gwas模拟数据计算scPagwas结果
    
```R
library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/simulation_gwas/simudata")
for(i in 1:100){
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = paste0("gwas_simulated_",i,".txt"),
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                     output.dirs="model_random_gwas", 
                     output.prefix=paste0("gwas_",i),
                     iters_singlecell = 10,
                     block_annotation = block_annotation,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     singlecell=T,
                     celltype=T
)
save(Pagwas,file=paste0("model_random_gwas_",i,".RData"))
}

```

服务器上运行
    
```bash
#!/bin/bash
#SBATCH -e randomgwas.err
#SBATCH -o randomgwas.out
#SBATCH -J randomgwas
#SBATCH -w in006
#SBATCH --mem=70000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate R4.2
Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/simulation_gwas/random_gwas.R
```

### 2.2.3 结果的整合

```R
library(scPagwas)
library(Seurat)
library(WebGestaltR)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/simulation_gwas")
percent_list<-c()
refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
pa_num<-list()
outputDirectory <- getwd()
for(i in 1:100){
    load(paste0("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/simulation_gwas/simudata/model_random_gwas_",i,".RData"))
    df<-Pagwas@meta.data
    #将df按照scPagwas.TRS.Score1的列从大到小排序，获得celltype列中前50%的数据中monocytes元素的占比
    df<-df[order(df$scPagwas.TRS.Score1,decreasing = T),]
    df<-df[1:round(nrow(df)/2),]
    #计算monocytes在celltypes列中的占比
    df<-table(df$celltype)
    percent_list[i]<-df["monocytes"]/sum(df)
    scPagwas_topgenes <- names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation, decreasing = T), ])[1:1000]
    result<-WebGestaltR(enrichMethod="ORA", organism="hsapiens",
        enrichDatabase="geneontology_Biological_Process_noRedundant", 
        interestGene=scPagwas_topgenes,
        interestGeneType="genesymbol", 
        referenceGeneFile=refFile,
        referenceGeneType="genesymbol", 
        isOutput=FALSE,
        projectName=paste0("random_",i))
        a1<-sum(result$FDR<0.01)
        a2<-sum(result$FDR<0.001)
        a3<-sum(result$FDR<0.0001)
        a4<-sum(result$FDR<0.00001)
        pa_num[[i]]<-c(a1,a2,a3,a4)
}
save(percent_list,file="random_gwas_percent_list.RData")
save(pa_num,file="random_gwas_pa_num.RData")


```

### 2.2.4 结果的可视化

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/simulation_gwas")
load("random_gwas_percent_list.RData")
load("random_gwas_pa_num.RData")
library(ggplot2)
library(ggpubr)
pa_num<-do.call(rbind,pa_num)
pa_num<-as.data.frame(pa_num)

data.frame(sensitivity=percent_list,pathway_num=pa_num$V1)->df
p<-ggplot(data = df, aes(x = sensitivity)) + geom_histogram(bins = 10, color = "black", fill = "white") + geom_density(alpha = 0.2, fill = "red") + geom_rug(color = "blue") + theme_bw()
#画数据分布图
pdf("random_gwas_sensitivity_density.pdf")
print(p)
dev.off()
png("random_gwas_sensitivity_density.png")
print(p)
dev.off()
```

## 2.3 对kegg通路增加冗余基因进行准确性测试

### 2.3.1 构建不断增加的冗余性通路

```R
library(scPagwas)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/simulation_kegg")
reduce_percent<-c()
reduce_length<-c()
#读入kegg通路基因
all_genes<-unique(unlist(Genes_by_pathway_kegg))
num<-round(length(unlist(Genes_by_pathway_kegg))/length(Genes_by_pathway_kegg))
len<-length(Genes_by_pathway_kegg)

for(i in c(0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,3.5,4,4.5,5)){
  print(i)
  if(i<1){
    for(j in 1234:1264){
      set.seed(j)
reduce_genes.by.kegg.pathway<-scPagwas::reduce_pathway(
  pathway_seed=names(Genes_by_pathway_kegg)[sample(1:length(Genes_by_pathway_kegg),10)],pathway_list=Genes_by_pathway_kegg,remove_proporion=i)
    }
  }else if(i==1){
    reduce_genes.by.kegg.pathway<-Genes_by_pathway_kegg
  }else{
    #当i>1时，我们将kegg通路中的基因进行扩增i倍
     reduce_genes.by.kegg.pathway<-Genes_by_pathway_kegg
    #从all_genes中不放回的随机选择长度为num的基因
    #选择(i-1)*len次
    times<-round((i-1)*len)
    for(j in 1:times){
      set.seed(j)
      sample_genes<-sample(all_genes,num,replace=F)
      #将sample_genes添加到Genes_by_pathway_kegg的list中
      reduce_genes.by.kegg.pathway<-append(reduce_genes.by.kegg.pathway,list(sample_genes))
    }
    
  }
  print(paste0("the length of reduce pathway is: ",length(reduce_genes.by.kegg.pathway)))
save(reduce_genes.by.kegg.pathway,file=paste0("reduce_genes.by.kegg.pathway_",i,".RData"))
percent<-length(unlist(reduce_genes.by.kegg.pathway))/length(unique(unlist(reduce_genes.by.kegg.pathway)))
print(paste0("the reduce percent is:", percent))
reduce_percent<-c(reduce_percent,percent)
reduce_length<-c(reduce_length,length(reduce_genes.by.kegg.pathway))
}

#输出reduce_percent和reduce_length组成的数据框
reduce_percent_length<-data.frame(reduce_percent,reduce_length)
write.csv(reduce_percent_length,"reduce_percent_length.csv")
```


### 2.3.2 根据上面得到的pathway数据计算模拟和真实数据的scPagwas结果

模拟数据

```R
library(scPagwas)
library(Seurat)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/simulation_kegg")
for(i in c(0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,3.5,4,4.5,5)){
  print(i)
  load(paste0("reduce_genes.by.kegg.pathway_",i,".RData"))
  Pagwas<-scPagwas_main(Pagwas = NULL,
                        gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                        Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds",
                        output.dirs="model_random_kegg", 
                        output.prefix=paste0("kegg_",i),
                        iters_singlecell = 10,
                        block_annotation = block_annotation,
                        Pathway_list=reduce_genes.by.kegg.pathway,
                        chrom_ld = chrom_ld,
                        singlecell=T,
                        celltype=T
  )
  save(Pagwas,file=paste0("Pagwas_kegg_",i,".RData"))
}
```

服务器
```shell
#!/bin/bash
#SBATCH -e rk.err
#SBATCH -o rk.out
#SBATCH -J rk
#SBATCH -w in006
#SBATCH --mem=70000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate R4.2
Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/simulation_kegg/1.r
```

### 2.3.3 结果可视化
  
```R
#计算随着冗余性增加，敏感性的变化
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/simulation_kegg")

sensitivity<-lapply(c(0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,3.5,4,4.5,5),function(i){
  load(paste0("Pagwas_kegg_",i,".RData"))
  df<-Pagwas@meta.data
    #将df按照scPagwas.TRS.Score1的列从大到小排序，获得celltype列中前50%的数据中monocytes元素的占比
    df<-df[order(df$scPagwas.TRS.Score1,decreasing = T),]
    df<-df[1:round(nrow(df)/2),]
    #计算monocytes在celltypes列中的占比
    df<-table(df$celltype)
    return(df["monocytes"]/sum(df))
})
#构建数据框
df<-data.frame(sensitivity=unlist(sensitivity),reduce_percent=c(0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,3.5,4,4.5,5))
library(ggplot)
#绘制点折现图
library(ggplot2)
p<-ggplot(df,aes(x=reduce_percent,y=sensitivity))+geom_point(color="#725E82")+geom_line(color="#725E82")+xlab("reduce_percent")+ylab("sensitivity")+theme_bw()
pdf("sensitivity_reducepathway.pdf")
print(p)
dev.off()
png("sensitivity_reducepathway.png")
print(p)
dev.off()

```

## 2.4 sclinker

```R
library(readr)
SNP<- read_table("/share/pub/dengcy/GWAS_Multiomics/covid19/COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt")
SNP2<-SNP[,c(1:4,7:9,11:13)]
SNP2<-SNP2[,c("CHR","POS","rsid","REF","ALT","all_inv_var_meta_beta","all_inv_var_meta_sebeta","all_inv_var_meta_p","all_meta_AF","all_meta_sample_N")]
colnames(SNP2)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf","N")
colnames(SNP2)<-c("chrom","pos","SNP","A1","A2","BETA","se","pval","maf")
write.table(SNP2,file= "/share/pub/dengcy/GWAS_Multiomics/covid19/COVID19_HGI_gwas_data_ldsc.txt",sep="\t",row.names=F,quote=F)


SNP<- read_table("/share/pub/dengcy/GWAS_Multiomics/ad_test/ieu-b-2.vcf.gz",skip=107)
colnames(SNP)<-c("CHROM","POS","ID","ALT","REF","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
value_list<-lapply(1:nrow(SNP),function(i){
  #index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  #index<-unlist(index)
  value <- strsplit(SNP$IEU[i], split = ":",fixed=T)
  value<-unlist(value)[1:4]
  names(value)<-c("ES","SE","LP","ID")
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","ID")
 value_df<-value_df[,c("ES","SE","LP")]
 gwas_data<-SNP[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 gwas_data$p <-  10^(-gwas_data$LP)
 colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
gwas_data<-gwas_data[,c("CHROM", "POS","ID","ALT","REF","ES","SE","LP")]
colnames(gwas_data)<-c("chrom","pos","SNP","A1","A2","BETA","se","pval")
write.table(gwas_data,file= "/share/pub/dengcy/GWAS_Multiomics/ad_test/AD_gwas_data_ldsc.txt",sep="\t",row.names=F,quote=F)
```

```shell
source activate ldsc
/share/pub/dengcy/software/ldsc-master/munge_sumstats.py \
--sumstats /share/pub/dengcy/GWAS_Multiomics/covid19/COVID19_HGI_gwas_data_ldsc.txt \
--merge-alleles /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/w_hm3.snplist \
--signed-sumstat BETA,0 \
--N  905674 \
--chunksize 500000 \
--out /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI

/share/pub/dengcy/software/ldsc-master/munge_sumstats.py \
--sumstats /share/pub/dengcy/GWAS_Multiomics/ad_test/AD_gwas_data_ldsc.txt \
--merge-alleles /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/w_hm3.snplist \
--signed-sumstat BETA,0 \
--N 63926 \
--chunksize 500000 \
--out /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/AD
```


```R
library(tidyverse)
library("Seurat")

###########################################
library(tidyverse)
library(data.table)
library(dplyr)
gene_coordinates0 <-
  read_tsv("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/NCBI37.3.gene.loc.extendedMHCexcluded",
           col_names = FALSE,col_types = 'cciicc') %>%
  mutate(start=ifelse(X3-50000<0,0,X3-50000),end=X4+50000) %>%
  select(X2,start,end,6,1) %>%
  rename(chr="X2", Gene="X6",EntrezGene="X1") %>%
  mutate(chr=paste0("chr",chr))


#gene_coordinates1 = data.frame(fread("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/ABCpaper_NasserFulcoEngreitz2020_Blood_AvgHiC.txt.gz"))

 df_pre = data.frame(fread("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/AllPredictions.AvgHiC.ABC0.015.minus150.withcolnames.ForABCPaper.txt.gz"))
  df_pre = df_pre[which(df_pre$class == "intergenic" | df_pre$clas == "genic"), ]
  tissuename ="BLD"
  tissuenames2 = as.character(read.table("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/ABC.listbloodQC.txt", header=F)[,1])
 tissue_ids = as.numeric(unlist(sapply(tissuenames2, function(x) return(grep(x, df_pre$CellType)))))
 df = df_pre[tissue_ids, ]
  df2 = cbind.data.frame(df$chr, df$start, df$end, df$TargetGene)
  colnames(df2) = c("chr", "start", "end", "Gene")

  write.table(final_bed, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


Enhancer = read.table("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/Roadmap_Enhancers_Blood.txt",header=F)
df3<-Enhancer[,1:4]
colnames(df3)<-c("chr","start","end","Gene")
df3 = rbind(df3,df2)
gene_coordinates<-df3

#将gene_coordinates0中的EntrezGene列加入到gene_coordinates中
gene_coordinates<-merge(gene_coordinates,gene_coordinates0[,c("Gene","EntrezGene")],by="Gene",all.x=TRUE)
#删除na
gene_coordinates<-gene_coordinates[!is.na(gene_coordinates$EntrezGene),]
save(gene_coordinates,file="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/ABC_gene_coordinates.RData")

###########################################
#计算每个基因的平均表达量，以及输出bed
library(tidyverse)
library(data.table)
library(dplyr)
library(Seurat)
adata1<-readRDS("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds")
Avg_exp_bed(scdata=adata1,celltype="celltype",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/modeldata/")
adata2<-readRDS("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_realgroundtruth.rds")
Avg_exp_bed(scdata=adata2,celltype="celltypes",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/realdata/")

adata1<-readRDS("/share/pub/dengcy/GWAS_Multiomics/ad_test/GSE160936.rds")
Avg_exp_bed(scdata=adata1,celltype="annotation",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/AD/")

#covid19
##moderate
load("/share/pub/dengcy/GWAS_Multiomics/covid19/scpagwas_moderate.v1.10.RData")
Avg_exp_bed(scdata=Pagwas,celltype="annotation",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/moderate/")
##Normal
load("/share/pub/dengcy/GWAS_Multiomics/covid19/scpagwas_Normal.v1.10.RData")
Avg_exp_bed(scdata=Pagwas,celltype="annotation",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/Normal/")
##severe
load("/share/pub/dengcy/GWAS_Multiomics/covid19/scpagwas_severe.v1.91.RData")
Avg_exp_bed(scdata=Pagwas,celltype="annotation",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/severe/")

##mild
load("/share/pub/dengcy/GWAS_Multiomics/covid19/scpagwas_mild.v1.10.RData")
Avg_exp_bed(scdata=Pagwas,celltype="annotation",bedfile="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/mild/")


Avg_exp_bed<-function(scdata,celltype,bedfile,gene_coordinates=gene_coordinates){
average_expression<-AverageExpression(scdata, group.by = celltype)
average_expression<-average_expression$RNA
if(!file.exists(bedfile)){
dir.create(bedfile)
}
ave_func(average_expression,path=bedfile)
}


ave_func<-function(average_expression,path){
average_expression<-na.omit(average_expression)
average_expression<-average_expression[apply(average_expression,1,sum)!=0,]
average_expression<-as.data.frame(average_expression)
top10_function(average_expression,path=path)

}
###top10的函数
top10_function <-function(exp,path){
exp$Gene <- rownames(exp)
exp <- exp %>% add_count(Gene) %>%
  filter(n==1) %>%
  select(-n) %>%
  gather(key = column,value=Expr,-Gene) %>%
  as.tibble()

#############################
###2.Each cell type is scaled to the same total number of molecules.
############################
exp <- exp %>%
  group_by(column) %>%
  mutate(Expr_sum_mean=Expr*1e6/sum(Expr))
##############################
###3.Specificity Calculation
#####################################
exp<- exp %>%
  group_by(Gene) %>%
  mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean)) %>%
  ungroup()
###################
###4.Get MAGMA genes
######################
#Only keep genes that are tested in MAGMA
exp2 <- inner_join(exp,gene_coordinates,by="Gene")
##################
###5.Get number of genes
#################
#Get number of genes that represent 10% of the dataset
n_genes <- length(unique(exp2$EntrezGene))
n_genes_to_keep <- (n_genes * 0.1) %>% round()
##################
###7.Write MAGMA/LDSC input files
#########################
exp2 %>% filter(Expr_sum_mean>1) %>% ldsc_bedfile("column",n_genes_to_keep,path)
print("sucess!")
}

write_group  = function(df,Cell_type,path) {
  df <- select(df,column,chr,start,end,EntrezGene)
  dir.create(paste0(path), showWarnings = FALSE,recursive = TRUE)
  write_tsv(df[-1],paste0(path,make.names(unique(df[1])),".bed"),col_names = F)
  return(df)
}
ldsc_bedfile <- function(d,Cell_type,n_genes_to_keep,path){
  d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity)
  d_spe %>% do(write_group(.,Cell_type,path))
}
```

### 3.5 第五步，计算annot文件

```shell
cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/
source activate mypy
bunzip2 w_hm3.snplist.bz2
tail -n +2 w_hm3.snplist | cut -f 1 > hm_snp.txt

cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/1000G_EUR_Phase3_baseline/
#gunzip baseline.*.annot.gz
for file in baseline.*.annot
do
awk 'NR > 1{print "chr"$1"\t"$2"\t"$2"\t"$3}' $file >> tmp.bed
done
sortBed -i tmp.bed > 1000genomes_phase3_SNPs.bed2
rm tmp.bed

######shell运行
#!/bin/bash
#SBATCH -e Tcells.err
#SBATCH -o Tcells.out
#SBATCH -J Tcells
#SBATCH -w in010
#SBATCH --mem=50000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate ldsc
#B NK DC monocytes Tcells

##########shell脚本
test=$1
f=$2
path_name="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/"
all_snps="1000G_EUR_Phase3_baseline/1000genomes_phase3_SNPs.bed2"
all_annotations="1000G_EUR_Phase3_baseline/"
plink_file="1000G_EUR_Phase3_plink/"
hapmap_snps="hm_snp.txt"
weights="1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
frq="1000G_Phase3_frq/1000G.EUR.QC."
echo $test
mkdir "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/"$test"/results/"
mkdir "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/"$test"/results/"
cd "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/"$test
cd "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/"$test

#for f in *.bed
#do
echo $f
intersectBed -c -a $path_name$all_snps -b $f > $f".1000genomes.intersect"
awk '{if($5!=0) print $4}' $f".1000genomes.intersect" > $f".1000genomes.intersect.snp"
mkdir $f"_tissue_dir"
rm $f".1000genomes.intersect"
cd $f"_tissue_dir"
for j in $path_name$all_annotations/*.annot
do
echo $j
file_name=`basename $j`
perl $path_name/fast_match2_minimal.pl ../$f".1000genomes.intersect.snp" $f $j > $file_name
done
gzip *annot
for i in {1..22}
do
/share/pub/dengcy/software/ldsc-master/ldsc.py --l2 --bfile $path_name$plink_file/1000G.EUR.QC.$i --ld-wind-cm 1 --print-snps $path_name$hapmap_snps --annot baseline.$i.annot.gz --out ./baseline.$i
done
cd ..
rm $f".1000genomes.intersect.snp"
```


### 3.6 第六步，计算ldsc结果

```shell
#!/bin/bash
#SBATCH -e Tcells.err
#SBATCH -o Tcells.out
#SBATCH -J Tcells
#SBATCH -w in010
#SBATCH --mem=50000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate ldsc
sumstats="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/monocytecount_prune.sumstats.gz"
path_name="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/"
weights="1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
frq="1000G_Phase3_frq/1000G.EUR.QC."
all_annotations="1000G_EUR_Phase3_baseline"

cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/modeldata/

for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done
mkdir log_get_pvalues


cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/realdata/

for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done
mkdir log_get_pvalues

sumstats="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/AD.sumstats.gz"
cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/AD/

for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done

sumstats="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI.sumstats.gz"
cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/moderate/

for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done



cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/mild/

for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done

cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/severe/
for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done

cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/Normal/
for f in *_tissue_dir
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
/share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
cd ..
done

```

计算pvalue
```R
library("tidyverse")
library("stringr")
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/modeldata/")
files <- list.files(".",pattern=".tissue_dir.results",full.names = TRUE,recursive=T)
#files <- files[grepl("age",files)]

d <- data_frame(filename=files) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames=gsub(".bed_tissue_dir.results","",basename(filename)),
           makenames=gsub(".bed_continuous_tissue_dir.results","",basename(makenames))) %>% unnest()
  d <- d %>%  filter(Category=="L2_1") %>% mutate(P=1-pnorm(`Coefficient_z-score`)) %>%
  mutate(Trait=sub("_.*","",makenames),Cell_Type=gsub("^_","",str_extract(makenames, "_.*"))) %>%
  select(Trait,Cell_Type,Enrichment,Enrichment_std_error,Enrichment_p,P) %>% arrange(Enrichment_p)
    
write_tsv(d,path="modeldata_cell_types_pvalues.txt")


setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/realdata")
files <- list.files(".",pattern=".tissue_dir.results",full.names = TRUE,recursive=T)
d <- data_frame(filename=files) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames=gsub(".bed_tissue_dir.results","",basename(filename)),
           makenames=gsub(".bed_continuous_tissue_dir.results","",basename(makenames))) %>% unnest()
  d <- d %>%  filter(Category=="L2_1") %>% mutate(P=1-pnorm(`Coefficient_z-score`)) %>%
  mutate(Trait=sub("_.*","",makenames),Cell_Type=gsub("^_","",str_extract(makenames, "_.*"))) %>%
  select(Trait,Cell_Type,Enrichment,Enrichment_std_error,Enrichment_p,P) %>% arrange(Enrichment_p)
write_tsv(d,path="realdata_cell_types_pvalues.txt")


setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/AD")
files <- list.files(".",pattern=".tissue_dir.results",full.names = TRUE,recursive=T)
#files <- files[grepl("age",files)]

d <- data_frame(filename=files) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames=gsub(".bed_tissue_dir.results","",basename(filename)),
           makenames=gsub(".bed_continuous_tissue_dir.results","",basename(makenames))) %>% unnest()
  d <- d %>%  filter(Category=="L2_1") %>% mutate(P=1-pnorm(`Coefficient_z-score`)) %>%
  mutate(Trait=sub("_.*","",makenames),Cell_Type=gsub("^_","",str_extract(makenames, "_.*"))) %>%
  select(Trait,Cell_Type,Enrichment,Enrichment_std_error,Enrichment_p,P) %>% arrange(Enrichment_p)
write_tsv(d,path="AD_cell_types_pvalues.txt")

for(i in c("mild","moderate","severe","Normal")){
setwd(paste0("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/",i))
files <- list.files(".",pattern=".tissue_dir.results",full.names = TRUE,recursive=T)
d <- data_frame(filename=files) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames=gsub(".bed_tissue_dir.results","",basename(filename)),
           makenames=gsub(".bed_continuous_tissue_dir.results","",basename(makenames))) %>% unnest()
  d <- d %>%  filter(Category=="L2_1") %>% mutate(P=1-pnorm(`Coefficient_z-score`)) %>%
  mutate(Trait=sub("_.*","",makenames),Cell_Type=gsub("^_","",str_extract(makenames, "_.*"))) %>%
  select(Trait,Cell_Type,Enrichment,Enrichment_std_error,Enrichment_p,P) %>% arrange(Enrichment_p)
  path=paste0(i,"_cell_types_pvalues.txt")
write_tsv(d,path=path)
}
```
整理结果

```shell
cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/
#for i in Normal severe moderate mild
#cd $i
#ehco $i
for f in *_tissue_dir
do
cd $f
#删除monocytecount_prune为开头的文件
rm -rf monocytecount_prune*
cd ..
done
```




```R
for(i in c("mild","moderate","severe","Normal")){
setwd(paste0("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/COVID19_HGI/",i))
df<-read.table(paste0(i,"_cell_types_pvalues.txt"),header=T)
#删除Cell_Type中的prune_前缀
df$Cell_Type<-gsub("prune_","",df$Cell_Type)
#输出
write.table(df,paste0("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/sclinker_",i,"_pvalues.txt"),sep="\t",quote=F,row.names=F)
}

for(i in c("modeldata","realdata","AD")){
setwd(paste0("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/",i))
df<-read.table(paste0(i,"_cell_types_pvalues.txt"),header=T)
#删除Cell_Type中的prune_前缀
df$Cell_Type<-gsub("prune_","",df$Cell_Type)
#输出
write.table(df,paste0("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/sclinker_",i,"_pvalues.txt"),sep="\t",quote=F,row.names=F)
}

```

## 2.5 EPIC方法计算

```R
#source epic的代码
source("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/package/R/EPIC.R")
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC")
gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt"
library(data.table)
gwas_data<-read.table(gwas_data,header=T)
colnames(gwas_data)<-c("chr", "pos", "rsid", "A1", "A2", "beta", "se", "P", "MAF", "N")
#计算Zscore: z=beta/se
#计算MAC: minor allele count, defined as the number of subjects with at least one observed mutation. MAC=2 * N * MAF.
gwas_data$Zscore <- gwas_data$beta/gwas_data$se
gwas_data$MAC <- 2* gwas_data$N * gwas_data$MAF
snp_to_gene.Demo <- map_snp_to_gene(gwas = gwas_data,anno_path = "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/package/inst/extdata/annotation_dictionary.txt")
snp_division_obj.Demo <- divide_common_rare(gwas = gwas_data, snp_to_gene = snp_to_gene.Demo)
common_snp_to_gene.Demo <- snp_division_obj.Demo$common_snp_to_gene
rare_snp_to_gene.Demo <- snp_division_obj.Demo$rare_snp_to_gene
library(data.table)
prunein_gene_pval.Demo <- read_in(gwas = gwas_data, snp_to_gene = common_snp_to_gene.Demo)
rare_gene_pval.Demo <- read_in(gwas = gwas_data, snp_to_gene = rare_snp_to_gene.Demo)
save.image(file = "model_data.Demo.RData")
save.image(file = "real_data.Demo.RData")
#Single-cell RNA-seq

library(biomaRt)
library(Seurat)
#model data
Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds"
#Single_data <- readRDS(Single_data)
#scRNA.object <- process_scRNA(SeuratObject = Single_data, meta_ct_col = "celltype")
#realdata
load("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_groundtruth_bmmc_monocyte_v10.RData")
#取Single_data和Pagwas_groundtruth基因的交集
genes<-intersect(rownames(Single_data),rownames(Pagwas_groundtruth))
Pagwas_groundtruth<-Pagwas_groundtruth[genes,]
load("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/grch37.RData")
scRNA.object <- process_scRNA(SeuratObject = Pagwas_groundtruth, meta_ct_col = "celltypes",grch37=grch37)
rm(Single_data)

process_scRNA <- function(SeuratObject, meta_ct_col,grch37) {
  gene.loc = read.table(file = "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/package/inst/extdata/gene.noMHC.loc", sep = "\t",
                        header = FALSE, stringsAsFactors = FALSE)

  HGNC_to_ENSE=getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', "strand"),
                     filters = 'hgnc_symbol',
                     values = rownames(SeuratObject),
                     mart = grch37)

  dictionary = as.data.frame(HGNC_to_ENSE)
  dictionary <- dictionary[dictionary$ensembl_gene_id %in% gene.loc$V1, ]

  dictionary$chromosome_name <- as.numeric(dictionary$chromosome_name)
  o <- order(dictionary$chromosome_name, dictionary$start_position, dictionary$end_position)
  dictionary <- dictionary[o,]

  rawCount <- GetAssayData(object = SeuratObject, assay = "RNA", slot = "counts")
  keep.idx <- which(!is.na(match(rownames(rawCount), dictionary$hgnc_symbol)))
  rawCount <- rawCount[keep.idx, ]
  rawCount <- rawCount[match(dictionary$hgnc_symbol, rownames(rawCount)), ]

  # RPKM
  gene.length <- (dictionary$end_position - dictionary$start_position)/1000
  library.size <- apply(rawCount, 2, sum)
  rpkm <- (10^6)*t(t(as.matrix(rawCount))/library.size)/gene.length
  cts <- sort(unique(SeuratObject@meta.data[[meta_ct_col]]))

  log2_rpkm_ct <- matrix(NA, nrow = nrow(rpkm), ncol = length(cts))
  rownames(log2_rpkm_ct) <- rownames(rpkm)
  colnames(log2_rpkm_ct) <- cts
  for (ct in cts) {
    cell.idx = rownames(SeuratObject@meta.data)[which(SeuratObject@meta.data[[meta_ct_col]] == ct)]
    log2_rpkm_ct[,ct] <- log2(1 + apply(rpkm[, cell.idx, drop = FALSE], 1, mean))
  }

  SeuratObject <- subset(SeuratObject, features = rownames(rpkm))

  ngene = 8000
  var.features <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = ngene)

  common.genes.ngene <- var.features@assays$RNA@var.features
  dictionary.ngene <- dictionary[dictionary$hgnc_symbol %in% common.genes.ngene, ]
  log2_rpkm_ct.ngene <- log2_rpkm_ct[rownames(log2_rpkm_ct) %in% common.genes.ngene, ]
  colnames(log2_rpkm_ct.ngene) <- gsub(" ", "_", colnames(log2_rpkm_ct.ngene))

  log2_rpkm_ct.ngene = cbind(log2_rpkm_ct.ngene, apply(log2_rpkm_ct.ngene, 1, mean))
  colnames(log2_rpkm_ct.ngene)[ncol(log2_rpkm_ct.ngene)] = "Average"
  rownames(log2_rpkm_ct.ngene) <- dictionary$ensembl_gene_id[match(rownames(log2_rpkm_ct.ngene), dictionary$hgnc_symbol)]

  log2_rpkm_ct.ngene = as.data.frame(log2_rpkm_ct.ngene)

  dictionary.ngene <- dictionary.ngene[,-1]
  colnames(dictionary.ngene) <- c("V1","V2","V3","V4","V5")

  return(list(scRNA.rpkm = log2_rpkm_ct.ngene, scRNA.loc = dictionary.ngene))
}


scRNA.rpkm <- scRNA.object$scRNA.rpkm
scRNA.loc <- scRNA.object$scRNA.loc
head(scRNA.rpkm)
library(genio)
file_bed <- file.path("/share/pub/dengcy/1000G_EUR_Phase3_plink", "all_chr_1000G.bed")
file_bim <- file.path("/share/pub/dengcy/1000G_EUR_Phase3_plink", "all_chr_1000G.bim")
file_fam <- file.path("/share/pub/dengcy/1000G_EUR_Phase3_plink", "all_chr_1000G.fam")
bim <- read_bim(file_bim)
fam <- read_fam(file_fam)
genotype <- read_bed(file_bed, bim$id, fam$id)

X_ref <- genotype[match(bim$id, rownames(genotype)),]
rm(bim,fam )
save.image(file="real_data.Demo.RData")
save.image(file="model_data.Demo.RData")

inter_dir <-"/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/modeldata/inter"
inter_dir <-"/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/realdata/inter"
save(prunein_gene_pval.Demo,X_ref,scRNA.loc,file="real_poet_sw.RData")
for(chr in 1:22){

calculate_POET_sw(genotype = X_ref, gene_pval = prunein_gene_pval.Demo, 
                  gene.loc = scRNA.loc, # For scRNA-seq data, specify gene.loc = scRNA.loc, which is returned by function process_scRNA()
                  chr =chr, type = "POET", inter_dir = inter_dir)
}
load("real_data.Demo.RData")
pruneinObject.Demo <- get_gene_chisq(gene_pval = prunein_gene_pval.Demo, 
                                     gene.loc = scRNA.loc, # For scRNA-seq data, specify gene.loc = scRNA.loc, which is returned by function process_scRNA() 
                                     type = "POET", inter_dir = inter_dir)
save.image(file = "model_data.Demo.RData")
save.image(file = "real_data.Demo.RData")

prunein_POET_gene_pval.Demo <- pruneinObject.Demo$gene_POET_pval
prunein_gene_corr.mat.Demo <- pruneinObject.Demo$gene_corr.mat
#4.2 Burden test for rare variants
rare_gene_pval.Demo <- get_burden(genotype = X_ref, gene_pval = rare_gene_pval.Demo)

#4.3 Joint analysis for common and rare variants
combine_gene_pval.Demo <- get_combined(prunein_gene_pval = prunein_gene_pval.Demo, rare_gene_pval = rare_gene_pval.Demo)
X_super.Demo <- construct_X_super(genotype = X_ref, rare_gene_pval = rare_gene_pval.Demo)
save.image(file = "model_data.Demo.RData")
save(X_ref,combine_gene_pval.Demo,scRNA.loc,X_super.Demo,file = "model_data.poet2.RData")

load( "model_data.poet2.RData")
for(chr in 1:22){
calculate_POET_sw(genotype = X_ref, gene_pval = combine_gene_pval.Demo, 
                  gene.loc = scRNA.loc, # For scRNA-seq data, specify gene.loc = scRNA.loc, which is returned by function process_scRNA()
                  chr =chr, type = "iPOET", X_super = X_super.Demo, inter_dir = inter_dir)
}
save(X_ref,combine_gene_pval.Demo,scRNA.loc,X_super.Demo,file = "model_data.poet2.RData")
save.image(file = "model_data.Demo.RData")
save.image(file = "real_data.Demo.RData")
combineObject.Demo <- get_gene_chisq(gene_pval = combine_gene_pval.Demo, 
                                     gene.loc = scRNA.loc, # For scRNA-seq data, specify gene.loc = scRNA.loc, which is returned by function process_scRNA() 
                                     type = "iPOET", inter_dir = inter_dir)
save.image(file = "model_data.Demo.RData")
combine_POET_gene_pval.Demo <- combineObject.Demo$gene_POET_pval
combine_gene_corr.mat.Demo <- combineObject.Demo$gene_corr.mat
#5. Prioritizing trait-relevant tissue(s) and cell type(s)
gtex.enrichment.Demo <- prioritize_relevance(gene_pval = prunein_POET_gene_pval.Demo, gene_corr.mat = prunein_gene_corr.mat.Demo, 
                                             gene_expr = scRNA.rpkm, gene.loc = scRNA.loc, chrs =1:22) 

gtex.enrichment.joint.Demo <- prioritize_relevance(gene_pval = combine_POET_gene_pval.Demo, gene_corr.mat = combine_gene_corr.mat.Demo, 
                                                   gene_expr = scRNA.rpkm, gene.loc = scRNA.loc, chrs = 1:22)
save(gtex.enrichment.joint.Demo,file = "epic_modeldata_result.RData")
save.image(file = "model_data.Demo.RData")
save.image(file = "real_data.Demo.RData")
enrichment_plot.Demo <- plot_relevance(gtex.enrichment.joint.Demo)

```


```shell
#!/bin/bash
#SBATCH -e epic2.err
#SBATCH -o epic2.out
#SBATCH -J epic2
#SBATCH -w in006
#SBATCH --mem=30000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate R4.2
Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/2.r
Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/1.r
```


model data细胞类型pvalue可视化
```R
#epic细胞类型结果
load("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/epic_modeldata_result.RData")
#scPagwas结果
load("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/CelltypeP/Pagwas_monocytecount_modeldata_celltype.RData")
#sclinker结果
sclinker_p<-read.table("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/modeldata/results/modeldata_cell_types_pvalues.txt",header = T)
#sclinker_p$Cell_Type中的内容为prune_monocytes等，需要切割留下_后面的细胞类型
sclinker_p$Cell_Type<-gsub(".*_","",sclinker_p$Cell_Type)
sclinker<-sclinker_p[,c("Cell_Type","P")]
rownames(sclinker)<-sclinker_p$Cell_Type
scPagwas<-Pagwas$bootstrap_results[c(2:6),c("annotation","bp_value")]
rownames(scPagwas)<-scPagwas$annotation

df<-data.frame(scPagwas=-log2(scPagwas[names(gtex.enrichment.joint.Demo),2]),sclinker=-log2(sclinker[names(gtex.enrichment.joint.Demo),2]),epic=-log2(gtex.enrichment.joint.Demo),celltypes=names(gtex.enrichment.joint.Demo))
library(reshape2)
melted_df <- melt(df, id.vars = "celltypes")
library(ggplot2)
#绘制柱状图
pdf("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/EPIC/modeldata_celltype_pvalue.pdf")
ggplot(melted_df, aes(x = celltypes, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values = c("scPagwas" = "red", "sclinker" = "blue","epic"="green"))+
  labs(x = "Cell types", y = "-log2(p-value)")+
  theme(legend.position = "none")
  dev.off()