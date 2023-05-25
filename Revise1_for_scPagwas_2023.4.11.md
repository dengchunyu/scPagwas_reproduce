# scPagwas的review意见修稿

## 一、第一个revierw的回复工作

### 1.1 计算时间和内存的测试

#### 1.1.1 测试数据的产生

构造10w,20w,30w,40w,50w,60w,70w,80w,90w,100w的gwas数据，用于scPagwas的计算

```R
#读取monocytecount例子gwas数据，随机抽取100w的snp
#导入读取大文件的R包
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/TimeMemoryTest")
library(data.table)
gwas_file ="monocytecount_gwas_data.txt"
#读取gwas数据
gwas_data = fread(gwas_file)
#分别随机抽取100w,200w,300w,400w,500w,600w,700w,800w,900w,1000w的gwas数据，并输出txt文件
for(i in 1:10){
  gwas_data_sample = gwas_data[sample(nrow(gwas_data),i*1000000),]
  write.table(gwas_data_sample,file=paste0("gwas_TimeMemoryTest_",as.character(i*1000000),".txt"),sep="\t",quote=F,row.names=F,col.names=T)
}

```

构造1w，2w,3w,4w,5w,6w,7w,8w,9w,10w的单细胞数据，用于scPagwas的计算

```R
#library(scPagwas)
#查看scPgwas的版本
#packageVersion("scPagwas")
library(ggplot2)
suppressMessages(library(Seurat))
suppressMessages(library("dplyr"))
#读取bmmc的例子数据
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/TimeMemoryTest")
Single_data ="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/TimeMemoryTest/NM_Healthy_pbmc.rds"
#读取单细胞数据
Single_data = readRDS(Single_data)
#分别随机抽取1w,2w,3w,4w,5w,6w,7w,8w,9w,10w的单细胞数据，并输出rds文件
for(i in 1:10){
    print(i)
    if(i!=10){
        Single_data_sample = Single_data[,sample(1:ncol(Single_data),i*10000)]
        saveRDS(Single_data_sample,file=paste0("Single_TimeMemoryTest_",i*10000,".rds"))
    }else{
  saveRDS(Single_data,file=paste0("Single_TimeMemoryTest_",i*10000,".rds"))
        }
}

```

#### 1.1.2 利用上面得到的数据，循环计算scPagwas

```R
Args <- commandArgs(T)
gwasN = print(Args[1])
singleN = print(Args[2])
library(Seurat)
library(dplyr)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/TimeMemoryTest")
#例子
# Set the path to the code folder
code_folder <- "/share/pub/dengcy/GWAS_Multiomics/scPagwasCode"
# Get a list of all R files in the folder
r_files <- list.files(code_folder, pattern = "\\.R$", full.names = TRUE)
# Source each R file in the list
lapply(r_files, source)
load("block_annotation.RData")
load("Genes_by_pathway_kegg.RData")
load("chrom_ld.RData")

mem_before<-gc()
result <- system.time(
scPagwas_main(Pagwas =NULL,
                     gwas_data =file.path(glue::glue("gwas_TimeMemoryTest_{gwasN}.txt")),
                     Single_data =file.path(glue::glue("Single_TimeMemoryTest_{singleN}.rds")),
                     output.prefix=paste0("TimeMemoryTest_",gwasN,"_",singleN),
                    singlecell=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     output.dirs="TimeMemoryTest",
                     assay="RNA",
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
)
mem_after<-gc()
mem_used<-mem_after-mem_before
sink(file="TimeMemoryTest.txt",append=T)
cat(gwasN, "\t")
cat(singleN, "\t")
cat(result[["elapsed"]], "\t")
cat( mem_used[2,6], "\t")
sink()

```

服务器上运行
    
```bash
#!/bin/bash
#SBATCH -e job.err
#SBATCH -o job.out
#SBATCH -J job
#SBATCH --mem=150000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate R4.2
#1
for i in 1e+05 2e+05 3e+05 4e+05 5e+05 6e+05 7e+05 8e+05 9e+05 1e+06
do
    for j in 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000
    do
        Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/TimeMemoryTest/TimeMemoryTest.R $i $j
    done
done
#第二次计算
for i in 2e+06 3e+06 4e+06 5e+06
do
    for j in 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000
    do
        Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/TimeMemoryTest/TimeMemoryTest.R $i $j
    done
done
```

#### 1.1.3 计算时间和内存的结果可视化

时间可视化

```R
#分别读取gwas为10w，20w，30w，40w，50w，60w，70w，80w，90w，100w，单细胞为1w的时间和内存
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/TimeMemoryTest")
TimeMemo<-read.table("TimeMemoryTest.txt")
#将TimeMemo的第一行，变为4列的矩阵，按行读取
TimeMemo<-as.data.frame(matrix(TimeMemo,ncol=4,byrow=T))
colnames(TimeMemo)<-c("gwasN","singleN","Time","Memory")
#将TimeMemo的每一列变为数值型的向量
TimeMemo$gwasN<-as.numeric(TimeMemo$gwasN)
#将TimeMemo的gwasN的数值，转换为非科学计数法
TimeMemo$gwasN<-format(TimeMemo$gwasN,scientific=F)
TimeMemo$gwasN<-factor(TimeMemo$gwasN,levels=c(" 100000"," 200000"," 300000"," 400000"," 500000"," 600000"," 700000"," 800000"," 900000","1000000"))
TimeMemo$singleN<-as.numeric(TimeMemo$singleN)
TimeMemo$Time<-as.numeric(TimeMemo$Time)
TimeMemo$Memory<-as.numeric(TimeMemo$Memory)
#输出TimeMemo
write.csv(TimeMemo,file="TimeMemoryTest.csv")
write.csv(TimeMemo,file="TimeMemoryTest2.csv")
#####
TimeMemo1<-read.csv("TimeMemoryTest.csv")
TimeMemo2<-read.csv("TimeMemoryTest2.csv")
#合并
TimeMemo<-rbind(TimeMemo1,TimeMemo2)
#输出TimeMemo
write.csv(TimeMemo,file="TimeMemoryTest_cbind.csv")
#画折线图加点图，singleN为横坐标，Time为纵坐标，每个gwasN画一条折现图,点的颜色为gwasN
#导入ggplot2包
library(ggplot2)
#将gwasN的数值，转换为非科学计数法
TimeMemo$gwasN<-format(TimeMemo$gwasN,scientific=F)
TimeMemo$gwasN<-factor(TimeMemo$gwasN,levels=c(" 100000"," 200000"," 300000"," 400000"," 500000"," 600000"," 700000"," 800000"," 900000","1000000","2000000","3000000","4000000","5000000"))
#画图
ggplot(TimeMemo,aes(x=singleN,y=Time,color=gwasN))+geom_line()+geom_point(size=1)+theme_bw()+theme(legend.position="right")+labs(x="Number of cells",y="Time(s)",color="Number of SNPs")
#导出图片
ggsave("TimeMemoryTest_Time.pdf",width=8,height=6)
ggsave("TimeMemoryTest_Time.png",width=8,height=6)

```
内存可视化
  
```R
#画折线图，singleN为横坐标，Memory为纵坐标，每个gwasN画一条折现图
#导入ggplot2包
library(ggplot2)
ggplot(TimeMemo,aes(x=singleN,y=Memory,color=gwasN))+geom_line()+geom_point(size=1)+theme_bw()+theme(legend.position="right")+labs(x="Number of cells",y="Memory(MB)",color="Number of SNPs")
#导出图片
ggsave("TimeMemoryTest_Memory.pdf",width=8,height=6)
ggsave("TimeMemoryTest_Memory.png",width=8,height=6)
#单独看100w的gwasN
TimeMemo_100w<-subset(TimeMemo,gwasN=="1000000")
ggplot(TimeMemo_100w,aes(x=singleN,y=Memory,color=gwasN))+geom_line()+geom_point(size=1)+theme_bw()+theme(legend.position="right")+labs(x="Number of cells",y="Memory(MB)",color="Number of SNPs")
#导出图片
ggsave("TimeMemoryTest_Memory_100w.pdf",width=8,height=6)
ggsave("TimeMemoryTest_Memory_100w.png",width=8,height=6)

```

### 1.2 重新计算图3的单个细胞的pvalue

以下是用来产生背景校正pvalue的函数
Compute p-value from empirical null
    For score T and a set of null score T_1,...T_N, the p-value is

        p= [1 + \Sigma_{i=1}^N 1_{ (T_i \geq T) }] / (1+N)

    If T, T_1, ..., T_N are i.i.d. variables following a null distritbuion,
    then p is super-uniform.

    The naive algorithm is N^2. Here we provide an O(N log N) algorithm to
    compute the p-value for each of the N elements in v_t

下面是重新计算pvalue的函数，我们将函数整合进scPagwas当中。

#### 1.2.1 计算单细胞pvalue

```R
#1.scPagwas_monocytecount_modeldata_v1.10
source("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/NMFtest/scPagwasCode/Get_CorrectBg_p.R")
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/CorrectBGp")
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")
v_raw_score<-Pagwas@meta.data$scPagwas.TRS.Score1
gene_list <-names(Pagwas@misc$gene_heritability_correlation[order(Pagwas@misc$gene_heritability_correlation, decreasing = T), ])[1:1000]
#计算背景校正pvalue
correct_pdf<-Get_CorrectBg_p(singledata=Pagwas,
v_raw_score=v_raw_score, 
n_iters=1000,
n_genes=1000,
gene_list=gene_list
)

#输出correct_pdf到csv文件
write.csv(correct_pdf,file="scPagwas_monocytecount_modeldata_v1.10.0_correct_pvalue.csv",row.names = F)
#2.Pagwas_groundtruth_bmmc_monocyte_v10
load("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_groundtruth_bmmc_monocyte_v10.RData")
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/CorrectBGp")
v_raw_score<-Pagwas_groundtruth@meta.data$scPagwas.TRS.Score1
#计算背景校正pvalue
gene_list <-names(Pagwas_groundtruth@misc$gene_heritability_correlation[order(Pagwas_groundtruth@misc$gene_heritability_correlation, decreasing = T), ])[1:1000]
#计算背景校正pvalue
correct_pdf<-Get_CorrectBg_p(singledata=Pagwas_groundtruth,
v_raw_score=v_raw_score,
n_iters=1000,
n_genes=1000,
gene_list=gene_list
)
#输出correct_pdf到csv文件
write.csv(correct_pdf,file="Pagwas_groundtruth_bmmc_monocyte_v10_correct_pvalue.csv",row.names = F)
```

#### 1.2.2 计算细胞类型pvalue

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/CorrectBGp")
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")
single_p<-read.csv("scPagwas_monocytecount_modeldata_v1.10.0_correct_pvalue.csv")
Merge_celltype_p<-function(single_p,celltype){
    celltype_p<-data.frame(celltype=celltype,pvalue=single_p)
    celltype_p<-aggregate(pvalue~celltype,celltype_p,p_merge)
    return(celltype_p)
}
p_merge<-function(pvalues){
zvalues <- -sqrt(2) * qnorm(pvalues/2)
ztotal <- mean(zvalues)
p_total <- pnorm(-abs(ztotal))
return(p_total)
}

Merge_celltype_p(single_p$adj_p,Pagwas$celltype)

```


### 1.3 图3绘制校正后pvalue的结果

#### 1.3.1 图3模拟数据和真实数据，单细胞pvalue的结果

导入图3a的模拟数据的结果
    
```R
library(Seurat)
library(scPagwas)
library(ggplot2)
library(ggpubr)
library(ggrepel)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/CorrectBGp/")
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")
#导入校正后的pvalue
correct_pdf<-read.csv("scPagwas_monocytecount_modeldata_v1.10.0_correct_pvalue.csv")
all_fortify_can <- fortify.Seurat.umap(Pagwas)
all_fortify_can$CorrectBG_adj_p <- correct_pdf$adj_p
all_fortify_can$CorrectBG_z <- correct_pdf$pooled_z
#第一版pvalue的结果
p_thre <- 0.05;size <- 0.1
plots_sigp1 <- ggplot() +
    geom_point(
    data = all_fortify_can[all_fortify_can$CellScaleqValue > p_thre, ],
    aes(x = UMAP_1, y = UMAP_2), size = size, alpha = 0.8,
    color = "#E4DCCF"
    ) +
    umap_theme() +
    # new_scale_color() +
    geom_point(
    data = all_fortify_can[all_fortify_can$CellScaleqValue <= p_thre, ],
    aes(x = UMAP_1, y = UMAP_2), color = "#EA5455", size = .1
    ) +
    umap_theme() +
    # new_scale_color() +
    ggtitle(paste0("ScaleRank q value <", p_thre, " significant cells"))
#第二版pvalue的结果
plots_sigp2 <- ggplot() +
    geom_point(
    data = all_fortify_can[all_fortify_can$CorrectBG_adj_p > p_thre, ],
    aes(x = UMAP_1, y = UMAP_2), size = size, alpha = 0.8,
    color = "#E4DCCF"
    ) +
    umap_theme() +
    # new_scale_color() +
    geom_point(
    data = all_fortify_can[all_fortify_can$CorrectBG_adj_p <= p_thre, ],
    aes(x = UMAP_1, y = UMAP_2), color = "#EA5455", size = .1
    ) +
    umap_theme() +
    # new_scale_color() +
    ggtitle(paste0("CorrectBG adjp value <", p_thre, " significant cells"))
#整合两个版本的图片到一起
plots_sigp <- cowplot::plot_grid(plots_sigp1, plots_sigp2, ncol = 2)
#保存图片
ggsave("scPagwas_monocytecount_modeldata_v1.10.0_correct_pvalue.pdf", plots_sigp, width = 8, height = 4)

#计算假阳性率
df<-all_fortify_can[,c("CorrectBG_adj_p","CellScaleqValue","celltype")]

df$CorrectBG_adj_p<-ifelse(df$CorrectBG_adj_p<=0.05,1,0)
df$CellScaleqValue<-ifelse(df$CellScaleqValue<=0.05,1,0)
df$celltype<-ifelse(df$celltype=="monocytes",1,0)
#CorrectBG_adj_p的假阳性率
sum(df$CorrectBG_adj_p==1&df$celltype==0)/sum(df$CorrectBG_adj_p==1)
[1] 0.03732304
#CellScaleqValue的假阳性率
sum(df$CellScaleqValue==1&df$celltype==0)/sum(df$CellScaleqValue==1)
[1] 0.07578558
```

导入图3e的真实数据的结果

```R
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/CorrectBGp")
load("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_groundtruth_bmmc_monocyte_v10.RData")
#导入校正后的pvalue
correct_pdf<-read.csv("Pagwas_groundtruth_bmmc_monocyte_v10_correct_pvalue.csv")
all_fortify_can <- fortify.Seurat.umap(Pagwas_groundtruth)
all_fortify_can$CorrectBG_adj_p <- correct_pdf$adj_p
all_fortify_can$CorrectBG_z <- correct_pdf$pooled_z
#第一版pvalue的结果
p_thre <- 0.05;size <- 0.05
plots_sigp1 <- ggplot() +
    geom_point(
    data = all_fortify_can[all_fortify_can$CellScaleqValue > p_thre, ],
    aes(x = UMAP_1, y = UMAP_2), size = size, alpha = 0.8,
    color = "#E4DCCF"
    ) +
    umap_theme() +
    # new_scale_color() +
    geom_point(
    data = all_fortify_can[all_fortify_can$CellScaleqValue <= p_thre, ],
    aes(x = UMAP_1, y = UMAP_2), color = "#EA5455", size = size
    ) +
    umap_theme() +
    # new_scale_color() +
    ggtitle(paste0("ScaleRank q value <", p_thre, " significant cells"))
#第二版pvalue的结果
plots_sigp2 <- ggplot() +
    geom_point(
    data = all_fortify_can[all_fortify_can$CorrectBG_adj_p > p_thre, ],
    aes(x = UMAP_1, y = UMAP_2), size = size, alpha = 0.8,
    color = "#E4DCCF"
    ) +
    umap_theme() +
    # new_scale_color() +
    geom_point(
    data = all_fortify_can[all_fortify_can$CorrectBG_adj_p <= p_thre, ],
    aes(x = UMAP_1, y = UMAP_2), color = "#EA5455", size =size
    ) +
    umap_theme() +
    # new_scale_color() +
    ggtitle(paste0("CorrectBG adjp value <", p_thre, " significant cells"))
#整合两个版本的图片到一起
plots_sigp <- cowplot::plot_grid(plots_sigp1, plots_sigp2, ncol = 2)
#保存图片
ggsave("Pagwas_groundtruth_bmmc_monocyte_v10_correct_pvalue.pdf", plots_sigp, width = 8, height = 4)

#计算假阳性率
df<-all_fortify_can[,c("CorrectBG_adj_p","CellScaleqValue","celltypes")]

df$CorrectBG_adj_p<-ifelse(df$CorrectBG_adj_p<=0.05,1,0)
df$CellScaleqValue<-ifelse(df$CellScaleqValue<=0.05,1,0)
df$celltypes<-ifelse(df$celltypes %in% c("12_CD14.Mono.2","11_CD14.Mono.1","13_CD16.Mono"),1,0)
#CorrectBG_adj_p的假阳性率
sum(df$CorrectBG_adj_p==1&df$celltypes==0)/sum(df$CorrectBG_adj_p==1)
#CellScaleqValue的假阳性率

sum(df$CellScaleqValue==1&df$celltypes==0)/sum(df$CellScaleqValue==1)
[1] 0.01647655
sum(df$CellScaleqValue==1&df$celltypes==0)/sum(df$CellScaleqValue==1)
[1] 0.05832229
```

#### 1.3.2 图3模拟数据和真实数据，细胞类型pvalue的结果

```R
library(scPagwas)
library(Seurat)
library(dplyr)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/CorrectBGp")
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")
 Pagwas_data<-scPagwas_main(Pagwas = Pagwas,
                     gwas_data =system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas"), 
                     output.prefix="test", 
                     output.dirs="scPagwastest_output",
                     block_annotation = block_annotation,
                     assay="RNA", 
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     singlecell=F, 
                     celltype=T
)
```

### 1.4 利用NMF进行通路活性的打分，测试NMF的效果和原来的方法的效果

#### 1.4.1 nmf的计算，在图3的数据上进行计算

```R
#library(NMF)
#library(scPagwas)
library(Seurat)
library(dplyr)
##导入修改后的scPagwas相关的代码函数
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/NMFtest")
#例子
# Set the path to the code folder
code_folder <- "/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/NMFtest/scPagwasCode"
# Get a list of all R files in the folder
r_files <- list.files(code_folder, pattern = "\\.R$", full.names = TRUE)
# Source each R file in the list
lapply(r_files, source)
load("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/NMFtest/Genes_by_pathway_kegg.RData")
load("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/NMFtest/block_annotation.RData")
load("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/NMFtest/chrom_ld.RData")
#NMF结果
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")
Pagwas_nmf<-scPagwas_main(Pagwas = NULL,
                          gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                          output.prefix = "nmf",
                          output.dirs = "monocytecount_modeldata_nmf",
                          block_annotation = block_annotation,
                          Single_data = Pagwas,
                          assay = "RNA",
                          chrom_ld = chrom_ld,
                          Pathway_list = Genes_by_pathway_kegg,
                          pa_method='NMF',
                          n_topgenes = 1000)
save(Pagwas_nmf,file="Pagwas_monocytecount_modeldata_nmf.RData")



###真实数据的nmf结果
load("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_groundtruth_bmmc_monocyte_v10.RData")
#运行时间计算，并输出时间结果
system.time(
Pagwas_nmf<-scPagwas_main(Pagwas = NULL,
                          gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                          output.prefix = "nmf",
                          output.dirs = "monocytecount_realdata_nmf",
                          block_annotation = block_annotation,
                          Single_data = Pagwas_groundtruth,
                          assay = "RNA",
                          Pathway_list = Genes_by_pathway_kegg,
                          pa_method='NMF',
                          Correct_BG_p=TRUE,
                          chrom_ld = chrom_ld,
                          iters_celltype = 200,
                          iters_singlecell = 100,
                          n_topgenes = 1000)
)
save(Pagwas_nmf,file="Pagwas_monocytecount_realdata_nmf.RData")

#user   system  elapsed 
#7053.546  781.146 2559.615
```

服务器上运行
    
```bash
#!/bin/bash
#SBATCH -e model.err
#SBATCH -o model.out
#SBATCH -J model
#SBATCH -w in004
#SBATCH --mem=100000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate R4.2
#1
Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/NMFtest/model_nmf.r
#!/bin/bash
#SBATCH -e real.err
#SBATCH -o real.out
#SBATCH -J real
#SBATCH -w in004
#SBATCH --mem=100000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate R4.2
Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/NMFtest/real_nmf.r

```

#### 1.4.2 模拟数据和真实数据中NMF结果和原来结果的效果比对

```R
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
##导入修改后的scPagwas相关的代码函数
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/NMFtest")
#模拟数据
load("Pagwas_monocytecount_modeldata_nmf.RData")
#可视化
all_fortify_can <- fortify.Seurat.umap(Pagwas_nmf)

 df<-Pagwas_nmf@meta.data
  #将df按照scPagwas.TRS.Score1的列从大到小排序，获得celltype列中前50%的数据中monocytes元素的占比
  df<-df[order(df$scPagwas.TRS.Score1,decreasing = T),]
  df<-df[1:round(nrow(df)/2),]
  #计算monocytes在celltypes列中的占比
  df<-table(df$celltype)
  Sensitivity<- df["monocytes"]/sum(df)
  ti<-paste0("Model data for NMF method","\n","Sensitivity: ",Sensitivity)

p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas.TRS.Score1), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D",
                             midpoint = median(all_fortify_can$scPagwas.TRS.Score1))+
        theme(aspect.ratio=1) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle(ti)

pdf(file="Pagwas_monocytecount_modeldata_nmf.pdf",width =6, height =6)
print(p1)
dev.off()
png(file="Pagwas_monocytecount_modeldata_nmf.png")
print(p1)
dev.off()

#真实数据
load("Pagwas_monocytecount_realdata_nmf.RData")
#可视化
all_fortify_can <- fortify.Seurat.umap(Pagwas_nmf)
 df<-Pagwas_nmf@meta.data
  #将df按照scPagwas.TRS.Score1的列从大到小排序，获得celltype列中前50%的数据中monocytes元素的占比
  df<-df[order(df$scPagwas.TRS.Score1,decreasing = T),]
  df<-df[1:round(nrow(df)/2),]
  #计算monocytes在celltypes列中的占比
  df<-table(df$celltype)
  Sensitivity<- sum(df[c('11_CD14.Mono.1','12_CD14.Mono.2','13_CD16.Mono')])/sum(df)
  ti<-paste0("Real data for NMF method","\n","Sensitivity: ",Sensitivity)

p1<- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =scPagwas.TRS.Score1), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D",
                             midpoint = median(all_fortify_can$scPagwas.TRS.Score1))+
        theme(aspect.ratio=1) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle(ti)
pdf(file="Pagwas_monocytecount_realdata_nmf.pdf",width =6, height =6)
print(p1)
dev.off()
png(file="Pagwas_monocytecount_realdata_nmf.png")
print(p1)
dev.off()
```


### 1.5 细胞类型pvalue的计算

要求模拟数据，真实数据，所有血细胞trait都要计算细胞类型的pvalue

#### 1.5.1 模拟数据和真实数据中细胞类型pvalue的计算

```R
library(scPagwas)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/CelltypeP")
#模拟数据
load("/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/Pagwas_monocytecount_modeldata_v1.10.0.RData")
#将celltype中的11_CD14.Mono.1,11_CD14.Mono.2 改为11_CD14.Mono

Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data =Pagwas,
                     output.prefix="celltype", 
                     output.dirs="model_monocytecount_celltype",
                     block_annotation = block_annotation,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     singlecell=F,
                     celltype=T
)
save(Pagwas,file="Pagwas_monocytecount_modeldata_celltype.RData")
Pagwas$bootstrap_results$bp_value<-p.adjust(Pagwas$bootstrap_results$bp_value,method = "fdr")
#绘图
pdf("model_monocytecount_celltype.pdf")
Bootstrap_estimate_Plot(bootstrap_results=Pagwas$bootstrap_results,
                        width = 9,
                        height = 6,
                        do_plot=T)
dev.off()
#输出Pagwas$bootstrap_results结果
df<-Pagwas$bootstrap_results[-1,]
write.table(df,file="model_monocytecount_celltypeP.txt",sep="\t",quote=F,row.names = F)

#真实数据
load("/share/pub/dengcy/GWAS_Multiomics/realgroundtruth/Pagwas_groundtruth_bmmc_monocyte_v10.RData")
#删除celltypes为11_CD14.Mono.1的单细胞
Pagwas_groundtruth<-Pagwas_groundtruth[,which(Pagwas_groundtruth$celltypes!="11_CD14.Mono.1")]
Pagwas_groundtruth$celltypes[which(Pagwas_groundtruth$celltypes=="12_CD14.Mono.2")]<-"CD14.Mono"
Idents(Pagwas_groundtruth)<-Pagwas_groundtruth$celltypes
Pagwas<-scPagwas_main(Pagwas = NULL,
                     gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt",
                     Single_data =Pagwas_groundtruth,
                     output.prefix="celltype", 
                     output.dirs="real_monocytecount_celltype",
                     block_annotation = block_annotation,
                     Pathway_list=Genes_by_pathway_kegg,
                     chrom_ld = chrom_ld,
                     singlecell=F,
                     celltype=T
)
save(Pagwas,file="Pagwas_monocytecount_realdata_celltype.RData")
Pagwas$bootstrap_results$bp_value<-p.adjust(Pagwas$bootstrap_results$bp_value,method = "fdr")
pdf("real_monocytecount_celltype.pdf")
Bootstrap_estimate_Plot(bootstrap_results=Pagwas$bootstrap_results,
                        width = 9,
                        height = 7,
                        do_plot=T)
dev.off()
#输出Pagwas$bootstrap_results结果
df<-Pagwas$bootstrap_results[-1,]
write.table(df,file="real_monocytecount_celltypeP.txt",sep="\t",quote=F,row.names = F)

```

#### 1.5.2 所有血细胞trait的细胞类型pvalue的计算

结果在眼视光服务器存储

```R
#导出结果中的细胞类型结果
library(gplots)
library(RColorBrewer)
library(scPagwas)
 library(ggplot2)
 suppressMessages(library(Seurat))
 suppressMessages(library("dplyr"))
library('ComplexHeatmap')
library(circlize)
traits<-c("basophilcount","eosinophilcount" ,"Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")

result_TRS<-lapply(traits,function(i){
   print(i)
load(paste0("/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/",i,"_Hema_bmmc_scPagwas_v1.9.1.RData"))
    a<-data.frame(celltypes=Pagwas$celltypes,trs=Pagwas$scPagwas.topgenes.Score1)
     a$celltypes<-factor(a$celltypes,levels=unique(a$celltypes))
    #a have tow column, celltypes and trs, celltypes are characters, compute the mean trs for each celltype in a dataframe
    df<-aggregate(trs ~ celltypes, a, mean)
    return(df$trs)
})
result_TRS<-as.data.frame(result_TRS)
colnames(result_TRS)<-traits
rownames(result_TRS)<- unique(Pagwas$celltypes)
#保存
save(result_TRS,file="/share/pub/dengcy/GWAS_Multiomics/compare/Hema_test/result_TRS.RData")
load("D:/OneDrive/GWAS_Multiomics/Compare/Hema_test2/result_TRS.RData")
load("D:/OneDrive/GWAS_Multiomics/Compare/scPagwas_bmmc_list.RData")
result_list<-result_list[rownames(result_TRS),]
result_list_adj<- apply(result_list,2, function(x) p.adjust(x,method = "BH"))
#输出结果，输出的是result_TRS和result_list_adj的合并结果，
#新建的数据list每个元素是一个数据框，数据框包含两列，第一列是result_TRS的值，第二列是result_list_adj的值
list<-lapply(1:10,function(i){
  a<-data.frame(result_TRS[,i],result_list_adj[,i])
  colnames(a)<-paste0(colnames(result_TRS)[i],c("_mean_TRS","_FDR"))
  return(a)
})
#合并list中的数据框
result<-do.call(cbind,list)
#输出csv
write.csv(result,file="D:/OneDrive/GWAS_Multiomics/Manuscripts/Revise_comments/CelltypeP/Figure4b_bloodtrait_celltype_p.csv")



#coul <- colorRampPalette(brewer.pal(8, "Oranges"))(25)
#设定不同的颜色起始
coul <- colorRampPalette(c("#F5E8C7", "#AC7088"))(6)
p2star <- function(p){
  symnum(p,cutpoints = c(0,0.001,0.01,0.05,1),
         symbols = c('***','**','*',''),na = NA)
}
hM <- apply(result_list_adj,2, function(x) as.character(p2star(x)))
rownames(hM)<-rownames(result_list)

#定义，result_TRS中相同位置的值>0.2时，hM相同位置的值为""
hM[which(result_TRS < 0.2)]<-""

result_TRS<-as.matrix(result_TRS)
#删除行名为"14_Unk","26_Unk"的行
result_TRS<-result_TRS[!(rownames(result_TRS) %in% c("14_Unk","26_Unk")),]
result_TRS<-result_TRS[c("10_cDC","09_pDC","16_Pre.B","18_Plasma","17_B","13_CD16.Mono","12_CD14.Mono.2","11_CD14.Mono.1","25_NK","24_CD8.CM","22_CD4.M","23_CD8.EM","19_CD8.N","20_CD4.N1","21_CD4.N2","06_CLP.1","15_CLP.2","05_CMP.LMPP","07_GMP","01_HSC","04_Early.Baso",
"02_Early.Eryth","03_Late.Eryth","08_GMP.Neut"),c("Lymphocytecount3","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume","WhiteBloodCellcount","neutrophilcount","eosinophilcount","basophilcount","monocytecount")]
hM<-hM[c("10_cDC","09_pDC","16_Pre.B","18_Plasma","17_B","13_CD16.Mono","12_CD14.Mono.2","11_CD14.Mono.1","25_NK","24_CD8.CM","22_CD4.M","23_CD8.EM","19_CD8.N","20_CD4.N1","21_CD4.N2","06_CLP.1","15_CLP.2","05_CMP.LMPP","07_GMP","01_HSC","04_Early.Baso",
"02_Early.Eryth","03_Late.Eryth","08_GMP.Neut"),c("Lymphocytecount3","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume","WhiteBloodCellcount","neutrophilcount","eosinophilcount","basophilcount","monocytecount")]
colnames(result_TRS)<-c("Lymphocytecount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume","WhiteBloodCellcount","Neutrophilcount","Eosinophilcount","Basophilcount","Monocytecount")
pdf("D:/OneDrive/GWAS_Multiomics/Manuscripts/Revise_comments/CelltypeP/Figure_scPagwas_bmmc_cellytpe_p.pdf",
    height =8,width =6)
heatmap.2(result_TRS,
          trace="none",#
          col=coul,#
          density.info = "none",
          #设置颜色映射最小值
            breaks = seq(0, 0.6, length.out = 7),
          key.xlab ='Correlation',
          key.title = "",
          cexRow = 1,cexCol = 1,#
          Rowv = F,Colv =F, #
          margins = c(10,10),
          cellnote = hM,notecol='black'
)
dev.off()

```


## 子函数

```R
umap_theme <- function() {
  theme_grey() %+replace%
    theme(
      panel.background = element_rect(fill = "white", colour = "#4D4D4D", size =0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_blank()
    )
}

fortify.Seurat.umap <- function(x) {
  xy1 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "umap"))
  colnames(xy1) <- c("UMAP_1", "UMAP_2")
  xy1$UMAP_1 <- as.numeric(xy1$UMAP_1)
  xy1$UMAP_2 <- as.numeric(xy1$UMAP_2)

  return(cbind(xy1, as.data.frame(x@meta.data)))
}
```

