### 1.1 test for time and memory

#### 1.1.1 test data
```R
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/TimeMemoryTest")
library(data.table)
gwas_file ="monocytecount_gwas_data.txt"
gwas_data = fread(gwas_file)
for(i in 1:10){
  gwas_data_sample = gwas_data[sample(nrow(gwas_data),i*1000000),]
  write.table(gwas_data_sample,file=paste0("gwas_TimeMemoryTest_",as.character(i*1000000),".txt"),sep="\t",quote=F,row.names=F,col.names=T)
}

```

```R
#library(scPagwas)
#packageVersion("scPagwas")
library(ggplot2)
suppressMessages(library(Seurat))
suppressMessages(library("dplyr"))
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/TimeMemoryTest")
Single_data ="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/TimeMemoryTest/NM_Healthy_pbmc.rds"
Single_data = readRDS(Single_data)

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

#### 1.1.2 run scPagwas

```R
Args <- commandArgs(T)
gwasN = print(Args[1])
singleN = print(Args[2])
library(Seurat)
library(dplyr)
setwd("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/TimeMemoryTest")
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

linux server
    
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
for i in 1e+05 2e+05 3e+05 4e+05 5e+05 6e+05 7e+05 8e+05 9e+05 1e+06 2e+06 3e+06 4e+06 5e+06
do
    for j in 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000
    do
        Rscript /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/TimeMemoryTest/TimeMemoryTest.R $i $j
    done
done
```