enrichment analysis for top 1000 gene for each method and each trait 
```R
Args <- commandArgs(T)
gwas = print(Args[1])
method = print(Args[2])
library(WebGestaltR)
os.file= file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/enrichment_result/"))
setwd(os.file)
enrich.file= file.path(glue::glue("/share/pub/dengcy/GWAS_Multiomics/MultiMethods/{method}/{gwas}.top1000.genes.csv"))
gene_df<-read.csv(enrich.file)
refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
outputDirectory <- getwd()
WebGestaltR(enrichMethod="ORA", organism="hsapiens",
enrichDatabase="geneontology_Biological_Process_noRedundant", interestGene=gene_df$gene_name,
interestGeneType="genesymbol", referenceGeneFile=refFile,
referenceGeneType="genesymbol", isOutput=TRUE,
outputDirectory=outputDirectory, projectName=gwas)
``````

## run cell markers for positive celltypes

```R
library(Seurat)
library(dplyr)

#read BMMC dataset
Single_data<-readRDS("/share/pub/dengcy/GWAS_Multiomics/singlecelldata/Seu_Hema_data.rds")

#merge the cluster of monocyte/lymphocyte/erythroid

levels(fc1)[grep('Mono',levels(fc1))] <- 'monocyte'
levels(fc1)[grep('14',levels(fc1))] <- 'monocyte'

levels(fc1)[grep('B',levels(fc1))] <- 'lymphocyte'
levels(fc1)[grep('18',levels(fc1))] <- 'lymphocyte'
levels(fc1)[grep('NK',levels(fc1))] <- 'lymphocyte'
levels(fc1)[grep('CD',levels(fc1))] <- 'lymphocyte'

levels(fc1)[grep('Eryth',levels(fc1))] <- 'erythroid'

Single_data@meta.data$'full_clustering2' <- fc1

#Findmarkers of monocyte,lymphocyte,erythroid
mono_marker <- FindMarkers(object = Single_data, ident.1 = 'monocyte',group.by = 'full_clustering2',logfc.threshold = 0) %>% filter(avg_log2FC > 0) 
lympho_marker <- FindMarkers(object = Single_data, ident.1 = 'lymphocyte',group.by = 'full_clustering2',logfc.threshold = 0) %>% filter(avg_log2FC > 0) 
erythroid_marker <- FindMarkers(object = Single_data, ident.1 = 'erythroid',group.by = 'full_clustering2',logfc.threshold = 0) %>% filter(avg_log2FC > 0) 

#save markers
write.csv(mono_marker,'./mono_marker.csv')
write.csv(lympho_marker,'./lympho_marker.csv')
write.csv(erythroid_marker,'./erythroid_marker.csv')

```

run the gene enrich analysis for trait-marker overlap genes

```R
library(dplyr)
library(readr)
library(WebGestaltR)

#read markers of monocyte, erythroid and lymphocyte
mono_marker <- read.csv('./mono_marker.csv')
lympho_marker <- read.csv('./lympho_marker.csv')
erythroid_marker <- read.csv('./erythroid_marker.csv')

#Run WebGestaltR for top1000 markers
WebGestaltR(enrichMethod="ORA", organism="hsapiens",
enrichDatabase="geneontology_Biological_Process_noRedundant", interestGene=mono_marker$X[1:1000],
interestGeneType="genesymbol", referenceGeneFile=refFile,
referenceGeneType="genesymbol", isOutput=TRUE,
outputDirectory='/share2/pub/zhouyj/zhouyj/scPAGWAS/blood/analysis', projectName=glue::glue("monocyte1000_enrich"))

WebGestaltR(enrichMethod="ORA", organism="hsapiens",
enrichDatabase="geneontology_Biological_Process_noRedundant", interestGene=erythroid_marker$X[1:1000],
interestGeneType="genesymbol", referenceGeneFile=refFile,
referenceGeneType="genesymbol", isOutput=TRUE,
outputDirectory='/share2/pub/zhouyj/zhouyj/scPAGWAS/blood/analysis', projectName=glue::glue("erythroid1000_enrich"))

WebGestaltR(enrichMethod="ORA", organism="hsapiens",
enrichDatabase="geneontology_Biological_Process_noRedundant", interestGene=lympho_marker$X[1:1000],
interestGeneType="genesymbol", referenceGeneFile=refFile,
referenceGeneType="genesymbol", isOutput=TRUE,
outputDirectory='/share2/pub/zhouyj/zhouyj/scPAGWAS/blood/analysis', projectName=glue::glue("lymphocyte1000_enrich"))


#read top 1000 significant pathway of lymphocyte
lympho_1000_enrich <- readr::read_tsv('./Project_lymphocyte1000_enrich/enrichment_results_lymphocyte1000_enrich.txt')
lympho_1000_enrich_0.05 <- lympho_1000_enrich %>% filter(FDR < 0.05)

#Find overlapped pathways of phenotypes related to lymphocyte
lympho_list <- c('Lymphocytecount3','LymphocytePercent')
lympho_result <- as.data.frame(matrix(0,nrow=2,ncol=12))
colnames(lympho_result) <- c('pheno','marker','magma','magma_overlap','scpagwas','scpagwas_overlap','smultixcan','smultixcan_overlap','spredixcan','spredixcan_overlap','TWAS','TWAS_overlap')
lympho_result$marker <- length(lympho_1000_enrich_0.05$geneSet)
for (i in c(1:2)){
    memo <- lympho_list[i]
    lympho_result$pheno[i] <- memo
    if(file.exists(glue::glue('./{memo}/Project_{memo}_magma1000/enrichment_results_{memo}_magma1000.txt'))){
        memo_magma_enrich <-readr::read_tsv(glue::glue('./{memo}/Project_{memo}_magma1000/enrichment_results_{memo}_magma1000.txt')) %>% filter(FDR <0.05)
        lympho_result$magma[i] <- length(memo_magma_enrich$geneSet)
        lympho_result$magma_overlap[i] <- length(intersect(lympho_1000_enrich_0.05$geneSet,memo_magma_enrich$geneSet))    
    }
    
    if(file.exists(glue::glue('./{memo}/Project_{memo}_scpagwas1000/enrichment_results_{memo}_scpagwas1000.txt'))){
        memo_scpagwas_enrich <- readr::read_tsv(glue::glue('./{memo}/Project_{memo}_scpagwas1000/enrichment_results_{memo}_scpagwas1000.txt')) %>% filter(FDR <0.05)
        lympho_result$scpagwas[i] <- length(memo_scpagwas_enrich$geneSet)
        lympho_result$scpagwas_overlap[i] <- length(intersect(lympho_1000_enrich_0.05$geneSet,memo_scpagwas_enrich$geneSet))
    }
    
    if(file.exists(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/smultixcan/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt'))){
        memo_smu_enrich <- readr::read_tsv(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/smultixcan/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt')) %>% filter(FDR <0.05)
        lympho_result$smultixcan[i] <- length(memo_smu_enrich$geneSet)
        lympho_result$smultixcan_overlap[i] <- length(intersect(lympho_1000_enrich_0.05$geneSet,memo_smu_enrich$geneSet))     
    }
    
    if(file.exists(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/spredixcan/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt'))){
        memo_spre_enrich <- readr::read_tsv(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/spredixcan/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt')) %>% filter(FDR <0.05)
        lympho_result$spredixcan[i] <- length(memo_spre_enrich$geneSet)
        lympho_result$spredixcan_overlap[i] <- length(intersect(lympho_1000_enrich_0.05$geneSet,memo_spre_enrich$geneSet)) 
     }
    
    if(file.exists(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/TWAS/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt'))){
        memo_twas_enrich <- readr::read_tsv(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/TWAS/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt')) %>% filter(FDR < 0.05)
        lympho_result$TWAS[i] <- length(memo_spre_enrich$geneSet)
        lympho_result$TWAS_overlap[i] <- length(intersect(lympho_1000_enrich_0.05$geneSet,memo_spre_enrich$geneSet))
    }
                                        
}

#read top 1000 significant pathway of monocyte
mono_1000_enrich <- readr::read_tsv('./Project_monocyte1000_enrich/enrichment_results_monocyte1000_enrich.txt')
mono_1000_enrich_0.05 <- mono_1000_enrich %>% filter(FDR < 0.05)

#Find overlapped pathways of phenotypes related to monocyte
mono_list <- c('monocytecount','basophilcount','eosinophilcount','neutrophilcount','WhiteBloodCellcount')
mono_result <- as.data.frame(matrix(0,nrow=5,ncol=12))
colnames(mono_result) <- c('pheno','marker','magma','magma_overlap','scpagwas','scpagwas_overlap','smultixcan','smultixcan_overlap','spredixcan','spredixcan_overlap','TWAS','TWAS_overlap')
mono_result$marker <- length(mono_1000_enrich_0.05$geneSet)
for (i in c(1:5)){
    memo <- mono_list[i]
    mono_result$pheno[i] <- memo
    if(file.exists(glue::glue('./{memo}/Project_{memo}_magma1000/enrichment_results_{memo}_magma1000.txt'))){
        memo_magma_enrich <-readr::read_tsv(glue::glue('./{memo}/Project_{memo}_magma1000/enrichment_results_{memo}_magma1000.txt')) %>% filter(FDR <0.05)
        mono_result$magma[i] <- length(memo_magma_enrich$geneSet)
        mono_result$magma_overlap[i] <- length(intersect(mono_1000_enrich_0.05$geneSet,memo_magma_enrich$geneSet))    
    }
    
    if(file.exists(glue::glue('./{memo}/Project_{memo}_scpagwas1000/enrichment_results_{memo}_scpagwas1000.txt'))){
        memo_scpagwas_enrich <- readr::read_tsv(glue::glue('./{memo}/Project_{memo}_scpagwas1000/enrichment_results_{memo}_scpagwas1000.txt')) %>% filter(FDR <0.05)
        mono_result$scpagwas[i] <- length(memo_scpagwas_enrich$geneSet)
        mono_result$scpagwas_overlap[i] <- length(intersect(mono_1000_enrich_0.05$geneSet,memo_scpagwas_enrich$geneSet))
    }
    
    if(file.exists(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/smultixcan/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt'))){
        memo_smu_enrich <- readr::read_tsv(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/smultixcan/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt')) %>% filter(FDR <0.05)
        mono_result$smultixcan[i] <- length(memo_smu_enrich$geneSet)
        mono_result$smultixcan_overlap[i] <- length(intersect(mono_1000_enrich_0.05$geneSet,memo_smu_enrich$geneSet))     
    }
    
    if(file.exists(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/spredixcan/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt'))){
        memo_spre_enrich <- readr::read_tsv(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/spredixcan/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt')) %>% filter(FDR <0.05)
        mono_result$spredixcan[i] <- length(memo_spre_enrich$geneSet)
        mono_result$spredixcan_overlap[i] <- length(intersect(mono_1000_enrich_0.05$geneSet,memo_spre_enrich$geneSet)) 
     }
    
    if(file.exists(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/TWAS/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt'))){
        memo_twas_enrich <- readr::read_tsv(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/TWAS/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt')) %>% filter(FDR < 0.05)
        mono_result$TWAS[i] <- length(memo_spre_enrich$geneSet)
        mono_result$TWAS_overlap[i] <- length(intersect(mono_1000_enrich_0.05$geneSet,memo_spre_enrich$geneSet))
    }
                                        
}

#read top 1000 significant pathway of erythroid
erythroid_marker_enrich <- readr::read_tsv('./Project_erythroid1000_enrich/enrichment_results_erythroid1000_enrich.txt')
erythroid_marker_enrich_0.05 <- erythroid_marker_enrich %>% filter(FDR < 0.05)

#Find overlapped pathways of phenotypes related to erythroid
erythroid_list <- c('Hemoglobinconcen','MeanCorpuscularHemoglobin','MeanCorpusVolume')
erythroid_result <- as.data.frame(matrix(0,nrow=3,ncol=12))
colnames(erythroid_result) <- c('pheno','marker','magma','magma_overlap','scpagwas','scpagwas_overlap','smultixcan','smultixcan_overlap','spredixcan','spredixcan_overlap','TWAS','TWAS_overlap')
erythroid_result$marker <- length(erythroid_marker_enrich_0.05$geneSet)
for (i in c(1:3)){
    memo <- erythroid_list[i]
    erythroid_result$pheno[i] <- memo
    if(file.exists(glue::glue('./{memo}/Project_{memo}_magma1000/enrichment_results_{memo}_magma1000.txt'))){
        memo_magma_enrich <-readr::read_tsv(glue::glue('./{memo}/Project_{memo}_magma1000/enrichment_results_{memo}_magma1000.txt')) %>% filter(FDR <0.05)
        erythroid_result$magma[i] <- length(memo_magma_enrich$geneSet)
        erythroid_result$magma_overlap[i] <- length(intersect(erythroid_marker_enrich_0.05$geneSet,memo_magma_enrich$geneSet))    
    }
    
    if(file.exists(glue::glue('./{memo}/Project_{memo}_scpagwas1000/enrichment_results_{memo}_scpagwas1000.txt'))){
        memo_scpagwas_enrich <- readr::read_tsv(glue::glue('./{memo}/Project_{memo}_scpagwas1000/enrichment_results_{memo}_scpagwas1000.txt')) %>% filter(FDR <0.05)
        erythroid_result$scpagwas[i] <- length(memo_scpagwas_enrich$geneSet)
        erythroid_result$scpagwas_overlap[i] <- length(intersect(erythroid_marker_enrich_0.05$geneSet,memo_scpagwas_enrich$geneSet))
    }
    
    if(file.exists(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/smultixcan/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt'))){
        memo_smu_enrich <- readr::read_tsv(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/smultixcan/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt')) %>% filter(FDR <0.05)
        erythroid_result$smultixcan[i] <- length(memo_smu_enrich$geneSet)
        erythroid_result$smultixcan_overlap[i] <- length(intersect(erythroid_marker_enrich_0.05$geneSet,memo_smu_enrich$geneSet))     
    }
    
    if(file.exists(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/spredixcan/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt'))){
        memo_spre_enrich <- readr::read_tsv(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/spredixcan/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt')) %>% filter(FDR <0.05)
        erythroid_result$spredixcan[i] <- length(memo_spre_enrich$geneSet)
        erythroid_result$spredixcan_overlap[i] <- length(intersect(erythroid_marker_enrich_0.05$geneSet,memo_spre_enrich$geneSet)) 
     }
    
    if(file.exists(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/TWAS/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt'))){
        memo_twas_enrich <- readr::read_tsv(glue::glue('/share/pub/dengcy/GWAS_Multiomics/MultiMethods/TWAS/enrichment_result/Project_{memo}/enrichment_results_{memo}.txt')) %>% filter(FDR < 0.05)
        erythroid_result$TWAS[i] <- length(memo_spre_enrich$geneSet)
        erythroid_result$TWAS_overlap[i] <- length(intersect(erythroid_marker_enrich_0.05$geneSet,memo_spre_enrich$geneSet))
    }
                                        
}

#combind overlap results
final_overlap <- rbind(lympho_result,erythroid_result,mono_result)
pure_overlap <- subset(final_overlap,select = c('pheno','magma_overlap','scpagwas_overlap','smultixcan_overlap','spredixcan_overlap','TWAS_overlap'))
pure_overlap$pheno <- c('Lymphocyte count','Lymphocyte percent','Hemoglobinconcen','Mean CorpuscularHemoglobin','MeanCorpus Volume','Monocyte count','Basophil count','Eosinophil count','Neutrophil count','WhiteBloodCell count')
gather_overlap <- gather(pure_overlap,method,value,magma_overlap,scpagwas_overlap,smultixcan_overlap,spredixcan_overlap,TWAS_overlap)
write.csv(gather_overlap,"./gather_overlap.csv")

```