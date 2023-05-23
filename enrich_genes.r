Args <- commandArgs(T)

gwas = print(Args[1])
method = print(Args[2])

gwas = "Lymphocytecount3"
method = "smultixcan"
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

