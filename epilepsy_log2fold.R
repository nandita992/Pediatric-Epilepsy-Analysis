"Project : Analysis of mRNA as therapeutics drugs" 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("GEOquery")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("WGCNA")

install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival", "BiocManager"))
BiocManager::install(c("GO.db", "preprocessCore", "impute"));

devtools::install_github("kevinblighe/CorLevelPlot")

devtools::install_github('kevinblighe/EnhancedVolcano')

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
require(GEOquery)
library(dplyr)
library(tidyr)
library(gridExtra)
library(EnhancedVolcano)

allowWGCNAThreads()     
"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94744"

# 1. Data Preparation ------------------------------------------------

'Step 1: Import the Data'
data <- read.delim('/Users/nanditapuri/Downloads/GSE193842_All_miRNA.UMI.TPM.list.txt', 
                   header=T)

' Step 2: Import PhenoData '
geo_id <- "GSE193842"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))
phenoData <- phenoData[,c(1,2,42)]


#2. Data Manipulation -------------------------------------------------

'Step 1: Change Pheno Data column names and manipulate it'

rownames(phenoData)<-NULL


phenoData <- phenoData %>%
  mutate(title = str_sub(title, 6,-4)) %>%
  select(1,3) %>%
  column_to_rownames(var='title') %>%
  dplyr::rename("disease_stage" = "disease state:ch1")

'Step 1: Change Data column names and manipulate it'

data <- data %>%
  gather(key='samples',value = 'counts',-X.ID) %>%
  mutate(samples= gsub('\\.', '-',samples)) %>%
  spread(key='samples',value='counts') %>%
  dplyr::rename("mRNA" = "X.ID") %>%
  column_to_rownames(var= 'mRNA')



#3 Outlier Detection ----------------------------------------

"Step 1: Use Hiearchichal Clustering for Outlier Detection"

htree <- hclust(dist(t(data)), method = "average")
plot(htree)


"Step 2: Use PCA for Outlier Detection"

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


'Step 3: Exclude outlier samples'

samples.to.be.excluded <- c('A328-24','A326-11')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]

phenoData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)


# 3. Deseq Object----------------------------------------------------------------------


"Step 1 : Sanity Check "

# making the rownames and column names identical
all(rownames(phenoData) %in% colnames(data.subset))
all(rownames(phenoData) == colnames(data.subset))

"Step 2: Create Deseq Object"
dds <- DESeqDataSetFromMatrix(countData = round(data),
                              colData = phenoData,
                              design = ~ disease_stage)

dds$disease_stage <- relevel(dds$disease_stage, ref = "healthy")


dds <- DESeq(dds)
res <- results(dds)

res <- as.data.frame(res)

write.csv(res,"/Users/nanditapuri/Documents/Codes/Epilepsy/log2fold.csv", row.names = TRUE)


'Step 3: Shortlist Up-regulated Genes according to log2foldchange and p-values'
resSigind_up = res[ which(res$pvalue < 0.05 & res$log2FoldChange >= 2), ]
resSigind_up
write.csv(resSigind_up,"/Users/nanditapuri/Documents/Codes/Epilepsy/up_regulated_genes_p_value.csv", row.names = TRUE)


'Step 4: Shortlist Down-regulated Genes according to log2foldchange and p-values'
resSigind_down = res[ which(res$pvalue < 0.05 & res$log2FoldChange < -1), ]
resSigind_down 
write.csv(resSigind_down,"/Users/nanditapuri/Documents/Codes/Epilepsy/down_regulated_genes_p_value.csv", row.names = TRUE)


#5 . Visualization -------------------------

'Visualization 1: Volcano Plot'

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'N061011 versus N61311',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0)



  

