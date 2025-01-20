#####################################
### Install the required packages ###
#####################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("apeglm")

#########################
### Load the packages ###
#########################

library("DESeq2")
library("ggplot2")
library("ggrepel") 
library("apeglm")

#############################
### Set working directory ###
#############################

#### This should be set to the directory where you have the data downloaded into

setwd("</path/to/data_directory>")

#################
### Read Data ###
#################

countData<- read.csv('data/GSE227234_RawCount_subset.csv', header = TRUE, sep = ",")
head(countData)

metaData <- read.csv('data/GSE227234_RawCount_metadata.csv', header = TRUE, sep = ",")
metaData

###############################
### DESeq object initiation ###
###############################

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~type, tidy = TRUE)

keep <- rowSums(counts(dds)) >= 10  ## Keep only the rows whose sum is greater than or equal to 10
dds <- dds[keep,]

dds$type <- relevel(dds$type, ref="control")

################################
### Running the DESeq object ###
################################

dds <- DESeq(dds)

###########################
### Results and Sumamry ###
###########################

res <- results(dds)
head(results(dds, tidy=TRUE))

summary(res)

#######################
### Contrast matrix ###
#######################

res <- results(dds, name="type_treated_vs_control")
res <- results(dds, contrast=c("type","control","treated"))
resultsNames(dds)

#####################
### LFC Shrinkage ###
#####################

resLFC <- lfcShrink(dds, coef="type_treated_vs_control", type="apeglm")
resLFC


###########################
### Extracting DE Genes ###
###########################

resOrdered <- res[order(res$padj),]   ##order the results in the ascending order
head(resOrdered)

summary(resOrdered)
summary(resOrdered, alpha=0.05)

### What is the number of significant Differentially Expressed Genes (DEGs) in each case?

### Identifying the significantly upregulated and downregulated genes

LFC_threshold <- 0
padj_threshold <- 0.05

res_df <- as.data.frame(res)
res_df$gene_symbol <- rownames(res_df)
res_df$diffexpressed <- "Unregulated" #labelling all genes as NO
res_df$diffexpressed[res_df$log2FoldChange>LFC_threshold & res_df$padj<padj_threshold ] <- "Upregulated"
res_df$diffexpressed[res_df$log2FoldChange<LFC_threshold & res_df$padj<padj_threshold ] <- "Downregulated"
res_df$delabel <- ifelse(res_df$gene_symbol %in% head(res_df[order(res_df$padj), "gene_symbol"], 30), res_df$gene_symbol, NA)


#############
### Plots ###
#############

### MA plot

plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

### Volcano plot 

ggplot(data = res_df , aes(x= log2FoldChange , y= -log10(pvalue), col = diffexpressed, label=delabel))+
  geom_point()+
  theme_minimal()+
  scale_color_manual(values = c('blue','gray','red'))+
  theme(text = element_text(size = 15))+
  geom_vline(xintercept=c(-0.1, 0.1), col="violet", linetype='dashed') +
  geom_hline(yintercept=-log10(0.05), col="violet", linetype='dashed') +
  coord_cartesian(ylim = c(0, 15), xlim = c(-5,5)) + # since some genes can have minuslog10padj of inf, we set these limits
  geom_text_repel(max.overlaps = Inf) # To show all labels 

## To visualize APOL1 in the volcano plot
ggplot(data = res_df , aes(x= log2FoldChange , y= -log10(pvalue), col = diffexpressed, label=delabel))+
  geom_point()+
  theme_minimal()+
  scale_color_manual(values = c('blue','gray','red'))+
  theme(text = element_text(size = 15))+
  geom_vline(xintercept=c(-0.1, 0.1), col="violet", linetype='dashed') +
  geom_hline(yintercept=-log10(0.05), col="violet", linetype='dashed') +
  coord_cartesian(ylim = c(0, 20), xlim = c(-20,5)) + 
  geom_text_repel(max.overlaps = Inf)


### PCA plots
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="type")


