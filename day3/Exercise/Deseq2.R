if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("apeglm")

library("DESeq2")
library("ggplot2")
library("apeglm")

setwd("/Users/revathy/Personal/CAGT/GenomicResourceModules/Spring2025")

countData<- read.csv('data/GSE227234_RawCount_subset.csv', header = TRUE, sep = ",")
head(countData)

metaData <- read.csv('data/GSE227234_RawCount_metadata.csv', header = TRUE, sep = ",")
metaData

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~type, tidy = TRUE)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$type <- relevel(dds$type, ref="control")

dds <- DESeq(dds)

res <- results(dds)
head(results(dds, tidy=TRUE))

summary(res)

res <- results(dds, name="type_treated_vs_control")
res <- results(dds, contrast=c("type","control","treated"))

resultsNames(dds)

resLFC <- lfcShrink(dds, coef="type_treated_vs_control", type="apeglm")
resLFC

resOrdered <- res[order(res$padj),]
head(resOrdered)

#reset par
par(mfrow=c(1,1))

## Volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3), ylim=c(0,20)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

## PCA plots
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="type")

## MA plot
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

## adjust p-value to 0.05
res05 <- results(dds, alpha=0.05)
summary(res05)

## Volcano plot for p=0.05
with(res05, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3), ylim=c(0,20)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res05, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res05, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

## MA plot for p=0.05
plotMA(res05, ylim=c(-2,2))

