# DESeq2 Tutorial for Differential Gene Expression Analysis

## Overview
This repository contains a tutorial for using the DESeq2 package in R to identify differentially expressed genes (DEGs) from count data. The tutorial includes:

- An R script (`Deseq2.R`) demonstrating DESeq2 workflow.
- Subsetted count data from the GEO dataset [GSE227234](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227234).

## About DESeq2
DESeq2 is a widely used R package for differential gene expression analysis of RNA-seq data. It employs statistical models to analyze count data, providing robust results even with small sample sizes. Key features include:

- Normalization of read counts.
- Shrinkage of effect size estimates for improved accuracy.
- Identification of statistically significant DEGs.
- Visualization tools for exploratory data analysis and result presentation.

For more information, visit the [DESeq2 documentation](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

## Glossary
- **DEG (Differentially Expressed Gene):** A gene that shows statistically significant differences in expression levels between experimental groups.
- **Log Fold Change (LFC):** The logarithmic ratio of gene expression between conditions.
- **P-value:** The probability of observing the results if the null hypothesis is true.
- **Adjusted P-value (padj):** A corrected p-value accounting for multiple testing to control the false discovery rate (FDR).
- **MA Plot:** A scatterplot of log ratios (M) versus average expression levels (A).
- **Volcano Plot:** A scatterplot showing the relationship between statistical significance (-log10(p-value)) and effect size (log2FoldChange).
- **PCA (Principal Component Analysis):** A dimensionality reduction method to visualize variation in gene expression across samples.

## Data
This tutorial uses a subset of count data from the GEO dataset [GSE227234](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227234). The data includes:

1. `GSE227234_RawCount_subset.csv`: A matrix of raw read counts for selected genes and samples.
2. `GSE227234_RawCount_metadata.csv`: Metadata describing sample groups and experimental conditions.

### Data Citation
> Yoshida, T., Latt, K. Z., Santo, B. A., Shrivastav, S., Zhao, Y., Fenaroli, P., ... & Kopp, J. B. (2023). APOL1 kidney risk variants in glomerular diseases modeled in transgenic mice. bioRxiv.
> GEO accession: GSE227234. Available at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227234

## Contents
- **`Deseq2.R`:** The main script demonstrating the DESeq2 workflow, including:
  - Data preprocessing and filtering.
  - Differential expression analysis.
  - Visualization of results (MA plot, volcano plot, PCA).

## Prerequisites
Ensure the following R packages are installed:
- `DESeq2`
- `ggplot2`
- `ggrepel`
- `apeglm`

Install these packages using BiocManager if needed:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("apeglm")
```

## Running the Tutorial
1. Clone the repository and set the working directory in `Deseq2.R` to the folder containing the data.
2. Run the R script step by step to:
   - Load and preprocess the data.
   - Perform differential expression analysis.
   - Generate visualizations for insights into the data and results.

## DESeq2 Citation
Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology, 15:550. [DOI:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

