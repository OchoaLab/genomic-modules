# plink2 terminal commands

## Before running

Check that your input genotype files are in the same folder as your `plink2` executable  and run the following commands from that same directory on your terminal.

<img width="323" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/957c4c44-a5ea-47c0-84b1-ad5646d1abe9">

And navigate to the data folder. Not sure how to do that? check this https://github.com/OchoaLab/genomic-modules/blob/main/day5/Before%20Class%3A%20Download_and_Install.md#downloading-examples

## Data description

The genotype data is the African ancestry subset of 1000 Genomes, for chromosome 22 only, and is saved in this set of three files:

- `African.pgen`: compressed matrix of genotypes, encoded as zeroes, ones, twos, or missing value, for every individual and variant pair.
- `African.psam`: plain text table describing individuals, including their phenotype (which we simulated for this exercise to exactly match high-risk genotypes according to the standard APOL1 kidney disease model).
- `African.pvar`: plain text table describing variants (their IDs, locations, reference and alternate alleles).

There're two additional files we will use for this exercise:

- `variant_ID_APOL1_G1_and_G2.txt`: list of true causal variants for our simulated APOL1 kidney disease risk phenotype
- `sample_population_info_1kgp.tsv`: table that relates samples to their population labels.

## Quality control (QC) steps

Run the following `plink2` command to filter out :

1. variants with high missing rate (missing rate > 2%), 
2. samples with high missing rate (>5%), and 
3. super rare variants (minor allele counts <10).

```bash
./plink2 --pfile African --geno 0.02 --mind 0.05 --mac 10 \
	--make-pgen --out African.filtered
```
The terminal will output a report that looks something like the following, and also write the filtered genotype files (`African.filtered.pgen` and corresponding `psam` and `pvar` files).  Below you can see that the input data has 677 individuals (samples) and 106,657 variants (SNPs).  Further, 112 of our individuals are "cases" (fall in high risk category), while the rest (565) are "controls" (standard risk for kidney disease).  You'll notice that no samples or variants were removed due to missingness, because the 1000 Genomes data has no missingness (this is unusual of typical datasets).  However, a large number of variants (74,736) was removed because the minor allele was too rare (less than 10 counts).
```
PLINK v2.00a5LM AVX2 Intel (16 May 2023)       www.cog-genomics.org/plink/2.0/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to African.filtered.log.
Options in effect:
  --geno 0.02
  --mac 10
  --make-pgen
  --mind 0.05
  --out African.filtered
  --pfile African

Start time: Thu Sep 21 09:57:55 2023
15863 MiB RAM detected, ~9727 available; reserving 7931 MiB for main workspace.
Using up to 12 threads (change this with --threads).
677 samples (0 females, 0 males, 677 ambiguous; 677 founders) loaded from
African.psam.
106657 variants loaded from African.pvar.
1 binary phenotype loaded (112 cases, 565 controls).
Calculating sample missingness rates... done.
0 samples removed due to missing genotype data (--mind).
677 samples (0 females, 0 males, 677 ambiguous; 677 founders) remaining after
main filters.
112 cases and 565 controls remaining after main filters.
Calculating allele frequencies... done.
--geno: 0 variants removed due to missing genotype data.
74736 variants removed due to allele frequency threshold(s)
(--maf/--max-maf/--mac/--max-mac).
31921 variants remaining after main filters.
Writing African.filtered.psam ... done.
Writing African.filtered.pvar ... done.
Writing African.filtered.pgen ... done.
End time: Thu Sep 21 09:57:55 2023
```

## Principal component analysis (PCA)

Next we calculate the top 10 principal components (PCs), needed both to visualize population structure and to control when performing association.
The below command generates the file `African.filtered.eigenvec`, which calculates the PCs for each sample (they are technically "eigenvectors"), and `African.filtered.eigenval` are their corresponding eigenvalues (which are unused in this example).
For brevity the terminal output is not shown.
```bash
./plink2 --pfile African.filtered --pca --out \
	African.filtered
```

## Run Association Test

Lastly, the following command performs the desired genetic association test using a generalized linear model (GLM), which for a binary trait corresponds to logistic regression by default.
The output association table `African.filtered.PHENO1.glm.logistic.hybrid` is a plain text table with loads of statistics for each variant that was tested.
The terminal output is again omitted here.
```bash
./plink2 --pfile African.filtered \
    --glm \
    --covar African.filtered.eigenvec \
	--parameters 1 \
	--out African.filtered
```







# Data visualization in google colab (with R)

Here is the link to google colab https://colab.research.google.com/. Go to this webpage and login or create a new account.

## Open the GWAS_R_Visual.ipynb in google colab
1. Open a new notebook by:

![image](https://github.com/OchoaLab/genomic-modules/assets/53951161/3b26c690-7320-4b03-9d86-ace053785a03)


2. In the window that opens, choose `Upload`, then you can either click on browse and find your GWAS_Visual_R.ipynb or drag this file directly into the window.

![image](https://github.com/OchoaLab/genomic-modules/assets/53951161/5c6db674-ef74-400d-864b-f62c279c8544)

## Change runtime type:

1. Click `Runtime` and click `Change runtime type`
   
![image](https://github.com/OchoaLab/genomic-modules/assets/53951161/6e2de2d9-b407-4d6d-88e7-e321c06604dd)

2. In the window that opens, choose `R` in `Runtime type` and click `Save`

![image](https://github.com/OchoaLab/genomic-modules/assets/53951161/1c7780eb-c53d-4ac4-8abe-7c939a25a596)


## Run R command for each chunk

We will focus on the google colab for all the remaining steps. You can come back and check the exercise answers provided below.









# Data visualization in R (Optional - Answers at the end)

## Set working directory
This will help your R to locate the files directly:


In `Session`, click `Set Working Directory` then choose `Choose Directory`:

<img width="850" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/ae50f232-0cde-4566-8665-9401fbf0e037">


Navigate to the `data` folder, which should be right outside of the `GWAS` folder, and click `Open`:

<img width="788" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/72795693-185d-4ec2-ad91-b83c05aa56f9">


You will see this command (the path will be different in your R):

```
setwd("/Users/yunchengduan/OneDrive - Duke University/Projects/GenomicModule/data")
```
`setwd` means set-working-directory and the path to the directory is in the `()`.



Check the current working directory:

![image](https://github.com/OchoaLab/genomic-modules/assets/53951161/1ee22a0c-b039-439f-b9b4-bf2501cf5e00)



## Create a R markdown file:
1. click `New File` and choose `R markdown`:

   ![image](https://github.com/OchoaLab/genomic-modules/assets/53951161/44d8016f-64d9-48ec-8b2c-6dda7e053a2a)


3. Name your R markdown file as `GWAS_Visual`, click `OK`:

   <img width="580" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/ed5e1a24-acc5-4676-b646-c11d11c03440">


Inside the R notebook, replace the existing commands with:
```
---
title: "GWAS_Visual"
output: html_document
date: "2023-10-09"
---
```

For the following steps, we need to insert a chunk and run chunk by chunk. To insert a chunk, click on the green botton with `+c` meaning insert a new code chunk, then choose R:

![image](https://github.com/OchoaLab/genomic-modules/assets/53951161/d964cc7e-0537-48d7-bad7-535df578010d)


run the following commands to load the two R packages we need:
```{r}
# if not installed run: 
# install.packages("tidyverse")
# install.packages("qqman")
library(tidyverse)
library(qqman)
```
Click the `run` button to run the commands in the chunk and the output will be shown right below the chunk:

![image](https://github.com/OchoaLab/genomic-modules/assets/53951161/0844ac0b-1e2a-413c-89dc-3a01f1c88db1)




## PCA plot

Load the sample information for 1000 Genomes (sample ID, population, and super population):
```{r}
info <- read_tsv("sample_population_info_1kgp.tsv")
info
```
Load PCA results:
```{r}
data <- read_tsv("African.filtered.eigenvec")
data
```
Use the `left_join()` function to map the sample ID to their corresponding population info:
```{r}
data <- left_join(data, info, by = c("#IID" = "Sample"))
# better way to inspect big tables:
View( data )
```
Plot the first two PCs, which reveal the principal axes of variation of allele frequencies in our African samples.  You should see the following plot of individuals in PCA space, colored according to their population.  PC1 distinguishes admixed populations (ASW, ACB, both of which have some European ancestry) from the rest of the populations, which are practically unadmixed.  In contrast, PC2 distinguishes far Western African ancestry (GWD and MSL) from central populations, which also share linguistic connections.
```{r}
# label by population info
ggplot(data, aes(x = PC1, y = PC2, color = Population)) + geom_point()
```
![image](https://github.com/OchoaLab/genomic-modules/assets/53951161/a445926a-ca23-421f-8d56-4d5a2f9489ed)


## Manhattan plot

Load the association table:
```{r}
result <- read_tsv("African.filtered.PHENO1.glm.logistic.hybrid")
result
```
Load the list of causal variants:
```{r}
causal_variants <- read_lines("variant_ID_APOL1_G1_and_G2.txt")
```
Edit the column names to match what the `manhattan` function prefers:
```{r}
colnames(result)[1:3] <- c("CHR","BP","SNP")
```
Make the Manhattan plot, which reveals the locations of associated variants.  Each point is a variant, points above the red line are usually significant, and points above the blue line are considered "suggestive".  Since we simulated the trait, we know which are the two true causal variants, which we highlighted in green.  Note that several non-causal variants are also significant, which are neighbors of the causal ones and appear significant due to linkage disequilibrium (LD).  The smaller peak on the far right is not significant and is not a true association either (it's important not to overinterpret such weaker peaks)!
```{r}
manhattan(result, highlight = causal_variants)
```
![image](https://github.com/OchoaLab/genomic-modules/assets/53951161/77b69b3e-347f-4e0c-8bde-d1756d57473f)


## Quantile-Quantile (QQ) plot

Lastly, we make the QQ plot, which tells if the association study as a whole appears accurate or not.
On the x-axis are the expected quantiles of p-values when there is no association (p-values should be uniform), and on the y-axis are the observed quantiles.
When the test is well calibrated, the data curve is near the y=x (red) line for large p-values (small -log10 p-values), and as p-values become smaller there is greater departure from the red line because of the true associations.
However, if the test is not well-calibrated (often called "inflated"), which happens for a variety of interesting reasons, the data curve departs from the red line much sooner than expected, in which case null p-values are not actually uniform and the significant results cannot be trusted (technically, the test is not correctly controlling the type I error).
This case looks quite good considering our small sample size and the skewed proportion of cases and controls (ideally it should be 1:1 for this test).
```{r}
qq(result$P, main = "QQ Plot")
```
![image](https://github.com/OchoaLab/genomic-modules/assets/53951161/0c4d12e6-6e74-4e51-99b4-f1350d6e3f72)


