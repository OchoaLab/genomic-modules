# check before running
your input plink files are in the same folder as your plink2 package
# QC steps
run the following plink2 codes to filter out (1) variants with high missing rate (missing rate > 2%), (2) samples with high missing rate (>5%), 
(3) super rare variants (minor allele counts <10), and (4) variants that are not satisfied HWE (Hardy-Weinberg Equilibrium).
```
./plink2 --pfile African --geno 0.02 --mind 0.05 --mac 10 --hwe 1e-50 \
	--make-pgen --out African.filtered
```
filter on relatedness
```
./plink2 --pfile African.filtered --king-cutoff 0.125 --make-pgen --out \
	African.filtered.pruned
```
# Population PCA Analysis
## Plink2 code
run PCA
```
./plink2 --pfile African.filtered.pruned --pca --out \
	African.filtered.pruned.pca
```
You will get an output with the PCs for each sample: `African.filtered.pruned.pca.eigenvec`

## R code
load libraries
```{r}
# if not installed run: 
# install.packages("tidyverse")
library(tidyverse)
```
load sample information (sample ID, population, and super population)
```{r}
info <- read_tsv("sample_population_info_1kgp.tsv")
info
```
load PCA results
```{r}
data <- read_tsv("African.filtered.pruned.pca.eigenvec")
data
```

use left_join() function to map the sample ID with population info
```{r}
data <- left_join( data, info, by = c("#IID" = "Sample") )
data
```

plot the first two PCs
```{r}
# label by population info
ggplot( data, aes( x = PC1, y = PC2, color = Population ) ) + geom_point()
```
You should be able to see a plot of individuals in PCA space, colored according to their population
<img width="818" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/f357a91d-a2f1-4115-95e6-5d37e6411631">



# Run Association Test
## Plink2 code
glm test with logistic regression (default)
```
./plink2 --pfile African.filtered.pruned \
    --glm \
    --covar African.filtered.pruned.pca.eigenvec \
	--parameters 1 \
	--out African
```
then you will get `African.PHENO1.glm.logistic.hybrid` as the output of the glm test

## R code
### general steps of checking glm results
load library
```{r}
library(tidyverse)
library(qqman)
```

load glm result file
```{r}
data <- read_tsv("African.PHENO1.glm.logistic.hybrid")
data
```
edit colum names for manhattan plot and qq plot
```{r}
colnames(data)[1:3] <- c("CHR","BP","SNP")
```
manhattan plot
```{r}
manhattan(data)
```
<img width="757" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/ab57745b-8418-4542-a668-5c3029d01d02">

qq plot
```{r}
qq(data$P, main = "QQ Plot")
```
<img width="750" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/976db84e-5200-4de3-89e3-87ffc49f1b40">

### check the known associated variants
load variant file
```{r}
causal_variants <- read_lines("variant_ID_APOL1_G1_and_G2.txt")
```
filter glm results to see the pvalues for the two variants
```{r}
filter(data, SNP %in% causal_variants )[c("SNP","P")]
```
<img width="503" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/f199f8d6-e522-46cc-8bc8-4d1f127fa2e1">
