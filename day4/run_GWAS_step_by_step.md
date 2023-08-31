# check before running
your input plink files are in the same folder as your plink2 package
# QC steps
run the following plink2 codes to filter out (1) variants with high missing rate (missing rate > 2%), (2) samples with high missing rate (>5%), 
(3) super rare variants (minor allele counts <10), and (4) variants that are not satisfied HWE (Hardy-Weinberg Equilibrium).
```
./plink2 --pfile all_populations --geno 0.02  --mind 0.05  --mac 10  --hwe 1e-50 \
--make-pgen --out all_populations.v_s_cnt_hwe_filter
```
# Population PCA Analysis
## Plink2 code
LD pruning
```
./plink2 --pfile all_populations.v_s_cnt_hwe_filter --indep-pairwise 1500 150 0.2 \
--make-pgen --out all_populations.v_s_cnt_hwe_filter.LD
./plink2 --pfile all_populations.v_s_cnt_hwe_filter --extract all_populations.v_s_cnt_hwe_filter.LD.prune.in  \
--make-pgen --out all_populations.v_s_cnt_hwe_filter.LD_pruned
```
filter on relatedness
```
./plink2 --pfile all_populations.v_s_cnt_hwe_filter.LD_pruned --king-cutoff 0.125 --out \
all_populations.v_s_cnt_hwe_filter.LD_pruned.kinship
./plink2 --pfile all_populations.v_s_cnt_hwe_filter.LD_pruned --keep \
all_populations.v_s_cnt_hwe_filter.LD_pruned.kinship.king.cutoff.in.id \
--make-pgen --out all_populations.v_s_cnt_hwe_filter.LD_pruned.Rel_pruned
```
run PCA
```
./plink2 --pfile all_populations.v_s_cnt_hwe_filter.LD_pruned.Rel_pruned --make-rel --pca --out \
all_populations.v_s_cnt_hwe_filter.LD_pruned.Rel_pruned.pca
```
You will get an output with the PCs for each sample:

all_populations.v_s_cnt_hwe_filter.LD_pruned.Rel_pruned.pca.eigenvec

## R code
load libraries
```{r}
# if not installed run: 
# install.packages("tidyverse")
library(tidyverse)
```
load sample information (sample ID, population, and super population)
```{r}
sample_population_info_1kgp <- read_tsv("sample_population_info_1kgp.tsv")
head(sample_population_info_1kgp)
```
load PCA results
```{r}
all_populations_v_s_cnt_hwe_filter_LD_pruned_Rel_pruned_pca <- read_tsv("all_populations.v_s_cnt_hwe_filter.LD_pruned.Rel_pruned.pca.eigenvec")
head(all_populations_v_s_cnt_hwe_filter_LD_pruned_Rel_pruned_pca)
```

use left_join() function to map the sample ID with population info
```{r}
pca_input <- left_join(all_populations_v_s_cnt_hwe_filter_LD_pruned_Rel_pruned_pca,
  sample_population_info_1kgp, by = c("#IID" = "Sample")) %>% na.omit()
head(pca_input)
```

plot the first two PCs
```{r}
# label by population info
ggplot(pca_input,mapping=aes(x=PC1,y=PC2,color = Population)) + geom_point()
# label by super population info
ggplot(pca_input,mapping=aes(x=PC1,y=PC2,color = Super_Population)) + geom_point()
```
You should be able to see two plots, one with population info, one with super population info
<img width="818" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/f357a91d-a2f1-4115-95e6-5d37e6411631">

<img width="787" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/830813c2-4733-42e6-9ac0-3551f8942fcc">


# Run Association Test
## Plink2 code
glm test with logistic regression (default)
```
plink2 --pfile all_populations.v_s_cnt_hwe_filter \
       --glm genotypic firth-fallback \
       --covar all_populations.v_s_cnt_hwe_filter.LD_pruned.Rel_pruned.pca.eigenvec \
       --parameters 1-8
```
then you will get plink2.PHENO1.glm.logistic.hybrid as the output of the glm test
rename it by
```
mv plink2.PHENO1.glm.logistic.hybrid all_populations.glm.logistic.hybrid
```

## R code
### general steps of checking glm results
load library
```{r}
library(tidyverse)
library(qqman)
```

load glm result file
```{r}
all_populations_glm_hybrid <- read_tsv("all_populations.glm.logistic.hybrid")
head(all_populations_glm_hybrid)
```
filter file for manhattan plot and qq plot
```{r}
results_filtered <- filter(all_populations_glm_hybrid, TEST == "ADD") %>% na.omit
nrow(results_filtered)
colnames(results_filtered)[1:3] <- c("CHR","BP","SNP")
```
manhattan plot
```{r}
manhattan(results_filtered)
```
<img width="757" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/ab57745b-8418-4542-a668-5c3029d01d02">

qq plot
```{r}
qq(results_filtered$P, main = "QQ Plot")
```
<img width="750" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/976db84e-5200-4de3-89e3-87ffc49f1b40">

### check the known associated variants
load variant file
```{r}
variant_ID_APOL1_G1_and_G2 <- read_lines("variant_ID_APOL1_G1_and_G2.txt")
```
filter glm results to see the pvalues for the two variants and the covariates (PC1 as example)
```{r}
filter(all_populations_glm_hybrid, ID %in% variant_ID_APOL1_G1_and_G2 & (TEST == "ADD" | TEST == "PC1"))[c("ID","TEST","P")]
```
<img width="503" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/f199f8d6-e522-46cc-8bc8-4d1f127fa2e1">

PC1 is highly associated with the phenotype meaning that there might be population stratification 
