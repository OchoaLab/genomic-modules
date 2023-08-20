# check before running

# QC Steps
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
--make-pgen --out call_populations.v_s_cnt_hwe_filter.LD_pruned
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
sample_population_info_1kgp <- read_delim("sample_population_info_1kgp.tsv", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)
head(sample_population_info_1kgp)
```
load PCA results
```{r}
all_populations_v_s_cnt_hwe_filter_LD_pruned_Rel_pruned_pca <- read_delim("all_populations.v_s_cnt_hwe_filter.LD_pruned.Rel_pruned.pca.eigenvec",
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)
head(all_populations_v_s_cnt_hwe_filter_LD_pruned_Rel_pruned_pca)
```

use left_join() function to map the sample ID with population info
```{r}
pca_input <- left_join(all_populations_v_s_cnt_hwe_filter_LD_pruned_Rel_pruned_pca,
  select(sample_population_info_1kgp, Sample, Population, Super_Population),
  by = c("#IID" = "Sample")) %>% na.omit()
head(pca_input)
```

plot the first two PCs
```{r}
# label by population info
ggplot(pca_input,mapping=aes(x=PC1,y=PC2,color = Population)) + geom_point() + xlab("PC1") + ylab("PC2")
# label by super population info
ggplot(pca_input,mapping=aes(x=PC1,y=PC2,color = Super_Population)) + geom_point() + xlab("PC1") + ylab("PC2")
```
You should be able to see two plots, one with population info, one with super population info
<img width="818" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/f357a91d-a2f1-4115-95e6-5d37e6411631">

<img width="787" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/830813c2-4733-42e6-9ac0-3551f8942fcc">


# Run Association Test
