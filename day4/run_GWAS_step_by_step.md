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


# Run Association Test
