# Plink Data
# Plink2 Downloading and Installation

See downloading instructions: https://www.cog-genomics.org/plink/2.0/

## Downloading Examples:

Mac OS

M1
```
wget https://s3.amazonaws.com/plink2-assets/alpha4/plink2_mac_arm64_20230813.zip
unzip plink2_mac_arm64_20230813.zip
./plink2
```

And you should see something like
```
PLINK v2.00a4.5 M1 (13 Aug 2023)               www.cog-genomics.org/plink/2.0/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3

  plink2 <input flag(s)...> [command flag(s)...] [other flag(s)...]
  plink2 --help [flag name(s)...]

Commands include --rm-dup list, --make-bpgen, --export, --freq, --geno-counts,
--sample-counts, --missing, --hardy, --het, --fst, --indep-pairwise, --ld,
--sample-diff, --make-king, --king-cutoff, --pmerge, --pgen-diff,
--write-samples, --write-snplist, --make-grm-list, --pca, --glm, --adjust-file,
--gwas-ssf, --score, --variant-score, --genotyping-rate, --pgen-info,
--validate, and --zst-decompress.

"plink2 --help | more" describes all functions.
```
Intel
```
wget https://s3.amazonaws.com/plink2-assets/alpha4/plink2_mac_20230813.zip
unzip plink2_mac_20230813.zip
./plink2
```

Windows 64-bit
```
wget https://s3.amazonaws.com/plink2-assets/alpha4/plink2_win64_20230813.zip
unzip plink2_win64_20230813.zip
./plink2
```

# R and Rstudio Installation
## R
## Rstudio
## R packages
```{r}
install.packages("tidyverse")
```
```{r}
install.packages("qqman")
```

