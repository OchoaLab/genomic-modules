# Plink Data
Download data folder under day5

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
# Using Plink2
We will run all the plink2 commands in shell (command line). If you are not familiar with shell or command line, that's totally fine! 

Check this introduction: https://datacarpentry.org/shell-genomics/01-introduction.html

Once you get familiar with command line, you can check the plink2 documentations for file formats and availble functions: https://www.cog-genomics.org/plink/2.0/

# R and Rstudio Installation
Video for installing R and Rstudio for windows: https://www.youtube.com/watch?v=TFGYlKvQEQ4

Video for installing R and Rstudio for macbook: https://www.youtube.com/watch?v=AEebOXiMyyI

Not familiar with R? Check this: https://datacarpentry.org/genomics-r-intro/00-introduction.html

## R
R can be downloaded from https://cran.r-project.org/
click on the the correct link to match your own computer system and download the corresponding R
![image](https://github.com/OchoaLab/genomic-modules/assets/53951161/4965473e-9a2d-4eb1-9fff-5c452488cccb)
## Rstudio
You can search for rstudio.com as shown in the video and you will be directed to the new website: https://posit.co/products/open-source/rstudio/
Rstudio can be downloaded and installed using a similar way described in the videos above. 
## R packages
```{r}
install.packages("tidyverse")
```
```{r}
install.packages("qqman")
```

