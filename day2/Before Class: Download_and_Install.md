To make the most of our lecture time, it'd be awesome if you could download and install everything using the provided instructions before we meet. This way, we can dive straight into the good stuff during the lecture.
All you'll need to do/try is to (1) download files with data, (2) download and install plink2 (3) download and install R, RStudio, and Rscript, and/or download our colab notebook.

# Plink Data
Download the data folder and Result folder under day2
## option 1
you can download the whole repo:
![image](https://github.com/user-attachments/assets/377ce83a-f464-44b9-98bb-ef54b2483069)


## option 2
Nevigate to genomic-modules/day5/data, and download the files in red frames one by one:
![image](https://github.com/user-attachments/assets/eb905486-26da-4d64-9914-0dc12f7e3e60)

Click on the file name and then click the download button:
![image](https://github.com/user-attachments/assets/72b07bb6-29a6-4371-a8ef-4a4d5aec46de)


You can use the same steps to download `Result` file under day5

# Plink2 Downloading and Installation
We will run all the plink2 commands in shell (command line). If you are not familiar with shell or command line, that's totally fine! 

Check this introduction: https://datacarpentry.org/shell-genomics/01-introduction.html (highly recommend)

See downloading instructions: https://www.cog-genomics.org/plink/2.0/

## Downloading Examples:

### Mac OS

1. Open a terminal on your macbook (you can find it in Applications -> Utilities -> Terminal)

<img src="https://github.com/OchoaLab/genomic-modules/assets/53951161/dc1beb98-7106-47a7-9914-4a58744a778a" width="800">

Alternatively, you can find it in the launchpad >> others 

2. Once you opened the terminal, nativegate to the directory where your `data` folder located using `cd` command. 

You can manually type the path to the directory, or you can drag the folder after the `cd` command. (Note that there should be a `space` between `cd` and the path.

<img src="https://github.com/OchoaLab/genomic-modules/assets/53951161/08086890-6f0b-4e87-961c-fd96c0fbe869" width="800">


and you will get this：

<img width="800" alt="image" src="https://github.com/OchoaLab/genomic-modules/assets/53951161/b4686559-d09e-41b2-ab9d-e38039ec375a">

3. Hit `return` and you should be in the directory which has the plink files. Note that the path to the directory on your terminal will look different to what I have in the picture above. 

4. Check where you're using `pwd` and the files under the directory using `ls`. 

  Again, it's totally fine if you're not familiar with those commands. Reading through this introduction https://datacarpentry.org/shell-genomics/01-introduction.html is highly recommended. It will help you get prepared for this session. 


5. Copy this in your terminal to download and install plink2

For M1-Chip
```
wget https://s3.amazonaws.com/plink2-assets/alpha4/plink2_mac_arm64_20230813.zip
unzip plink2_mac_arm64_20230813.zip
./plink2
```

For Intel-Chip
```
wget https://s3.amazonaws.com/plink2-assets/alpha4/plink2_mac_20230813.zip
unzip plink2_mac_20230813.zip
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


**Get error message? command not found: wget? Maybe you don't have wget installed in your laptop.**

1. Copy https://s3.amazonaws.com/plink2-assets/alpha4/plink2_mac_arm64_20230813.zip into your browser if your computer has a M1 chip (or https://s3.amazonaws.com/plink2-assets/alpha4/plink2_mac_20230813.zip if you're using intel-chip), hit `return` and in this way you can download plink2 as well.

2. Drag the zip file you just downloaded into the `data` folder which contains your plink files
   
3. Copy this in your terminal and install plink2
   
For M1-Chip 
```
unzip plink2_mac_arm64_20230813.zip
./plink2
```

For Intel-Chip
```
unzip plink2_mac_20230813.zip
./plink2
```



### windows
#### Windows 64-bit
```
wget https://s3.amazonaws.com/plink2-assets/alpha4/plink2_win64_20230813.zip
unzip plink2_win64_20230813.zip
./plink2
```

# Using Plink2

You can check the plink2 documentations for file formats and availble functions: https://www.cog-genomics.org/plink/2.0/


# R and Rstudio Installation 

Video for installing R and Rstudio for macbook: https://www.youtube.com/watch?v=AEebOXiMyyI

Video for installing R and Rstudio for windows: https://www.youtube.com/watch?v=TFGYlKvQEQ4

Not familiar with R? Check this: https://datacarpentry.org/genomics-r-intro/00-introduction.html

## R
R can be downloaded from https://cran.r-project.org/
click on the the correct link to match your own computer system and download the corresponding R
![image](https://github.com/OchoaLab/genomic-modules/assets/53951161/4965473e-9a2d-4eb1-9fff-5c452488cccb)
## Rstudio
You can search for rstudio.com as shown in the video and you will be directed to the new website: https://posit.co/products/open-source/rstudio/
Rstudio can be downloaded and installed using a similar way described in the videos above. 

## R packages

Before installing the packages, we will need to create a R project. If you are not familiar with what R project is and how to create a R project, check this： https://datacarpentry.org/genomics-r-intro/00-introduction.html#create-an-rstudio-project

### Create R project
1. Double click on the Rstudio to open it
   
2. To create a project, go to the File menu, and click `New Project…`. In the window that opens, click `New Directory`.
   
   <img src="https://github.com/OchoaLab/genomic-modules/assets/53951161/8e292186-9ade-4a4c-9556-951597af2aa9" width="500">

   Then click `New Project `

   <img src="https://github.com/OchoaLab/genomic-modules/assets/53951161/b2f0137e-c7a0-4e20-8ddd-befe5cc0be3d" width="500">


3. For “Directory name:” enter `GWAS`. For “Create project as subdirectory of”, click `Browse…` and then nevigate to where you downloaded the plink files and plink2 then clink `Choose`.

   <img src="https://github.com/OchoaLab/genomic-modules/assets/53951161/118a0585-69c0-4165-a5e1-ed3ba3b0ee16" width="500">


### Installation

You can type the following commands into the console directlty or in a Rscript or Rmarkdown. We will use Rmarkdown to run all the r commands in the exercise.
```{r}
install.packages("tidyverse")
```
```{r}
install.packages("qqman")
```
To check if the libraries are installed

```{r}
library("tidyverse")
library("qqman")
```

# Download R for co-lab (Optional)

If your installation of R/Rstudio is unsuccessful, you can still use google co-lab to run R commands. You can find the R code in `GWAS_Visual_R.ipynb` under day5 in the repo. Try to download it use the same way described above (if you chose option 1 in the last step, you already have it!)
