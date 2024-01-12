# Setup on your own machine

## install R and R Studio
Use the following [link](https://rstudio-education.github.io/hopr/starting.html). Note, you need to install **both** R and RStudio.

## Install packages
Within R, install all required packages.
```R
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("enrichR")
install.packages("devtools")
install.packages("knitr")
install.packages("patchwork")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("ComplexHeatmap")
BiocManager::install("gemma.R")
```
When you install packages with `BiocManager::install`, you may be asked whether you want to `Update packages all/some/none [a/s/n]`. You can put `n` in the console, to not update packages.

NOTE: You only have to install packages ONCE. From then on, you only have to load them using the functions `package()` or `require()`.

## Download data
Download the data using the following code:
```R
download.file("https://github.com/csbg/Hands-on-Biomedical-Data/raw/main/data/data.RDS", "data.RDS")
download.file("https://github.com/csbg/Hands-on-Biomedical-Data/raw/main/data/design.RDS","design.RDS")
download.file("https://github.com/csbg/Hands-on-Biomedical-Data/raw/main/data/gmap.RDS", "gmap.RDS")
```

## IMPORTANT NOTES!!!
- Define the paths for data (like `readRDS("data.RDS")`) based on the path on your computer (for example: readRDS("/home/nfortelny/Desktop/HandsOn/data.RDS")).
- Sometimes the `download.file` does not work properly. If the files can be found but not properly read by `readRDS()`, try to download them manually from this link: https://github.com/csbg/Hands-on-Biomedical-Data/tree/main/data. (Manually means to click on them and then click on "download".)
- Specific commands may differ, based on which packages you install. This will most likely only affect very recent packages (like `gemma.R`) - if at all. See instructions on day 1.

## Start exercises
You are now ready to start the exercises.