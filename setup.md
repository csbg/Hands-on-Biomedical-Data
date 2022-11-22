# Setup on your own machine

## install R and R Studio
Use the following (link)[https://rstudio-education.github.io/hopr/starting.html]. Note, you need to install **both** R and RStudio.

## Install packages
Within R, install all required packages.
```R
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("enrichR")
install.packages("devtools")
install.packages("knitr")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("ComplexHeatmap")
BiocManager::install("gemma.R")
```

## Download data
Download the data using the following code
```R
download.file("https://github.com/csbg/Hands-on-Biomedical-Data/raw/main/data/data.RDS", "data.RDS")
download.file("https://github.com/csbg/Hands-on-Biomedical-Data/raw/main/data/design.RDS","design.RDS")
download.file("https://github.com/csbg/Hands-on-Biomedical-Data/raw/main/data/gmap.RDS", "gmap.RDS")
```

## Differences to Ricarda-based instructions
- Do not define `.libPaths("/home/handson/R/x86_64-pc-linux-gnu-library/4.0/")`, which is the path on Ricarda where packages are installed.
- Define the paths for data (like `readRDS("/home/handson/data/data.RDS")`) based on the path on your machine. We recommend downloading and saving all filese in folder, then no paths need to be defined, and you can just run `readRDS("data.RDS")` directly.
- Specific commands may differ, based on which packages you install. This will likely only affect more recent packages (like `gemma.R`).
