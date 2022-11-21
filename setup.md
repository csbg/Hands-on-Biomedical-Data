# Setup on your own machine

## install R and R Studio
Use the following (link)[https://rstudio-education.github.io/hopr/starting.html]. Note, you need to install both R and RStudio.

## Install packages
```R
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("enrichR")
install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("ComplexHeatmap")
BiocManager::install("gemma.R")
```

## Download data
Data can be downloaded under this (link)[https://filesender.aco.net/?s=download&token=cb05cfa9-68c8-4fd1-b0df-5d475132768c]