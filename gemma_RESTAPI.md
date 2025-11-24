# Downloading data from GEMMA using the REST-API

GEMMA provides approaches to download data and metadata for any dataset from the webbrowser.
This is called [representational state transfer (REST) application programming interface (API)](https://gemma.msl.ubc.ca/resources/restapidocs/#/).
As part of this, we will create URLs (Uniform Resource Locator). For example, for the dataset `GSE2018`:
* This URL downloads the design (metadata): `https://gemma.msl.ubc.ca/rest/v2/datasets/GSE2018/design`
* This URL downloads the expression data: `https://gemma.msl.ubc.ca/rest/v2/datasets/GSE2018/data/processed`

## Setup

Install packages
```{r}
install.packages("httr")
install.packages("jsonlite")
```

Load packages
```{r}
library(httr)
library(jsonlite)
library(tidyverse)
```

Define urls
```{r} 
gemma.url <- "https://gemma.msl.ubc.ca/rest/v2/datasets"
gse <- "GSE2018"
```

## Download

If this doesn't work it can be skipped by using the web browser instead.

Download metadata of the URL `https://gemma.msl.ubc.ca/rest/v2/datasets/GSE2018/design`.
```{r}
url.metadata <- paste(gemma.url, gse, "design", sep="/")
response <- GET(url.metadata)
writeBin(content(response, "raw"), paste0(gse, "_metadata.tsv"))
```

Download data of the URL `https://gemma.msl.ubc.ca/rest/v2/datasets/GSE2018/data/processed`
```{r}
url.data <- paste(gemma.url, gse, "data/processed", sep="/")
response <- GET(url.data)
writeBin(content(response, "raw"), paste0(gse, "_data.tsv"))
```

## Parse data

First the metadata ("obs" in AnnData).
```{r}
metadata <- read_tsv(paste0(gse, "_metadata.tsv"), comment = "#")
head(metadata)
metadata <- metadata |> 
  column_to_rownames("ExternalID")
head(metadata)
```

Next the data itself.
```{r}
data <- read_tsv(paste0(gse, "_data.tsv"), comment = "#")
str(data)
```

The data contains data and gene information. We will first format a gene map ("var" in AnnData).
```{r}
gmap <- data |> 
  select(c(Probe, GeneSymbol, NCBIid))
str(gmap)
```

We will next format the data matrix ("X" in AnnData)
```{r}
data <- data |> 
  column_to_rownames("Probe") |> 
  select(matches(gse)) |> 
  as.matrix()
str(data)
```

The data has weird sample names. We will clean those.
```{r}
colnames(data)
colnames(data) <- gsub("^.+?(GSM\\d+).+$", "\\1", colnames(data))
stopifnot(setequal(colnames(data), row.names(metadata)))
data <- data[,row.names(metadata)]
```
