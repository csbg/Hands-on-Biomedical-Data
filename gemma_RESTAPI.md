# Downloading data from GEMMA using the REST-API

GEMMA provides approaches to download data and metadata for any dataset from the webbrowser.
This is called [representational state transfer (REST) application programming interface (API)](https://gemma.msl.ubc.ca/resources/restapidocs/#/).
As part of this, we will create URLs (Uniform Resource Locator). For example, for the dataset `GSE2018`:
* This URL downloads the design (metadata): `https://gemma.msl.ubc.ca/rest/v2/datasets/GSE2018/design`
* This URL downloads the expression data: `https://gemma.msl.ubc.ca/rest/v2/datasets/GSE2018/data/processed`

Note: You can use these URLs in your web browser but to obtain the data in the correct format, we will use R to downlad the files (see below).

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
gse <- "GSE21713"
```

## Download

Download metadata of the URL `https://gemma.msl.ubc.ca/rest/v2/datasets/GSE2018/design`.
```{r}
url.metadata <- paste(gemma.url, gse, "design", sep="/")
response <- GET(url.metadata, add_headers(Accept = "text/tab-separated-values"))
writeBin(content(response, "raw"), paste0(gse, "_metadata.tsv"))
```

Download data of the URL `https://gemma.msl.ubc.ca/rest/v2/datasets/GSE2018/data/processed`
```{r}
url.data <- paste(gemma.url, gse, "data/processed", sep="/")
response <- GET(url.data)
writeBin(content(response, "raw"), paste0(gse, "_data.tsv"))
```

**Note:** This will download "processed" data, i.e. they are not counts and already normalized to log2CPM values.

## Parse metadata

First the metadata ("obs" in AnnData).
```{r}
metadata <- read_tsv(paste0(gse, "_metadata.tsv"), comment = "#")
head(metadata)
```

Summarize the experiment
```{r}
metadata |> 
  count(genotype, batch)
```

## Parse data

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

## Map metadata to data

Check if the names the same in metadata and data
```{r}
row.names(metadata)
colnames(data)
metadata$Bioassay
setdiff(colnames(data), metadata$Bioassay) # checking that all sample names in the data are also in the metadata
```

Let's make some names that work in both.
```{r}
metadata <- metadata |> 
   mutate(ID = gsub("^.+(Biomat_\\d+).+$", "\\1", Bioassay)) |> # extracts the commmon name from both
   remove_rownames() |>
   column_to_rownames("ID") # First for the metadata
print(row.names(metadata))

colnames(data) <- gsub("^.+(Biomat_\\d+).+$", "\\1", colnames(data)) # Then for the data
colnames(data)

setdiff(colnames(data), row.names(metadata)) # checking again
```

Make sure this worked out.
```{r}
stopifnot(setequal(colnames(data), row.names(metadata)))
```

And finally, we order the columns in data according to the metadata. This ensures the same samples in both and the same order.
```{r}
data <- data[,row.names(metadata)]
```

## Make nicer sample names

The sample names above are not informative. Here we create our own sample names.
```{r}
metadata <- metadata |> 
  mutate(genotype2 = genotype) |> 
  mutate(genotype2 = str_replace(genotype2, "MYCN.+", "MYCN")) |> 
  mutate(sample_name = str_replace(genotype2, "wild.+", "WT")) |> 
  mutate(sample_name = make.unique(paste0(genotype2, "_", batch))) |> 
  mutate(sample_name = str_replace(sample_name, "wild type genotype", "WT")) |> 
  remove_rownames() |>
  mutate(ID = sample_name) |> 
  column_to_rownames("ID")
```

Change data sample names (in columns)
```{r}
colnames(data) <- row.names(metadata)
```

## Use readable gene names
```{r}
gmap <- gmap  |> 
  mutate(GeneSymbol_u = make.unique(GeneSymbol))   |> 
  column_to_rownames("Probe")
row.names(data) <- gmap[row.names(data), "GeneSymbol_u"]
```
