

# Old Gemma commands (older versions of Gemma)

Let's load the required packages:
```
require(gemma.R)
require(tidyverse)
```

First, we will look for a dataset of interest. As an example, we will here look for datasets related to neuroblastoma. You can of course look for any type of topic you are interested in. 
```R
searchDatasets("neuroblastoma", limit = 100, taxon = "human") |>
  filter(geeq.batchCorrected == TRUE) |>
  select(ee.ShortName, ee.Name, ee.ID, ee.Accession, ee.Samples)
```

Next, pick on dataset and explore the description.
```R
gse <- "GSE21713"

getDatasetsInfo(gse) |>
  select(ee.ShortName, ee.Name, ee.ID, ee.Description)
```

Next, have a look at the metadata table for this dataset.
```R
metadata <- getDatasetDesign(gse)
str(metadata)
head(metadata)
with(metadata, table(batch, genotype))
```

Download the expression data.
```R
e <- getDatasetExpression(gse)
str(e)
colnames(e)
e <- as.data.frame(e)
```

The expression data is a data.frame / data.table object. We want to convert this to a matrix.
```R
# row.names of the metadata table are the sample names. Here we check whether they are all present in the expression matrix.
stopifnot(all(row.names(metadata) %in% colnames(e)))

# Next, let's get only the columns corresponding to sample names, make a matrix, and add gene symbols as row.names.
dataMT <- as.matrix(e[,row.names(metadata)])
str(dataMT)
row.names(dataMT) <- e$GeneSymbol
str(dataMT)
```

In the examples from day 2-5, we need to voom transform data (log2CPM). In GEMMA, this has already been done.
```R
boxplot(dataMT, las=2)
```