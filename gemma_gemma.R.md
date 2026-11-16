#For the new version of gemma, please use this code, update gse <- "dataset number" . Example gse <- "GSE1232"

 
<!-- 
In the following, we will look for interesting datasets to use. The code below depends on the gemma.R versions. If this fails, you can also identify a relevant dataset by going to the [website](https://gemma.msl.ubc.ca/browse/#/). Then skip this subsection (`explore datasets`) and continue below to obtain the data.

To explore datasets, we will look for a dataset of interest. As an example, we will here look for datasets related to neuroblastoma. You can of course look for any type of topic you are interested in.
```R
get_datasets("neuroblastoma", limit = 100, taxa = "human") |>
  filter(geeq.batchCorrected == TRUE) |>
  select(taxon.Name, taxon.ID, experiment.Accession, experiment.SampleCount)

```

The above commands are for `gemma.R` version 2.0.0 and can be different for your version. You can find out the version of gemma.R using the function `sessionInfo()` or `packageVersion("gemma.R")`. If you run the following command, it will tell you the column names (Gemma fields) in your version. You can then use the right ones instead of `taxon.Name`, `taxon.IDF`, `experiment.Accession`, and `experiment.SampleCount` above.
```R
colnames(get_datasets("neuroblastoma", limit = 100, taxa = "human"))
```

In particular, the different fields that describe a dataset (for example: `geeq.batchCorrected`, `taxon.Name`) may not have the same names in your versions. If you use an older version of `gemma.R`, see [here](day1_OldGemma.md) for older syntax. For newer version, see the vignettes for gemma versions [3.0.14](https://bioconductor.org/packages/3.19/bioc/html/gemma.R.html) and [3.2.0](https://bioconductor.org/packages/3.20/bioc/html/gemma.R.html). Go to the right website and click on the vignette for `Accessing curated gene expression data with gemma.R` to see examples how to run this in your version of `gemma.R`.

Next, pick one dataset and explore the description.
```R
gse <- "GSE21713"

get_datasets(gse)|>
  select(experiment.ShortName, experiment.Name, experiment.ID, experiment.Description)
```

Use column names accordingly.
```R
get_datasets(gse)|>
dplyr::select(experiment.shortName, experiment.name, experiment.ID, experiment.description)
``` -->




# Using Gemma.R

We can use Gemma.R, which enables us to easily download datasets from GEMMA. `gemma.R` is a package that enables us to download data and metadata for many expression datasets. Unfortunately, it is not as well developed as tidyverse and ggplot2, it may thus be a bit difficult to use. Ask for help if needed and - if `gemma.R` doesn't work - answer the finale exercise from today based on what you can do on the GEMMA website.

First, load the required packages:
```R
require(gemma.R)
require(tidyverse)
```

Store the dataset number in gse. Here we will use the following example:
```R
gse <- "GSE21713"
```

Now, download the metadata for this dataset. 
```R
metadata <- get_dataset_samples(gse) |>
  make_design('text') |> 
  select(-factorValues)
```

If the above function does not work, try the following from an older gemma version.
```R
metadata <- get_dataset_design(gse)
```

Now, let's explore this metadata.
```R
str(metadata)
head(metadata)
with(metadata, table(block, genotype))
```

Download the expression data.
```R
e <- get_dataset_processed_expression(gse)
str(e)
colnames(e)
e <- as.data.frame(e)
```

The expression data is a data.frame / data.table object. We want to convert this to a matrix. Rownames of the metadata table are the sample names. Here we check whether they are all present in the expression matrix.
```R
stopifnot(all(row.names(metadata) %in% colnames(e)))
```

Next, let's get only the columns corresponding to sample names, make a matrix, and add gene symbols as row.names.
```R
dataMT <- as.matrix(e[,row.names(metadata)])
str(dataMT)
row.names(dataMT) <- e$GeneSymbol
str(dataMT)
```

In the examples from day 2-5, we need to voom transform data (log2CPM). In GEMMA, this has already been done.
```R
boxplot(dataMT, las=2)
```


## Naming issues?

Note: If the command `stopifnot(all(row.names(metadata) %in% colnames(e)))` fails, it could be that the names slightly differ between the expression data matrix and the modelmatrix. You then may need to transform these names beforehand. The following should be sufficient:

For the metadata `metadata`, we run the function `make.names()`, which removes uncommon symbols.
```R
row.names(metadata) # This will show you the current row.names
make.names(row.names(metadata)) # This will "clean" the row.names
row.names(metadata) <- make.names(row.names(metadata)) # This will overwrite the row.names with the "cleaned" ones
row.names(metadata) # Now this should show you the clean row.names
``` 

For the expression data matrix `e`, we first remove empty spaces `" "` with the function `gsub()` and then run `make.names()`.
```R
colnames(e) # original colnames
gsub(" ", "", colnames(e))  # removed spaces (if there were any)
make.names(gsub(" ", "", colnames(e))) # otherwise "clean" colnames
colnames(e) <- make.names(gsub(" ", "", colnames(e))) # Overwrite
colnames(e) # clean names
``` 

* You don't need to run the above if the names already match.
* If the names are well designed (only use normal letters and digits) the above code will not change them and no change is required
* Use the `stopifnot()` statement after renaming to make sure the issue is solved
* The `gsub()` command is not required for the row names of `metadata` because row names of a `data.frame` are not allowed to have spaces
