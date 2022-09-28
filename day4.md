# Interaction effects
In this exercise, we will test differences along two axes:
* interferon alpha and PBS
* spleen and liver
We will thus be able if there is an interaction between those two effects, i.e. whether the effect of interferon treatment (stimulated versus unstimulated) differs between the two organs spleen and liver.

## Setup
First load packages.
```R
.libPaths("/home/handson/R/x86_64-pc-linux-gnu-library/4.0/")
require(tidyverse)
require(limma)
require(patchwork)
require(pheatmap)
<!-- require(ComplexHeatmap)
require(enrichR) -->
```

Load the data.
```R
data <- readRDS("/home/handson/data/data.RDS")
design <- readRDS("/home/handson/data/design.RDS")
gmap <- readRDS("/home/handson/data/gmap.RDS")
```

## Subset data
We will only work with liver and spleen fibroblasts (Gp38 positive) that were treated with interferon alpha, and compare them to those cultivated only in phosphate buffered saline..
* Filter the design table accordingly
* Subset the data matrix by selecting only the columns that are in the filtered design table
* How many samples do we end up with?
* Use `?stopifnot` to make sure the data matrix has as many columns as the design table.


## Correlation analysis

### Correlation heatmap

* Use the correlation function in R `?cor` to correlate the samples in the data matrix. 
* Next, generate a heatmap of the resulting correlation heatmap using the function `?pheatmap`

### MDS projection
Finally, use MDS projection, modify the code from the last exercise, to also show the organ of each sample using the `shape` aesthetic in ggplot.
```R
data.frame(cmdscale(dist(2-corMT),eig=TRUE, k=2)$points) %>%
  add_column(stimulus = design$stimulus) %>%
  rownames_to_column("sample") %>%
  mutate(sn = gsub("^.+?_(\\d)$", "\\1", sample)) %>%
  ggplot(aes(x=X1,y=X2)) + 
  geom_point(aes(color=stimulus)) +
  geom_text(aes(label=sn)) +
  theme_bw()
```

## Differential expression and data normalization
In the next step we will compare stimulated to PBS control samples, liver to spleen, and the interaction of both effects.

### Setup up the model matrix
Make sure the correct references (PBS for stimulus and Liver. for the organ) are used by generating a heatmap of the model matrix. Then setup the model matrix as follows:
```R
model.matrix(~stimulus*organ, data=design)
```

### Normalize data
Now use limma voom to normalize the data.
```R
dataVoom <- voom(data, design=model, plot = TRUE) # insert your model matrix with design=model
```

### Perform differential expression
After having normalized the data we can fit the differential expression model. 
```R
limmaFit <- lmFit(dataVoom, design=model)
limmaFit <- eBayes(limmaFit)
```

Now let's look at which coefficients we get
```R
head(coef(limmaFit))
```

Next, we extract the results from these models.
```R
limmaRes <- list() # start an empty list
for(coefx in colnames(coef(limmaFit))){ # run a loop for each coefficient
	print(coefx)
  limmaRes[[coefx]] <- topTable(limmaFit, coef=coefx,number = Inf) # topTable returns the statistics of our genes. We then store the result of each coefficient in a list.
}
limmaRes <- bind_rows(limmaRes, .id = "coef") # bind_rows combines the results and stores the name of the coefficient in the column "coef"
limmaRes <- filter(limmaRes, coef != "(Intercept)") # then we keep all results except for the intercept
```

### Fit contrast
The above model has four coefficients:
* Intercept
* stimulusIFNa: interferon alpha versus PBS (based on liver samples)
* organSpleen: spleen versus liver (based on PBS samples)
* Interaction effect ("stimulusIFNa:organSpleen")
Now, we want to also quantify the IFNa stimulus effect in spleen. To do so we fit a contrast, specifically summing up the IFNa effect PLUS the interaction:
```R
colnames(coef(limmaFit))
stopifnot(all(colnames(coef(limmaFit)) == c("(Intercept)", "stimulusIFNa", "organSpleen", "stimulusIFNa:organSpleen"))) # make sure we have the right names, otherwise we have to adapt the next line
contrast.mt <- cbind(IFNa_Spleen = c(0,1,0,1)) # we add the 2nd and 4th coefficient.
row.names(contrast.mt) <- colnames(coef(limmaFit))
contrast.mt
limmaFit.contrast <- contrasts.fit(limmaFit,contrast.mt)
limmaFit.contrast <- eBayes(limmaFit.contrast)
limmaRes.contrast <- topTable(limmaFit.contrast, coef=colnames(contrast.mt),number = Inf) %>%
  rownames_to_column("ensg") %>%
  mutate(coef=colnames(contrast.mt))
limmaRes <- rbind(limmaRes.contrast, limmaRes) # add this coefficient to the result table
table(limmaRes$coef)
```

Now, we will clean up the table using regular expressions:
```R
limmaRes$gene <- gmap[limmaRes$ensg,]$external_gene_name # here we add the gene symbol
limmaRes <- limmaRes %>%
  mutate(coef = str_replace(coef, "organ", "")) %>% # remove "organ"
  mutate(coef = str_replace(coef, "stimulus", "")) %>% # remove "stimulus"
  mutate(coef = str_replace(coef, "^IFNa$", "IFNa_Liver")) %>% # rename "IFNa" to "IFNa_Liver"
  mutate(coef = str_replace(coef, "^IFNa\\:Spleen$", "Interaction")) # Name interaction
table(limmaRes$coef)
```

## Data interpretation
The steps below are identical in terms of code to the example from yesterday. 

### Vulcano plot
Draw a vulcano plot from the `limmaRes` object. Use `?facet_wrap` or `?facet_grid` to separate the plots by the stimulus (coefficient).

### P-value distribution
Draw a p-value distribution using `geom_histogram`, separate the plot using facets, and again look at the `AveExpr`.

### Number of hits
Now, count the number of genes that are tested `?count`. Then, create a new table `limmaResSig` where you retain only those genes that significantly change between conditions, thus filtering on the `adj.P.Val`. Consider also filtering lowly expressed genes based on the above plots (p-value distribution).


## Visualizing results
A key element of any statistical analysis is to visualize results (differential genes) to assess whether the statistics obtained match the data. 

### Visualizing one gene
* Pick one gene with significant interaction effects and a large absolute (negative or positive) log fold change from `limmaResSig`.
* Now create a table that we can use to plot this gene. To this end, modify the table `design` by adding the normalized expression of your gene of interest, taken from `dataVoom$E`, as a new column.
* Generate a plot, where the x-axis is the stimulus (IFNa or PBS) and the y-axis is the expression of the gene.
* Add facets to separate the two organs
* Look at the log fold changes for all coefficients. Do the observed differences on this plot fit to the log fold change?
* Note: You don't have to write the log fold change on the plot.

Example plot:
<img src="03_03_Interaction/One.gene.png" width="50%">

### Visualizing multiple genes
Now let's make the following plot, which shows the expression data (left) and the statistical results (right) for the top 5 genes from each comparison.
<img src="03_03_Interaction/Coef_HM.png" width="100%">

The steps below are outlined in detail. 

#### get the genes of interest
Based on the significant hits in `limmaResSig`, group (`?group_by`) the hits by the coefficient `coef`, then get the top 5 genes by logFC (`top_n`), extract the ENSEMBL IDs from the column `ensg` using `?pull`, and store the result in a new object `goi.all`. 

#### plot statistical results
Next plot all statistical results for the genes above. This plot is the same as yesterday.
```R
(p.coef <- limmaRes %>%
  filter(ensg %in% goi.all) %>%
  mutate(gene = gmap[ensg,]$external_gene_name) %>%
  ggplot(aes(y=gene, x=gsub("stimulus", "", coef), color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw())
```

#### plot expression data
First we will collect the expression data of each gene, writing a for loop over all genes, and storing the data.frame for each gene in a list.
```R
dat.list <- list()
for(gg in goi.all){
  dat.list[[gg]] <- design %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample") %>%
    remove_rownames()
}
```

Next, we combine the above list of data.frame into one data.frame using `?bind_rows`, and then plot this data as a heatmap. The code below is the same as yesterday. Now, also separate samples by the organ.
```R
(p.vals <- bind_rows(dat.list, .id="ensg") %>%
  mutate(gene = gmap[ensg,]$external_gene_name) %>%
  mutate(stimulus = as.character(stimulus)) %>%
  ggplot(aes(x=sample, y=gene, fill=E)) + 
  geom_tile() +
  facet_grid(. ~ stimulus, space ="free", scales = "free") +
  scale_fill_gradient2(low="blue", high="red"))
```

#### Final plot
Finally, we combine the two plots as below, using the "patchwork" package. This command should show you the plot shown above.
```R
p.vals + p.coef
```