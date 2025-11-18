# A more complex example
In this exercise, we will repeat the analysis from the introductory exercise, but analyze all six stimuli compared to PBS controls.

**NOTE:** In this exercise, make sure to **split the genes between the different stimuli** when reporting numbers of genes, performing, enrichment, etc. It can be interesting to report which gene is regulated by *any* stimulus but ultimately we want to compare the genes affected by *each* stimulus!

## Setup
First load packages.
```R
require(tidyverse)
require(limma)
require(patchwork)
require(ComplexHeatmap)
require(msigdbr)
```

Load the data.
```R
data <- readRDS("data.RDS")
metadata <- readRDS("design.RDS")
gmap <- readRDS("gmap.RDS")
MSigDB <- readRDS("MSigDB.rds")
```

## Subset data
We will only work with liver fibroblasts (Gp38 positive) but use all stimuli.
* Filter the metadata table accordingly
* Subset the data matrix by selecting only the columns that are in the filtered metadata table
* Use `?stopifnot` to make sure the data matrix has as many columns as the metadata table has rows.

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 3.1:`
How many samples do we end up with?

## Correlation analysis

### Correlation heatmap
Generate the heatmap of the correlation matrix as discussed on day 2.

### MDS projection
Finally, use MDS projection, using the same code from the last exercise:
```R
data.frame(cmdscale(dist(2-corMT),eig=TRUE, k=2)$points) |>
  add_column(stimulus = metadata$stimulus) |>
  rownames_to_column("sample") |>
  mutate(sn = str_replace(sample, "^.+?_(\\d)$", "\\1")) |> # This shortens the sample names to just the number at the end
  ggplot(aes(x=X1,y=X2)) + 
  geom_point(aes(color=stimulus)) +
  geom_text(aes(label=sn)) +
  theme_bw()
```

## Differential expression and data normalization
In the next step we will compare stimulated to PBS control samples.

### Setup up the model matrix
Just like yesterday, we compare stimulated to unstimulated samples. Make sure the correct reference (PBS) is used by generating a heatmap of the model matrix.
```R
model.matrix(~stimulus, data=metadata)
```

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 3.2:`
Provide the heatmap of the model matrix. Which coefficient from the model matrix compares which groups?

### Normalize data
Now use limma voom to normalize the data.
```R
dataVoom <- voom(data, design=your_model_matrix, plot = TRUE) # insert your model matrix
```

### Perform differential expression
After having normalized the data we can fit the differential expression model. 
```R
limmaFit <- lmFit(dataVoom, design=your_model_matrix)
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
	# topTable returns the statistics of our genes. We then store the result of each coefficient in a list.
	# The rownames (ENSEMBL Gene IDs) are stored in the column with the name "ensg"
  limmaRes[[coefx]] <- topTable(limmaFit, coef=coefx,number = Inf) |>
		rownames_to_column("ensg")
}
limmaRes <- bind_rows(limmaRes, .id = "coef") # bind_rows combines the results and stores the name of the coefficient in the column "coef"
limmaRes <- filter(limmaRes, coef != "(Intercept)") # then we keep all results except for the intercept
```

## Data interpretation

### Vulcano plot
Draw a vulcano plot from the `limmaRes` object.
Use `?facet_wrap` or `?facet_grid` to separate the plots by the stimulus (coefficient).

### P-value distribution
Draw a p-value distribution using `geom_histogram`, separate the plot using facets, and again look at the `AveExpr`.

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 3.3:`
Provide the vulcano plots and p-value distributions in your protocol.

### Number of hits
Now, count the number of genes that are tested using `?count`. Then, create a new table `limmaResSig` where you retain only those genes that significantly change between conditions, thus filtering on the `adj.P.Val`. Consider also filtering lowly expressed genes based on the above plots (p-value distribution).

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 3.4:`
Report the number of tested and significant genes for each comparison. Report the code (should just be a few lines) you use to answer the questions.


## Visualizing results
A key element of any statistical analysis is to visualize results (differential genes) to assess whether the statistics obtained match the data. 

### Visualizing one gene

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 3.5:`
* Pick one gene from one comparison with significant effects and a large absolute (negative or positive) log fold change from `limmaResSig`.
* Explain / report your code on how you picked the gene.
* Now create a table that we can use to plot this gene. To this end, modify the table `metadata` by adding the normalized expression of your gene of interest, taken from `dataVoom$E`, as a new column.
* Generate a plot, where the x-axis is the stimulus (six stimuli and PBS) and the y-axis is the expression of the gene.
* Look at the log fold changes for all six stimuli. Do the observed differences on this plot fit to the log fold change?
* Note: You don't have to write the log fold changes on the plot but you should provide them in your protocol.

Example plot:
<img src="03_02_Complex/One.gene.png" width="50%">

### Visualizing multiple genes
Now let's make the following plot, which shows the expression data (left) and the statistical results (right) for the top 5 genes from each comparison.
<img src="03_02_Complex/Coef_HM.png" width="100%">

The steps below are outlined in detail. Make sure you understand the code, as you will have to modify it tomorrow.

#### get the genes of interest
Based on the significant hits in `limmaResSig`, group (`?group_by`) the hits by the coefficient `coef`, then get the top 5 genes by logFC (`slice_max()`), extract the ENSEMBL IDs from the column `ensg` using `?pull`, and store the result in a new object `goi.all`. What is the content of `goi.all`?

#### plot statistical results
Next plot all statistical results for the genes above.
```R
(p.coef <- limmaRes |>
  filter(ensg %in% goi.all) |>
  mutate(gene = gmap[ensg,]$external_gene_name) |>
  ggplot(aes(y=gene, x=str_remove(coef, "stimulus"), color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw())
```

#### plot expression data
First we will collect the expression data of each gene, writing a for loop over all genes, and storing the data.frame for each gene in a list.
```R
dat.list <- list()
for(gg in goi.all){
  dat.list[[gg]] <- metadata |>
    mutate(E=scale(dataVoom$E[gg,])) |>
    rownames_to_column("sample") |>
    remove_rownames()
}
```

Next, we combine the above list of data.frame into one data.frame using `?bind_rows`, and then plot this data as a heatmap.
```R
(p.vals <- bind_rows(dat.list, .id="ensg") |>
  mutate(gene = gmap[ensg,]$external_gene_name) |>
  mutate(stimulus = as.character(stimulus)) |>
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

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 3.6:`
Report this plot in your protocol.

## Enrichment analysis
Enrichment analysis help in interpreting long lists of genes. By measuring whether certain gene sets are enriched in our list of differential genes (often called hit list), enrichment analysis informs us on the involvement of biological pathways (among others) in the processes studied.

#### Perform enrichment analysis for each coefficient
Below is a loop over the individual coefficients (comparisons). Within each iteration of this loop, we will perform enrichment analysis for all genes significant in each coefficient. 

As a reminder, this were the instructions from yesterday:
* First, filter all genes with `logFC > 0` from the table of significant genes and store them in the object `goi` (note, this will overwrite the value of this object defined previously - so if you are going back to the previous exercise, you wil have to redefine the object).
* Next convert the ENSEMBL IDs to gene symbols: `goi <- gmap[goi,]$external_gene_name |> unique()`
* Next prepare the background/universe set
```R
universe <- gmap$gene_unique[match(limmaRes$ensg, rownames(gmap))] |> unique()
universe <- unique(universe[!is.na(universe) & universe != ""])
```
* Specify the database of genesets
You already downloaded MSigDB 
* Next perform enrichment analysis using the following function and store the results in the objec `fisher_tbl`.
Fisher's exact test is performed against each pathway in the MSigDB database, 

Now you will run this in a loop, with one iteration for each comparison (coefficient)
* Note: you will have to use the variable `coefx` INSIDE of the loop to get the right genes for each iteration.
* Look at the value of `unique(limmaResSig$coef)`
* Run the loop below once to see what it does. Then you have to edit the code as described to make it work

``` R
fisher_list <- list()
for(coefx in unique(limmaResSig$coef)){
    
    # What is coefx now?
    print(coefx)

	# Extract genes of interests (GOI) for a given coefficient (see yesterday's example)
	# goi <- .... # YOUR INPUT NEEDED HERE (remember to provide gene names using gmap)
	
	# Add code here to perform enrichment analysis (see yesterday's example)
	fisher_tbl <- map_df(names(MSigDB), function(pw) {
 	 pw_genes <- MSigDB[[pw]]
     sig_genes <- goi
     overlap <- intersect(sig_genes, pw_genes)
		a <- length(overlap)
		b <- length(sig_genes) - a
  		c <- length(pw_genes) - a
  		d <- length(universe) - (a+b+c)
  
  	ft <- fisher.test(matrix(c(a,b,c,d),2), alternative = "greater")
  
  	tibble(
    	pathway = pw,
    	pvalue = ft$p.value,
    	odds_ratio = as.numeric(ft$estimate),
    	DE_in_pathway = a,
    	overlapping_genes = paste(overlap, collapse = ",")
  	)
	})
	fisher_tbl <- fisher_tbl |> mutate(p_adj = p.adjust(pvalue, "BH"))

  # Store results in the list
  # fisher_list[[coefx]] <- fisher_tbl # JUST UNCOMMENT HERE
}
```

Finally, we combine the list (each entry is one coefficient) into one long table:
```R
fisher_tbl <- bind_rows(fisher_list, .id="coef")
```

#### Plot enrichments

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 3.7:`
Now generate the following plot:

<img width="1880" height="1950" alt="top_10_pathways" src="https://github.com/user-attachments/assets/574d614a-f286-4b1e-90fa-866b6980eaa9" />


Note: The plot only includes entries with: p_adj < 0.05 & odds_ratio > 6

#### Plot genes related to the enrichments
Now we will extract the genes underlying the above enrichments:
```
goi.enr <- fisher_tbl |>
  filter(p_adj < 0.05 & odds_ratio > 6) |>
  pull("overlapping_genes") |>
  str_split(",") |> 
  unlist() |>
  unique()
```

Then we will extract the statistics for these genes:
```R
limmaRes |>
  mutate(gene = gmap[ensg,]$external_gene_name) |>
  filter(gene %in% goi.enr)
```

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 3.8:`
From this table, generate the following plot and discuss whether the pathways are expected.

<img src="03_02_Complex/Enrichments.genes.png" width="50%" height="100%">


## Final questions
Answer the questions below by again showing the relevant plots from the exercises above.

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 3.9:`
Looking at the correlation heatmap and MDS plot - do you see strong effects and clear differences between groups?

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 3.10:`
Which treatment has the largest effects? How do you judge this? Which treatment effects most genes? Which treatment has the largest (positive or negative) log fold changes? Are those two consistent?
