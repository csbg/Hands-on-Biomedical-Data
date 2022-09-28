# Introduction in differential expression analysis


## Setup
First load packages.
```R
require(tidyverse)
require(limma)
require(ComplexHeatmap)
require(enrichR)
```

Then load the data.
```R
data <- readRDS("data.RDS")
design <- readRDS("design.RDS")
gmap <- readRDS("gmap.RDS")
```

For all three objects above, answer the following questions:
* What type of object is it?
* How many rows and columns are this object?
* What information is contained in rows and columns?

## Subset data
In today's exercise, we will only work with liver fibroblasts (Gp38 positive) that were treated with interferon alpha, and compare them to those cultivated only in phosphate buffered saline. To subset the dataset accordingly we need the following steps:
* Filter the design table accordingly
* Subset the data matrix by selecting only the columns that are in the filtered design table

After subsetting, the design table should only contain 6 rows and the data matrix only 6 columns.

## Correlation analysis

### Correlation heatmap

* Use the correlation function in R `?cor` to correlate the samples in the data matrix. 
* Next, generate a heatmap of the resulting correlation heatmap using the function `?Heatmap`

### MDS projection
Finally, we will use the calculated correlations to project the samples on 2 dimensions. The entire code for this step is shown below. This will:
* Transform correlations into distance measures by calculating `2-correlation`
* Run the function `?cmdscale` to get a 2-dimensional projection
* Modify the table to add information relevant for plotting
* Plot the samples on two dimensions
* --> don't forget: You can execute parts of the code to better understand what it does!
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
In the next step we will compare interferon-treated to PBS control samples.

### Setup up the model matrix
The model matrix, also called "design matrix", defines which group will be set as the intercept and which comparisons will be performed. In R, the function `?model.matrix` is used for this purpose. In our case, we only compare stimulated to unstimulated samples, so we can use:
```R
model.matrix(~stimulus, data=design)
```
* Make a heatmap of the resulting model.matrix
* Now figure out which condition is taken as the control / reference / intercept. This should be PBS. 
* If the reference is not the right one, use `factor`, `?relevel`, and `?mutate` to change the  factor levels.

### Normalize data
After defining the design matrix, we can use limma voom to normalize the data.
```R
dataVoom <- voom(data, design=model, plot = TRUE) # insert your model matrix with design=model
```

Now let's look at the data before and after normalization. The original data is in the object `data`, the normalized data is in `dataVoom$E`.
* What types of objects are those two?
* Use `?boxplot` to plot the distributions of the first few (30) genes of the original matrix. Try both of the following approaches. What's the difference?
	* `boxplot(data[1:30,])`
	* `boxplot(t(data[1:30,]))`
* Now use `dataVoom$E` to plot a boxplot of the first 30 genes.
* Now lot's look at the density of one gene using three approaches. What's the difference?
	* plot(density(data[8,]))
	* plot(density(dataVoom$E[8,]))
	* plot(density(log2(data[8,])))

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
  limmaRes[[coefx]] <- topTable(limmaFit, coef=coefx,number = Inf) # topTable returns the statistics of our genes. We then store the result of each coefficient in a list.
}
limmaRes <- bind_rows(limmaRes, .id = "coef") # bind_rows combines the results and stores the name of the coefficient in the column "coef"
limmaRes <- filter(limmaRes, coef != "(Intercept)") # then we keep all results except for the intercept
```

### Analyze results
