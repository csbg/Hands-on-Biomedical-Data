# Assignment - analyze a dataset of your choosing

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 5:` All of the results here should be reported in your protocol.

## Gemma
* As seen on day 1, the Gemma database contains expression datasets with curated expression and meta-data. 
* The database can also be browsed using this [link](https://gemma.msl.ubc.ca/).
* Under this [link](https://gemma.msl.ubc.ca/expressionExperiment/showAllExpressionExperiments.html), we can search the available expression datasets.
* To access datasets with R, see the [Gemma Tutorial on day1](day1.md).
* Please select datasets with > 10 samples and consider those where you can include interaction effects
  

<img src="img/Gemma.png" width="100%">

## Instructions
* For the final exercise of this (R) part of the practical, explore a dataset from Gemma using the function taught in this course.
* Be aware that based on the metadata table of the dataset, you will have to modify the model matrix (`model.matrix`) accordingly to correctly interpret the coefficients (for example, think about batch and interaction coefficients in addition to the main comparisons).
* Also - make sure you define an appropriate intercept term.
* While it is more elegant to analyze the full dataset, you can also subset the dataset if it is too large and/or complex.
  

## Exercise details:

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 5.1` (2 point)
How many samples do you have. What are the conditions? Will youfilter the dataset and focus on a subset of data?

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 5.2` (3 points)
Plot correlation heatmaps and MDS plots. Discuss the signal to noise ratio.

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 5.3` (3 points)
Report your model matrix (heatmap) and describe which conditions you are going to compare.

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 5.4` (3 points)
Provide your vulcano plot and p-value distribution. Do you see significant differences? In which conditions (if you have many). How many significant genes do you find?

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 5.5` (3 points)
Plot the (normalized) data for one and more significant genes and discuss whether the data fits to the statistics (log fold change and p-value) returned by Limma.

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 5.6` (3 points)
Perform enrichment analysis. Are the results expected? (If you do not have significant genes, use the 50 or 100 genes with the biggest effects (based on log fold changes).)

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 5.7` (3 points)
Discuss 3 significant genes (or those with large changes) in the context of the dataset studied. Plot the data for those genes and plot or report their statistics (log fold change and p-value).
