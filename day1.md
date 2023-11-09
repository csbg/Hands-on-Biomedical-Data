# R programming and visualizations
In this part, we will get started with R and Rstudio. Throughout this course, you will run R commands from these files in Rstudio. You can type them but you can also copy/paste them. 

If you run RStudio in a browser (RICARDA or Galaxy), one way is to open this window on side of your screen, and open you R studio in a different window on the other side of the screen.
If you can use another browser, you may be able to use `Strg + TAB` to switch bettween the two windows using your keyboard (and don't have to click between tabs).

### R studio basics
<img src="img/RStudio.png" width="100%">

For this part of the practical (next 5 days), you will be writing R code in R scripts:
* When you log in you have started an `R session`
* When you execute commands in the  console (bottom left) the output and potential error or warning messages will appear in the console 
* When you also create variables, they get stored in your `environment` (see top right)
* You can restart your session by quitting the session (top right, red button)
* You can also clear objects from your environment (also called "workspace") by clicking the small yellow broom icon on the top right under "Environment"
* Plots and help packages appear in the bottom right. There is also a file browser that you can use to open files.
* We will almost always work with R scripts (top left) and use `Strg + ENTER` to "send" individual commands or sections of commands (marking them with the mouse) to the console, where they will be executed
* ***While you can execute commands from your script in any order, make sure your script runs through from top to bottom if started from scratch - so do not use variables before you define them!***

## Loading packages

Once packages are installed, you only need to load them in your current session, load the following packages:
```R
require(tidyverse)
require(pheatmap)
```

## Variables

Now let's explore some basic variables
```R
a <- 1
b <- 3
typeof(a)
typeof(b)
str(a)
str(b)
a
b
```

## Functions
The sum function can be applied to numbers 
```R
sum(a,b)
a + b
sum(5,6)
c <- sum(a,b)
(c <- sum(a,b))
str(c)
sum(a,c)
```

The sum function cannot be applied to characters. The following code will create an error. We thus put `#+ error=TRUE` in the field so that the generation of the markdown file does not get interrupted.
```R
#+ error=TRUE
a <- "abc"
b <- "hello"
sum(a,b)
typeof(a)
typeof(b)
str(a)
str(b)
```

Caution: Some data structures can be automatically interpreted
```R
a <- TRUE
b <- FALSE
sum(a,b)
sum(a,a)
typeof(a)
typeof(b)
str(a)
str(b)
```

## Vectors
```R
vec <- c(1,2,3,4,1,2,3)
names(vec) <- LETTERS[1:4]
str(vec)
sum(vec)
mean(vec)
vec[2]
vec[3]
vec["A"]
vec[c("A", "D")]
unique(vec)
1:10
-5:5
vec <- 10:1
vec[2:3]
```

Now create a vector of the letters B,C,D by replacing `???` with the right code. The `stopifnot()` statement will test if this worked correctly.
```R
your.vector <- LETTERS[???]
stopifnot(your.vector == c("B", "C", "D"))
```

## Lists
```R
list_x <- list("a", 1, "b", "xyz", TRUE)
str(list_x)
list_x[[1]]
str(list_x[[1]])
str(list_x[1])
list_x[[2]]
list_x[2:4]
str(list_x[2:4])

vec_x <- c("a", 1, "b", "xyz", TRUE)
#' R has converted all non-characters to characters!
str(vec_x) 
```

## Loops and conditions
```R
for(x in 1:6){
  print(paste("x =", x))
  if(x > 3){
    print("...x is greater than three")  
  }
  if(x == 2){
    print("...x equals two")
  }
  if(x != 4){
    print("...x is not four")
  }
  if(x %% 2 == 0){
    print("...x is even")
  } else {
    print("...x is uneven")
  }
  if(x %in% c(3,5)){
    print("...x is three or five")
  }
  print("-------")
}
```

## Matrices
For this exercise, first load the 'data.RDS'. See [setup instructions](README.md) for details.
```R
m <- readRDS("data.RDS") # you probably will have to modify this line
```

The following functions summarize the matrix
```R 
str(m)
dim(m)
head(m)
row.names(m)
colnames(m)
```

We can subset the matrix by defining the rows or the columns: `matrix[ROWS,COLUMNS]`. The following code yields the first 20 rows of the matrix. We create a vector of numbers 1 to 20 `1:20` and then get those rows from the matrix `m[1:20,]`.
```R
1:20
dim(m)
dim(m[1:20,])
```

The next code yields all columns that match the string `"Liver_Fibroblasts"`. We first get a boolean `TRUE` and `FALSE` vector matching the string to the column names, and then get those columns from the matrix.
```R
colnames(m)
grepl("Liver_Fibroblasts", colnames(m))
dim(m)
dim(m[, grepl("Liver_Fibroblasts", colnames(m))])
```

Now we will subset the matrix in both the rows and columns:
```R
m <- m[1:20, grepl("Liver_Fibroblasts", colnames(m))]
```

Then we will rename the columns and the rows.
```R
colnames(m) <- gsub("^Liver_Fibroblasts_(.+)_RNA_(\\d)$", "\\1_\\2", colnames(m))
row.names(m) <- paste0("g", 1:nrow(m))
str(m)
m
```

The `t()` command transposes the matrix (switches rows and columns).
```R
t(m)
dim(m)
dim(t(m))
```

To calculate correlation, we only need to run the function `cor()`. Next, we will make a heatmap of the data.
```R
cor(m, method="spearman")
pheatmap(cor(m, method="spearman"))
```

This looks nicer if the diagnole is converted to NAs.
```R
cMT <- cor(m, method="spearman")
diag(cMT) <- NA
pheatmap(cMT)
```

Now subset the matrix to the first 30 rows and the first 10 columns by replacing `???` for the correct code. The `stopifnot()` statement will test if this worked correctly.
```R
matrix.check <- m[???,???]
dim(matrix.check)
stopifnot(all(dim(matrix.check) == c(30,10)))
```


## Data frames and dplyr

We will work with the `starwars` dataset, an example dataset of starwars characters. This is not a `data.frame` but a `tibble`, but they can be treated using similar functions, and we will convert the table to a data.frame later.
```R
starwars
str(starwars$name)
```

Useing `stopifnot()` to see if the results are what we expect. Here we test whether each name occurs only once, if the number of unique names is the same as the number of rows of the table.
```R
stopifnot(length(unique(starwars$name))==nrow(starwars))
```

Let's keep only columns that are not of type `list`, then we can convert the table to a `data.frame`:
```R
sw <- starwars %>%
  select(where(function(x) !is.list(x))) %>%
  as.data.frame()
```

Use the function `count()` to count the number of rows (characters) by their `homeworld`.
```R
sw %>% count(homeworld)
```

Now we will add a new column called `firstname`
```R
sw <- sw %>% 
  mutate(firstname = str_remove(name, " .+$")) 
```

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 1.1:`
Assess the gender balance of this table, using `count()` and `gender`.

Count the number of characters by their `skin_color`. Next, run this code:
```R
sw %>%
  pull("skin_color") %>%
  str_split(", ") %>%
  unlist() %>%
  table() %>%
  sort()
```
* `pull` extracts a column from the `data.frame`
* `str_split` splits strings (characters) by some character
* `unlist` transforms a list into a vector
* `table` is simlar to `count` and counts the occurances
* `sort` orders a vector

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 1.2:`
What is the difference between the two approaches to count by skin color?

Print the names of everyone over 2m (height greater than 200) by fixing the following code (replace `???` with the correct code).
```R
sw %>%
  filter(???) %>%
  pull(???)
```

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 1.3:`
Who is taller than 2m?

Now let's take a detailed look at filters and conditions.
```R
starwars %>%
	filter(hair_color == "blond") %>%
	filter(sex == "male")

starwars %>%
	filter(hair_color == "blond" & sex == "male")
	
starwars %>%
	filter(hair_color == "blond" | sex == "male")
```
![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 1.5:`
What is the difference between these 3 approaches?

In the data.frame `sw`, use the function `mutate()` to calculate the body mass index (BMI) for all characters using the formula `mass/(height/100)^2`, storing this in the new column `bmi` .

Look at the distribution of these values. Note, below are three ways of extracting a column from a `data.frame`
```R
# The dollar sign can be used with names
quantile(sw$bmi, na.rm = TRUE)

# data.frames are lists, thus the same syntax with two squared brackets works
quantile(sw[["bmi"]], na.rm = TRUE)

# the function 'pull' can also extract columns
quantile(pull(sw, "bmi"), na.rm = TRUE)
```

## Plotting
We will use ggplot to visualize results. The following is a basic plot.
```R
ggplot(sw, aes(x=mass, y=height)) + geom_point()
```

You can also store a plot in a variable. 
```R
px <- ggplot(sw, aes(x=mass, y=height))
```

Note that here the `geom_point()` was not called. So the following is an empty plot:
```R
px
```

Add points to this plot:
```R
px + geom_point()
```

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 1.5:`
Insert the obtained plot and the following modifications (until the next section) in your protocol.

Add `geom_hex()` to the plot.

In addition to the points above, add labels to the points:
* Add `geom_text()` to the plot, in addition to `geom_point()`. To do so, you will have to choose which text will be displayed. Do so by using `firstname` in the `label` aesthetic.
* Modify to only plot `geom_text()` for characters with mass greater than 1000. You can do this by setting the `data` parameter in `geom_text`
* Modify the `color` and `shape` in `geom_point`

Use `geom_histogram()` to plot a histogram of the `height`. Note: you only need an `x` aesthetic.

Modify the plot to show the density `geom_density()` and the empirical cumulative distribution function `stat_ecdf()` instead.

Now let's plot the BMI of some individuals:
```R
sw %>% 
  group_by(gender) %>%
  top_n(3, bmi) %>%
  ggplot(aes(x=fct_reorder(firstname, bmi), y=bmi, fill=gender)) + 
  geom_bar(stat="identity")
```

Modify the above plot to add facets with `facet_grid()`. There is a lot of white space. Remove white space by change the `space` and `scales` parameters.


## Factors

We will make a simple plot to compare the hight between sexes:
```R
ggplot(sw, aes(x=sex,y =height)) + geom_violin()
```

What if you only want to show males and females? Filter accordingly. Note: the `%in%` statement tests whether the values in the column `sex` are *in* the vector `c("male", "female")`. This is similar to using an "or" statement with `|` like this `sex == "male" | sex == "female"`.
```R
sw %>%
  filter(sex %in% c("male", "female")) %>%
  ggplot(aes(x=sex,y =height)) + geom_violin()
```

Now, assume we want to show the males first. For this, we will need to use a factor. Factors are vectors with an order, so-called *levels*. 
```R
str(sw$sex)

sw %>%
	filter(sex %in% c("male", "female")) %>%
	pull(sex) %>%
	str()

sw %>%
	filter(sex %in% c("male", "female")) %>%
	mutate(sex = factor(sex)) %>%
	pull(sex) %>%
	str()

sw %>%
	filter(sex %in% c("male", "female")) %>%
	mutate(sex = factor(sex, levels=c("male", "female"))) %>%
	pull(sex) %>%
	str()
```

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 1.6:`
Update the code below to show the violin plot with the males in the first violin and the females in the second one. Copy it in your protocol.
```R
sw %>%
  filter(sex %in% c("male", "female")) %>%
	mutate(sex = factor(???)) %>%
  ggplot(aes(x=sex,y =height)) + geom_violin()
```


# Downloading a dataset
For this practical, on day 5, you will analyze a gene expression dataset of your choosing using R, based on the examples from days 2-4. In the following, we will use the Gemma.R package that enables us to easily download datasets from GEMMA, a database where datasets have been manually curated.

Let's load the required packages:
```
require(gemma.R)
require(tidyverse)
```

First, we will look for a dataset of interest. As an example, we will here look for datasets related to neuroblastoma. You can of course look for any type of topic you are interested in. 
```R
searchDatasets("neuroblastoma", limit = 100, taxon = "human") %>%
  filter(geeq.batchCorrected == TRUE) %>%
  select(ee.ShortName, ee.Name, ee.ID, ee.Accession, ee.Samples)
```

Next, pick on dataset and explore the description.
```R
gse <- "GSE21713"

getDatasetsInfo(gse) %>%
  select(ee.ShortName, ee.Name, ee.ID, ee.Description)
```

Next, have a look at the design table for this dataset.
```R
d <- getDatasetDesign(gse)
str(d)
head(d)
with(d, table(batch, genotype))
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
# row.names of the design table are the sample names. Here we check whether they are all present in the expression matrix.
stopifnot(all(row.names(d) %in% colnames(e)))

# Next, let's get only the columns corresponding to sample names, make a matrix, and add gene symbols as row.names.
dataMT <- as.matrix(e[,row.names(d)])
str(dataMT)
row.names(dataMT) <- e$GeneSymbol
str(dataMT)
```

In the examples from day 2-5, we need to voom transform data (log2CPM). In GEMMA, this has already been done.
```R
boxplot(dataMT, las=2)
```

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise 1.7:`
Now, explore another term (other than "neuroblastoma") and another dataset (other than GSE21713). Report what you find in your protocol.

For more details on the final assignment see the [instructions for day 5](day5.md).