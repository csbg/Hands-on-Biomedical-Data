# Setup on your own machine

## Getting R studio setup
* Rstudio and R are already installed
* Visit ricarda.came.sbg.ac.at
* Enter your username and password
* Click on the top left `File > R Script` to create a new R script (see image below)
* Write `print("hello world")` in the R script
* Press `Strg + ENTER`, then the code should get executed in the console


## Install packages
Packages are already installed. Type (copy/paste) the following in the beginning of each of your R scripts.
This will enable you to then load packages from the HOME directory of user "handson", where the packages required for this course are installed.
```R
.libPaths("/home/handson/R/x86_64-pc-linux-gnu-library/4.0/")
```

## Access data
To load the data needed, type (copy/paste) the following in the beginning of each of your R scripts.
```R
data <- readRDS("/home/handson/data/data.RDS")
design <- readRDS("/home/handson/data/design.RDS")
gmap <- readRDS("/home/handson/data/gmap.RDS")
```

## Start exercises
You are now ready to start the exercises.