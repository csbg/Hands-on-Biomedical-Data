# R programming and visualizations
In this part, we will get started with R and Rstudio. Throughout this course, you will run R commands from these files in Rstudio. You can type them but you can also copy/paste them. One way is to open this window on side of your screen, and open you R studio in a different window, on the other side of the screen. If you can use another browser, you may be able to use `Strg + TAB` to switch bettween the two windows using your keyboard (and don't have to click between tabs).

## Getting R studio setup
* Visit ricarda.came.sbg.ac.at
* Enter your username and password
* Click on the top left `File > R Script` to create a new R script (see image below)
* Write `print("hello world")` in the R script
* Press `Strg + ENTER`, then the code should get executed in the console

<img src="img/RStudio.png" width="100%">

## R studio basics
For this part of the practical (next 5 days), you will be writing R code in R scripts:
* When you log in you have started an `R session`
* When you execute commands in the  console (bottom left) the output and potential error or warning messages will appear in the console 
* When you also create variables, they get stored in your `environment` (see top right)
* You can restart your session by quitting the session (top right, red button)
* You can also clear objects from your environment (also called "workspace") by clicking the small yellow broom icon on the top right under "Environment"
* Plots and help packages appear in the bottom right
* We will almost always work with R scripts (top left) and use `Strg + ENTER` to "send" individual commands or sections of commands (marking them with the mouse) to the console, where they will be executed
* Save **one R script for each day**.
***While you can execute commands from your script in any order, make sure your script runs through from top to bottom if started from scratch - so do not use variables before you define them!***


## R Markdown
For this course, you will create reports of your R code and session using R Markdown. Write the following in your R script:
```
a <- 1
b <- 3
a + b
```
* Now go to `File > Knit Document` (top right). This will generate a HTML file. 
* Now go to `Files` (bottom right) and open the HTML file in the Web Browser.




## Loading packages
Type (copy/paste) the following into your R script. This will enable you to load packages from the HOME directory of user "handson", where the packages required for this course are installed. In your research, you will want (and have) to install your own packages, which in R just requires the command `install.packages(PACKAGE_NAME)`. For this course we chose to install packages centrally in order to (1) not duplicated installations and (2) make sure everyone uses the same packages
```R
.libPaths("/home/handson/R/x86_64-pc-linux-gnu-library/4.0/")
```

Once packages are installed (which they are already), you only need to load them in your current session, load the following packages:
```
require(tidyverse)
require(pheatmap)
```
