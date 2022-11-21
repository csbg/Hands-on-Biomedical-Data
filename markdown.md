# Evaluation with R Markdown
* For the evaluation, we will consider the R Markdown HTML files that you will generate from your script.
* Save **one R script for each day** to not mix exercises and R sessions from different days. (Only exception is today, where we will generate two scripts, one to explore R Markdown (see below) and one to learn R.)
* Save the R scripts as `day1.R`, `day2.R`, `day3.R`, `day4.R`, and `day5.R`.
* After running R Markdown this will generatethe following HTML files: `day1.html`, `day2.html`, `day3.html`, `day4.html`, and `day5.html`. We will collect these files from the server (your user's home directory) and use them for the evaluation.
* In order to be able to create HTML files this way, you first have to run `.libPaths("/home/handson/R/x86_64-pc-linux-gnu-library/4.0/")` in your console.

### R Markdown
For this course, you will create reports of your R code and session using R Markdown. Create an R script, save it ("Markdown_intro.R"), and then write the following:
```R
a <- 1
b <- 3
a + b
```
* Now go to `File > Knit Document` (top right). This will generate a HTML file (this may raise an error). In this case, you can just press "Cancel" and follow the next step.
* Now go to `Files` (bottom right) and open the HTML file ("Markdown_intro.html") in the Web Browser.
* Now let's format the markdown a little bit, write the following in your Markdown file and then again generate the HTML file (using `File > Knit Document`).
```R
# This is a "normal" comment inside of code blocks
str(a)

#' this is a comment that will be shown outside of code blocks
str(b)

#'# This is a major heading (level 1)
#'## This is a heading of level 2
#'### This is a heading of level 3

#' Now let's plot something
plot(1:10)

#' Let's change the size of the plot
#+ plot1, fig.width=10
plot(1:10)

#' Let's change the size of the plot again
#+ plot2, fig.width=5, fig.height=10
plot(1:10)
```