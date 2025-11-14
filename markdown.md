# Evaluation and R Markdown
* Ideally you can create **HTML Markdown files** from your **R scripts** (*.R) (note: You can but do not need to make Markdown files *.Rmd). 
* Save **one R script for each day** to not mix exercises and R sessions from different days. (Only exception is the code below to explore R Markdown (see below).
* Save your R scripts as `prelearning.R`, `day1.R`, `day2.R`, `day3.R`, `day4.R`, and `day5.R`.
* These R scripts can then be compiled into HTML Markdown files: `prelearning.html`, `day1.html`, `day2.html`, `day3.html`, `day4.html`, and `day5.html`. We will collect these files from the server (your user's home directory) and use them for the evaluation.
<!-- * On Ricarda, in order to be able to create HTML files this way, you first have to run `.libPaths("/home/handson/R/x86_64-pc-linux-gnu-library/4.0/")` in your console. -->
<!-- * **On RIcarda, we can re-run your R scripts and do not need HTML files. Please do not run Knitr on Ricarda as it overwhelms the server!** -->

### R Markdown
For this course, you will create reports of your R code and session using R Markdown. 

To try this out, create an R script, save it ("Markdown_intro.R"), and then write the following:
```R
a <- 1
b <- 3
a + b
```
* To compile the HTML Markdown file, now go to `File > Knit Document` **or** `File > Compile Report` (this depends on your version of R Studio and R). This will generate (and usually also open and show) a HTML file. If opening the file throws and error, you can just press "Cancel" and follow the next step.
* Go to `Files` (bottom right) and open the HTML file ("Markdown_intro.html") in the Web Browser. You now have created a Markdown HTML file.
* Now let's format the Markdown a little bit. Write the comments and code below in your R script file. This will add sections and plots including formatting. Then again generate the HTML file (using `File > Knit Document` or `File > Compile Report`).

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