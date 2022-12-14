# Hands-on-Biomedical-Data
Practical exercises for the course "Hands-on Biomedical Data - Resources and Analysis Tools"

## Content
* Day 1: [R programming and visualizations](day1.md)
* Day 2: [Introduction in differential expression analysis](day2.md)
* Day 3: [A more complex example](day3.md)
* Day 4: [Interaction effects](day4.md)
* Day 5: [Assignment - analyze a dataset of your choosing](day5.md) (The exercises for day 5 are not be part of the course in 2022. Only the first 4 days need to be completed.)

Ideally you can run the practicals on your own machine, see here for [installation instructions](setup.md). If you do so, please see the [Markdown instructions](markdown.md).

## Evaluation
* For the evaluation, we will use R Scripts files that you will generate.
* Save **one R script for each day** to not mix exercises and R sessions from different days.
* Save the R scripts as `day1.R`, `day2.R`, `day3.R`, `day4.R`, and `day5.R`.
* If you work on your own machine, create HTML files using `File > Knit Document` in R. Please **DO NOT** use `File > Knit Document` on Ricarda.
* While you can execute commands from your script in any order, make sure your finally submitted script runs through from top to bottom if started from an empty environment!
* Use comments `#` to answer questions and comment on your code and results.
* We will evaluate all scripts a week after Day 5 has finished (Dec 2nd, 2022).

## The example dataset

In this part of the practical, we will study transcriptomics data of structural cells in mice upon cytokine stimulation from [Krausgruber & Fortelny et al., Nature, 2020](https://doi.org/10.1038/s41586-020-2424-4).
<br/>
<br/>
<img src="img/StructuralImmunity.png" width="50%">

## Basic analyses steps covered
* Quality control using sample correlations and dimensionality reduction
* Data normalization
* Differential expression
* Model diagnostics and quality control
* Plotting of results
* Interpretation of top genes
* Gene set enrichment analysis

## How to get help?
* Most commands should be explained in this practical.
* If you do not understand certain functions, type the question mark plus the function name in R, e.g.: "?median".
* If you need additional commands, Google is your friend. 
* Also consult this [list of function names](functions.md), which contains key functions relevant for this course.
* **Don't forget to raise you hand if lost!**