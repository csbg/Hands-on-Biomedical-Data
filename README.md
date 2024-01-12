# Hands-on-Biomedical-Data
Practical exercises for the course "Hands-on Biomedical Data - Resources and Analysis Tools"

## Comments up front
Make sure you read instructions in detail. Especially getting the Setup right.

## Exercises
* Prelearning: [R 4 Data Science](day0_prelearning.md)
* Day 1: [Basic R programming and visualizations](day1.md)
* Day 2: [Introduction in differential expression analysis](day2.md)
* Day 3: [A more complex example](day3.md)
* Day 4: [Interaction effects](day4.md)
* Day 5: [Assignment - analyze a dataset of your choosing](day5.md)

## Setup
You can run the practicals on:
* [Your personal computer](setup_PersonalLaptops.md)(ideal case)
* [Within Galaxy](setup_GALAXY.md)
* [On the PLUS Server RICARDA](setup_RICARDA.md)

## Evaluation
For the evaluation, you will get points based on the exercises indicated like this:
![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Exercise X:`

* The evaluation is based on a protocol that you will prepare.
* In this protocol your should address all exercises. Each exercise counts for 1 point unless otherwise stated.
* Usually exercises are just one or two plots. If you are asked to respond to questions, max. 2-3 sentences per exercise should be sufficient.
* You can should ideally use Markdown (see instructions below) to create the protocol. This contains code and plots together and makes it very easy to track and evaluate your progress.
* Alternatively, you can copy/paste plots and write answers in Powerpoint, Word or similar (convert and submit a PDF file). In this case, you also have to submit R Scripts which document your code.

## Instructions
* Save **one R script for each day** to not mix exercises and R sessions from different days.
* Save the R scripts as `day1.R`, `day2.R`, `day3.R`, `day4.R`, and `day5.R`.
* While you can execute commands from your script in any order, make sure your finally submitted script runs through from top to bottom if started from an empty environment!
* Submit all files through Blackboard.
* The deadline is (XXX).

## Markdown
* If you work on your personal computer, you can combine code and answers to questions using Markdown. See the following: [Markdown instructions](markdown.md)
* If you do use Markdown (on your personal computer), create HTML files using `File > Knit Document` in R. Please **DO NOT** use `File > Knit Document` on Ricarda.


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
