# Setup on GALAXY
## Galaxy for prelearning exercises
* If you want to use galaxy for pre-learning exercise, you could go to this link
https://usegalaxy.eu/login/start?redirect=None

Create an account/login

Follow the instructions from Step 3 onwards from below (Please note that the getting the shared data mentioned below is not important for the pre-learning exercise, but only during the course)

## Getting R studio setup
* Rstudio and R are already installed
* Visit
 
https://usegalaxy.eu/join-training/hands-on-plus

This URL is active during your training (from 2024-11-18 to 2024-12-20). Once it is over, the link will not be usable anymore but the users can still access their data at usegalaxy.eu.)


**Step 1**. Sign-up and activate your account (email) if you haven't done already

**Step 2**. Login to your account
 Click on Shared data --> Histories --> type in  **Hands_on_data** (owner-aarathy)-->![image](https://github.com/csbg/Hands-on-Biomedical-Data/assets/96147982/32dac7c2-904d-494e-97ea-454aee0f80bf)
-->click on it-->import this history (on the top left corner)

**Step 3**. Click on the Galaxy icon to get back to home. You will now have the data in your history. Notice that each data set in your history has a number.

 ![image](https://github.com/user-attachments/assets/4f676bd4-cbd1-4bd2-8ae2-acb9ad812e4d)
 
**Step 4**. To open Rstudio, on the top-left corner, search for Rstudio
--> run tool -->open or click on the eye icon
![image](https://github.com/user-attachments/assets/b2ea7c1a-6f14-459e-a151-9301e5c20903)


* Write `print("hello world")` in the R script
* Press `Strg + ENTER`, then the code should get executed in the console


## Install packages
* Many packages are already installed. 
* If not installed, to install packages try
```R
install.packages("package name")
```
* If it doesn't work with install.packages command, change to the Tab "Terminal", next to Console and Background Jobs.

![image](https://github.com/user-attachments/assets/11072545-faa7-44fa-a1fd-a95a42f455ba)




* Now you type in "conda install r-tidyverse" or "conda install bioconductor-limma" etc to install the respective packages.
```R
#type in the Rscript and press `Strg + ENTER`
require("package name") #to load the packages that are already installed
```
## Access data
To load the data needed, type (copy/paste) the following in the beginning of each of your R scripts.

name_to_be_stored <- readRDS(get_gx(**number associated to the respective dataset in your history**))
![image](https://github.com/user-attachments/assets/edb3a8b8-f614-40f2-9940-b68c48ee9e93)


```R
conda install r-tidyverse
conda install r-pheatmap
conda install r-enrichR
conda install r-devtools
conda install r-knitr
conda install r-patchwork

conda install bioconductor-limma
conda install bioconductor-ComplexHeatmap
conda install bioconductor-gemma.R
```

```R
# In this case
data <- readRDS(gx_get(3))
design <- readRDS(gx_get(2))
gmap <- readRDS(gx_get(1))
```

## Start exercises
You are now ready to start the exercises.
