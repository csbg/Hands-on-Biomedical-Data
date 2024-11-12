# Setup on GALAXY
## Galaxy for prelearning exercises
* If you want to use galaxy for pre-learning exercise, you could go to this link
https://usegalaxy.eu/login/start?redirect=None

Create an account/login

Follow the instructions from Step 3 onwards from below (Please note that the getting the shared data mentioned below is not important for the pre-learning exercise, but only during the course)

## Getting R studio setup
* Rstudio and R are already installed
* Visit 
https://eur05.safelinks.protection.outlook.com/?url=https%3A%2F%2Fusegalaxy.eu%2Fjoin-training%2Fplus-course-r&data=05%7C01%7Caarathy.ravisundarjosegeetha%40plus.ac.at%7Ceeb033492e0c4826baa108dbe6b22b1b%7C158a941a576e4e87993db2eab8526e50%7C1%7C0%7C638357425131951475%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C3000%7C%7C%7C&sdata=U4uhYkF0iVDtbGP3jpv3ZUm1Vev%2FmfBZ57YPDF6da7s%3D&reserved=0

Will be updated soon!!
This URL is active during your training (from 2023-11-20 to 2023-12-12). Once it is over, the link will not be usable
anymore but the users can still access their data at usegalaxy.eu.

**Step 1**. Sign-up and activate your account (email) if you haven't done already

**Step 2**. Login to your account
 Click on Shared data --> Histories --> type in  **Hands_on_data** (owner-aarathy)-->![image](https://github.com/csbg/Hands-on-Biomedical-Data/assets/96147982/32dac7c2-904d-494e-97ea-454aee0f80bf)
-->click on it-->import this history (on the top left corner)

**Step 3**. Click on the Galaxy icon to get back to home. You will now have the data in your history. Notice that each data set in your history has a number.

  ![image](https://github.com/csbg/Hands-on-Biomedical-Data/assets/96147982/ce98f67a-7a33-4401-9103-6d6dcb4d0f63)


**Step 4**. To open Rstudio, on the top-left corner, search for Rstudio
--> run tool -->open or click on the eye icon
  ![image](https://github.com/csbg/Hands-on-Biomedical-Data/assets/96147982/e5fc5127-2597-4522-9d57-77f2dadca9dc)

* Write `print("hello world")` in the R script
* Press `Strg + ENTER`, then the code should get executed in the console


## Install packages
* Many packages are already installed. 
* If not installed, to install packages try
```R
install.packages("package name")
```
* If it doesn't work with install.packages command, change to the Tab "Terminal", next to Console and Background Jobs.
* Now you type in "conda install r-tidyverse" or "conda install bioconductor-limma" etc to install the respective packages.
```R
#type in the Rscript and press `Strg + ENTER`
library("package name") #to load the packages that are already installed
```
## Access data
To load the data needed, type (copy/paste) the following in the beginning of each of your R scripts.

name_to_be_stored <- readRDS(get_gx(**number associated to the respective dataset in your history**))
```R
# In this case
data<-readRDS(gx_get(3))
design <- readRDS(gx_get(2))
gmap <- readRDS(gx_get(1))
```

## Start exercises
You are now ready to start the exercises.
