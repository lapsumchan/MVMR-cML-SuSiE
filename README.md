### Installation  

First, install `MVMR-cML-SuSiE` from GitHub using `devtools`:  

    # Install devtools if not already installed
    # install.packages("devtools") 
    devtools::install_github("lapsumchan/MVMR-cML-SuSiE")
    
Next, we need to install a few dependencies [`TwoSampleMR`](https://mrcieu.github.io/TwoSampleMR/), [`MRcML`](https://github.com/xue-hr/MRcML) and [`MVMRcML`](https://github.com/ZhaotongL/MVMR-cML):  

    install.packages("remotes")
    remotes::install_github("MRCIEU/TwoSampleMR")

    install.packages("devtools")
    devtools::install_github("xue-hr/MRcML")
    devtools::install_github("ZhaotongL/MVMR-cML")
    
which should install within a couple of minutes on a standard machine. Please refer to the above hyperlinked Github pages for help if the corresponding `R` packages cannot be installed properly.

# Demo

We first load all the source code dependencies:

```
library(TwoSampleMR)
library(MRcML)
library(MVMRcML)
```
