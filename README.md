# MVMR-cML-SuSiE

### Overview
This repository provides a demonstration on how to use the MVMR-cML-SuSiE `R` source code.

# System Requirements

## Software Requirements

### OS Requirements

The source code has been tested on Microsoft's Windows 10 operating system and Linux (Ubuntu 18.04). The source code should be compatible with Windows, Mac, and Linux operating systems.

Before using the MVMR-cML-SuSiE source code, users should have `R` version 4.3.0 or higher, and several packages installed.

### Installation  

First, we need to install a few dependencies [`TwoSampleMR`](https://mrcieu.github.io/TwoSampleMR/), [`susieR`](https://github.com/stephenslab/susieR), [`MRcML`](https://github.com/xue-hr/MRcML) and [`MVMRcML`](https://github.com/ZhaotongL/MVMR-cML):  

    install.packages("remotes")
    remotes::install_github("MRCIEU/TwoSampleMR")

    install.packages("susieR")

    install.packages("devtools")
    devtools::install_github("xue-hr/MRcML")
    devtools::install_github("ZhaotongL/MVMR-cML")
    
which should install within a couple of minutes on a standard machine. Please refer to above the hyperlinked Github pages for help if the corresponding `R` packages cannot be installed properly.

# Demo

We first load all the source code dependencies:

```
library(TwoSampleMR)
library(susieR)
library(MRcML)
library(MVMRcML)
```

and the source code containing all the main functions:

```
source("main.R")
```

We will illustrate our function via the same UK Biobank (UKB) metabolite example used in our manuscript, and the outcome of interest is Alzheimer's disease (AD). The summary statistics of the 249 UKB metabolites by Borges et al.<sup>[1]</sup> are available from the OpenGWAS database<sup>[2]</sup> with `met-d` prefix:

```
ao <- available_outcomes()

# Use grep to find ids that start with "met-d"
metd.idx <- grep("^met-d", ao$id)
exposure.ids <- ao$id[metd.idx]
```

As for the AD summary statistics, we will be using the largest AD cohort by Bellenguez et al.<sup>[3]</sup>, which is also available in [OpenGWAS](https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST90027158/). The corresponding ID is given by:
```
outcome.id <- "ebi-a-GCST90027158"
```

We also need to prepare a vector of sample sizes corresponding to each exposures in `sample.sizes`. However, only the minimum sample size will be used for cML. Thus for illustrative purpose, it suffice to set:

```
sample.sizes <- 85934
```

which is the minimum sample size coming from the outcome AD (but not exposure).

With the above three `R` objects, we are ready to run step 1 of MVMR-cML-SuSiE, which provides a vector of *p*-values for all 249 metabolite exposures:

```
step1.res <- mvmr.cml.susie.step1(exposure.ids, outcome.id, sample.sizes)
```

which, upon finishing looks like this:
```
head(step1.res)
[1] 0.521785682 0.472388096 0.012445326 0.005743464 0.747003373 0.023165548
```

Notice that it may take a while to run UVMR-cML on 249 metabolites in real time. The end results of this step are provided in this Github for convenience and can be loaded using:
```
step1.res <- readRDS("step1res.RDS")
```

With this, we can extract a list of metabolite exposure which we wish to further investigate:
```
subset.idx <- which(step1.res < 0.05 / 27)
exposure.ids.subset <- exposure.ids[subset.idx]
```
which is stored in `exposure.ids.subset` and will be used in step 2 of MVMR-cML-SuSiE. In this case, we used a Bonferroni correction cutoff of `0.05 / 27`, where 27 corresponds to the number of principal components (PCs) explaining 95% variability of the metabolite. This gives us 43 metabolite exposures remaining to work with. Notice the cutoff being used is study-specific and should be varied accordingly.
### References

[1] Borges, Maria Carolina, et al. "Role of circulating polyunsaturated fatty acids on cardiovascular diseases risk: analysis using Mendelian randomization and fatty acid genetic association data from over 114,000 UK Biobank participants." BMC medicine 20.1 (2022): 1-14.

[2] Elsworth, Ben, et al. "The MRC IEU OpenGWAS data infrastructure." BioRxiv (2020): 2020-08.

[3] Bellenguez, Céline, et al. "New insights into the genetic etiology of Alzheimer’s disease and related dementias." Nature genetics 54.4 (2022): 412-436.
