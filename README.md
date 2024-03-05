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

We also need to prepare a vector of sample sizes corresponding to each exposures in `sample.sizes`. However, only the minimum sample size amongst GWASs (exposure + outcome) will be used for cML. Thus for illustrative purpose, it suffice to set:

```
sample.sizes <- rep(85934, 249)
```

where 85934 is the minimum sample size coming from the outcome AD (but not exposure).

With the above three `R` objects, we are ready to run step 1 of MVMR-cML-SuSiE, which provides a vector of *p*-values for all 249 metabolite exposures:

```
step1.res <- mvmr.cml.susie.step1(exposure.ids, outcome.id, sample.sizes)
```

which, upon finishing looks like this:
```
head(step1.res)
[1] 0.521785682 0.472388096 0.012445326 0.005743464 0.747003373 0.023165548
```

Notice that it may take a while to run UVMR-cML on 249 metabolites in real time, depending on the traffic of OpenGWAS. The end results of this step are provided in this Github for convenience and can be loaded using:
```
step1.res <- readRDS("step1res.RDS")
```

With this, we can extract a list of metabolite exposure which we wish to further investigate:
```
subset.idx <- which(step1.res < 0.05 / 27)
exposure.ids.subset <- exposure.ids[subset.idx]
```
which is stored in `exposure.ids.subset` and will be used in step 2 of MVMR-cML-SuSiE. In this case, we used a Bonferroni correction cutoff of `0.05 / 27`, where 27 corresponds to the number of principal components (PCs) explaining 95% variability of the metabolite. This gives us 43 metabolite exposures remaining to work with. Notice the cutoff being used is study-specific and should be varied accordingly.

Similar to step 1, we need the exposure IDs (`exposure.ids.subset`), in this case, a vector containing 43 `met-d` IDs, `outcome.id` and `sample.sizes.subset` (which should be same length as `exposure.ids.subset` corresponding to their sample sizes). Nonetheless for simplicity, we use

```
sample.sizes.subset <- rep(85934, 43)
```
again, as it contains the minimum sample size. Then we can run step 2 of MVMR-cML-SuSiE using

```
step2.res <- mvmr.cml.susie.step2(exposure.ids.subset, outcome.id, sample.sizes.subset)
```

and this step should complete within half an hour on a standard computer. Nonetheless, for convenience we also provide the end results which can be loaded using:
```
step2.res <- readRDS("step2res.RDS")
```

This step provides the OpenGWAS harmonized data needed for multivariable Mendelian randomization (MVMR) stored in `step2.res$mvdat`, of which the matrix associated with the exposure has dimension 187 x 43:
```
dim(step2.res$mvdat$exposure_beta)
[1] 187  43
```
In other words, there are 187 instrumental variables (IVs) for the 43 exposures after data harmonization. Notice this harmonized data has not been filtered for invalid IVs. The indices of invalid IVs which are stored in `step2.res$invalid.idx`:
```
step2.res$invalid.idx
[1]  11  28  62  83  89  91  99 122 131 144
```

In addition, the initial estimates for exposures used for iterative SuSiE algorithm in step 3 are stored in `step2.res$theta.vec`.

Finally, in order to run the iterative SuSiE algorithm, we just need the (43 + 1 = 44) x 44 genetic correlation matrix, which was obtained using bivariate linkage disequilibrium score (LDSC) regression<sup>[4]</sup]. The entire genetic 250 x 250 correlation matrix is provided here:
```
rho.mat <- matrix(0, 250, 250)
rho.mat[1:249,1:249] <- readRDS("metdrho.RDS")
rho.mat[250,250] <- 1
```
Notice we need to subset the indices accordingly to obtain the 44 x 44 matrix:
```
rho.mat <- rho.mat[c(subset.idx,250),c(subset.idx,250)]
```
Finally, we can obtain the SuSiE results using:
```
step3.res <- mvmr.cml.susie.step3(step2.res$mvdat, step2.res$invalid.idx, step2.res$theta.vec, rho.mat)
```
which resembles that of a SuSiE output. The most relevant output is 10 x 43 posterior inclusion probability (PIP) matrix, where 10 is the number of assumed signal clusters in SuSiE (default). Below shows six rows and columns of this 10 x 43 PIP matrix:
```
head(step3.res$alpha)
     met-d-ApoA1 met-d-ApoB_by_ApoA1   met-d-Gln met-d-GlycA  met-d-IDL_C met-d-IDL_CE
[1,] 0.001298447         0.001179607 0.949821743 0.004922351 0.0007830185  0.000784462
[2,] 0.005830548         0.048773557 0.002128769 0.002302330 0.0028553836  0.002940338
[3,] 0.023255814         0.023255814 0.023255814 0.023255814 0.0232558140  0.023255814
[4,] 0.023255814         0.023255814 0.023255814 0.023255814 0.0232558140  0.023255814
[5,] 0.023255814         0.023255814 0.023255814 0.023255814 0.0232558140  0.023255814
[6,] 0.023255814         0.023255814 0.023255814 0.023255814 0.0232558140  0.023255814
```
noticing that Gln (glutamine) has 0.949821743 (95.0%) PIP for the first signal cluster.
The key is that

# TLDR

Step 1 of MVMR-cML-SuSiE narrows down the set of promising metabolites from 249 to 43 after applying Bonferroni correcction. The harmonized data for MVMR analysis is stored in `mvdat`, while the UVMR-estimates for the 43 exposures and the list of invalid IVs (their corresponding indices) identified in step 2 are stored in `theta.vec` and `invalid.idx` of the `step2.res` list, respectively.

### References

[1] Borges, Maria Carolina, et al. "Role of circulating polyunsaturated fatty acids on cardiovascular diseases risk: analysis using Mendelian randomization and fatty acid genetic association data from over 114,000 UK Biobank participants." BMC medicine 20.1 (2022): 1-14.

[2] Elsworth, Ben, et al. "The MRC IEU OpenGWAS data infrastructure." BioRxiv (2020): 2020-08.

[3] Bellenguez, Céline, et al. "New insights into the genetic etiology of Alzheimer’s disease and related dementias." Nature genetics 54.4 (2022): 412-436.

[4] Bulik-Sullivan, Brendan, et al. "An atlas of genetic correlations across human diseases and traits." Nature genetics 47.11 (2015): 1236-1241.
