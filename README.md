# MVMR-cML-SuSiE

### Overview
This repository provides a demonstration on how to use the MVMR-cML-SuSiE `R` source code. The default option is dependent on the OpenGWAS database<sup>[1]</sup>, but users can also provide their own harmonized Mendelian randomization (MR) data (see TLDR towards the bottom of this README).

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
    
which should install within a couple of minutes on a standard machine. Please refer to the above hyperlinked Github pages for help if the corresponding `R` packages cannot be installed properly.

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

We will illustrate our function via the same UK Biobank (UKB) metabolite example used in our manuscript, and the outcome of interest is Alzheimer's disease (AD). The summary statistics of the 249 UKB metabolites by Borges et al.<sup>[2]</sup> are available from the OpenGWAS database<sup>[1]</sup> with `met-d` prefix:

```
ao <- available_outcomes()

# Use grep to find ids that start with "met-d"
metd.idx <- grep("^met-d", ao$id)
exposure.ids <- ao$id[metd.idx]

# Make the ids in alphabetical order (to match the ordering of summary statistics given by TwoSampleMR package, i.e., mvdat in step 2 below)
exposure.ids <- sort(exposure.ids)
```

As for AD summary statistics, we will be using the largest AD cohort by Bellenguez et al.<sup>[3]</sup>, which is also available in [OpenGWAS](https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST90027158/). The corresponding ID is given by:
```
outcome.id <- "ebi-a-GCST90027158"
```

Notice that we need the minimum sample size amongst GWASs (exposure + outcome) for cML. Thus, we need to prepare a vector of sample sizes corresponding to each exposures in `sample.sizes`. This has been prepared in the file `metdn.RDS` so we just need to load it:
```
sample.sizes <- readRDS("metdn.RDS")
```

With the above three `R` objects, we are ready to run step 1 of MVMR-cML-SuSiE, which provides a vector of *p*-values for all 249 metabolite exposures:

```
step1.res <- mvmr.cml.susie.step1(exposure.ids, outcome.id, sample.sizes)
```

which, upon finishing looks like this:
```
head(step1.res)
[1] 0.401999338 0.172905493 0.963813629 0.026442776 0.013637944 0.001762121
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

Similar to step 1, we need the exposure IDs (`exposure.ids.subset`), a vector containing 43 `met-d` IDs, `outcome.id` and `sample.sizes.subset` (which should be same length as `exposure.ids.subset` corresponding to their sample sizes). We can obtain the sample sizes corresponding to the 43 exposures using:

```
sample.sizes.subset <- sample.sizes[subset.idx]
```

Now, we can run step 2:
```
step2.res <- mvmr.cml.susie.step2(exposure.ids.subset, outcome.id, sample.sizes.subset)
```

and this step should complete within half an hour on a standard computer. Again for convenience we provide the end results which can be loaded using:
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

In addition, the initial causal estimates for exposures used for the iterative SuSiE algorithm in step 3 are stored in `step2.res$theta.vec`.

Finally, in order to run the iterative SuSiE algorithm, we just need the (43 + 1 = 44) x 44 genetic correlation matrix, which was obtained using bivariate linkage disequilibrium score (LDSC) regression<sup>[4]</sup>. The entire genetic 250 x 250 correlation matrix is provided here:
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
which resembles that of a SuSiE output. The most relevant output is the 10 x 43 posterior inclusion probability (PIP) matrix, where 10 is the number of assumed signal clusters in SuSiE (default). Below shows six rows and columns of this 10 x 43 PIP matrix:
```
head(step3.res$alpha)
     met-d-ApoA1 met-d-ApoB_by_ApoA1   met-d-Gln met-d-GlycA  met-d-IDL_C met-d-IDL_CE

[1,] 0.001298515         0.001179678 0.949818626 0.004922675 0.0007830701 0.0007845139
[2,] 0.005830479         0.048772676 0.002128803 0.002302367 0.0028554352 0.0029403935
[3,] 0.023255814         0.023255814 0.023255814 0.023255814 0.0232558140 0.0232558140
[4,] 0.023255814         0.023255814 0.023255814 0.023255814 0.0232558140 0.0232558140
[5,] 0.023255814         0.023255814 0.023255814 0.023255814 0.0232558140 0.0232558140
[6,] 0.023255814         0.023255814 0.023255814 0.023255814 0.0232558140 0.0232558140
```
Note that Gln (glutamine) has 0.949818626 (95.0%) PIP for the first signal cluster. We noticed when there are no clear signal (likely due to overspecification on the number of clusters) from an assumed cluster (i.e., cluster 3 to 10), SuSiE assigns a uniform probability to all exposures, in this case, 0.0232558140 = 1/43. Thus, we can utilize the following code:
```
sort(unique(which(step3.res$alpha > 1/43, arr.ind = TRUE)[,1]))
[1] 1 2
```
to find out there are in total 2 signal clusters. Fortunately, SuSiE is robust to the overspecification of the number of signal clusters<sup>[5]</sup>. If we re-run step 3 by assuming 2 clusters, we pretty much obtain the same results:
```
step3.res.v2 <- mvmr.cml.susie.step3(step2.res$mvdat, step2.res$invalid.idx, step2.res$theta.vec, rho.mat, 2)
head(step3.res.v2$alpha)
     met-d-ApoA1 met-d-ApoB_by_ApoA1   met-d-Gln met-d-GlycA  met-d-IDL_C met-d-IDL_CE
[1,] 0.001298390         0.001179498 0.949823273 0.004922889 0.0007830073 0.0007844499
[2,] 0.005830347         0.048772690 0.002128719 0.002302289 0.0028553540 0.0029403107
```

With this, we can identify the exposures associated with each signal cluster:
```
which(step3.res$alpha[1,] > 1/43)
met-d-Gln 
        3 

which(step3.res$alpha[2,] > 1/43)
met-d-ApoB_by_ApoA1       met-d-L_LDL_P     met-d-M_VLDL_FC      met-d-M_VLDL_P 
                  2                   9                  18                  19 
    met-d-M_VLDL_PL       met-d-S_LDL_P      met-d-S_VLDL_C     met-d-S_VLDL_CE 
                 20                  22                  23                  24 
    met-d-S_VLDL_FC      met-d-S_VLDL_L      met-d-S_VLDL_P     met-d-S_VLDL_PL 
                 25                  26                  27                  28 
       met-d-VLDL_C       met-d-VLDL_CE        met-d-VLDL_P     met-d-XS_VLDL_L 
                 29                  30                  31                  36 
    met-d-XS_VLDL_P    met-d-XS_VLDL_PL 
                 37                  38 
```
and re-run MVMR-cML<sup>[6]</sup> for each of the 1 x 18 = 18 probable models. Using 1 of the 18 models with `met-d-Gln` and `met-d-ApoB_by_ApoA1` as an example:

```
exposure.dat.sub1 <- mv_extract_exposures(c("met-d-ApoB_by_ApoA1", "met-d-Gln"))
outcome.dat.sub1 <- extract_outcome_data(snps=exposure.dat.sub1$SNP, outcomes = outcome.id)
mvdat.sub1 <- mv_harmonise_data(exposure.dat.sub1, outcome.dat.sub1)

K.sub1 <- dim(mvdat.sub1$exposure_beta)[1]

ids.sub1 <- c(2, 3) # index 2 corresponds to met-d-ApoB_by_ApoA1 and 3 corresponds to met-d-Gln
rho.mat.sub1 <- rho.mat[c(ids.sub1,44), c(ids.sub1,44)]

Sig.inv.l.sub1 <- invcov_mvmr(se_bx = mvdat.sub1$exposure_se,
                         se_by = mvdat.sub1$outcome_se,
                         rho_mat = rho.mat.sub1)

n <- min(sample.sizes.subset[ids.sub1])

MVcML.res.sub1 <- MVmr_cML_DP(b_exp = mvdat.sub1$exposure_beta,
                         b_out = as.matrix(mvdat.sub1$outcome_beta),
                         se_bx = mvdat.sub1$exposure_se,
                         Sig_inv_l = Sig.inv.l.sub1, n = n, num_pert = 100,
                         K_vec = 0:(K.sub1-3))

MVcML.BIC.SE.sub1 <- MVcML_SdTheta(b_exp = mvdat.sub1$exposure_beta,
                              b_out = as.matrix(mvdat.sub1$outcome_beta),
                              Sig_inv_l = Sig.inv.l.sub1,
                              theta = MVcML.res.sub1$BIC_theta,
                              zero_ind = setdiff(1:length(mvdat.sub1$outcome_beta), MVcML.res.sub1$BIC_invalid))

MVcMLBIC.pval.sub1 <- pnorm(-abs(MVcML.res.sub1$BIC_theta/MVcML.BIC.SE.sub1))*2
```
we obtain the causal effect estimates for glutamine and ratio of apolipoprotein B to apolipoprotein A1 (ApoB/A1) to be both -0.12:
```
MVcML.res.sub1$BIC_theta
           [,1]
[1,] -0.1232354
[2,] -0.1182489
```
and this indicates both glutamine and ApoB/A1 has a protective effect for AD. Moreover, this result is highly significant:
```
MVcMLBIC.pval.sub1
             [,1]
[1,] 9.175066e-07
[2,] 2.487480e-05
```
For more detailed usage of the MVMR-cML `R` package, please refer to the Github page of [`MVMR-cML`](https://github.com/ZhaotongL/MVMR-cML).

# TLDR
Users can provide their own version of harmonized data. For step 1, we require length `L` lists of summary statistics coefficients (beta) and standard errors (se) for both the exposures and outcomes. This is basically providing the univariable MR (UVMR) harmonized data for each exposure (and the outcome summary statistics corresponding to the IVs used). Notice that the set of `m` IVs should be independent (can be achieved by LD clumping), as this is a requirement for both UVMR-cML<sup>[7]</sup>, MVMR-cML<sup>[6]</sup> as well as our method (which builds upon on the former two methods). In our case, `L = 249`. In addition, we also require a vector of length `L + 1` containing the sample sizes for each of the `L` exposures and the outcome GWAS (last element). The `metdn.RDS` file contains sample sizes for the 249 UKB exposures, while 487511 is the sample size for the AD (outcome) GWAS. Below shows all 5 objects required for step 1 if the users were to provide their own data:

```
sample.sizes <- readRDS("metdn.RDS")
sample.sizes <- c(sample.sizes, 487511)

beta.exposure.ls <- readRDS("beta.exposure.ls.RDS")
se.exposure.ls <- readRDS("se.exposure.ls.RDS")
beta.outcome.ls <- readRDS("beta.outcome.ls.RDS")
se.outcome.ls <- readRDS("se.outcome.ls.RDS")
```

which upon running
```
step1.res <- mvmr.cml.susie.step1(sample.sizes = sample.sizes, beta.exposure.ls = beta.exposure.ls, se.exposure.ls = se.exposure.ls, beta.outcome.ls = beta.outcome.ls, se.outcome.ls = se.outcome.ls, use.openGWAS = FALSE)
```
with the `use.openGWAS` option as `FALSE` should yield identical results as the OpenGWAS dependent version in the README.

Based on step 1, it should suggest a subset of `L.star` exposures that are further analysis-worthy:
```
subset.idx <- which(step1.res < 0.05 / 27)
sample.sizes.subset <- sample.sizes[subset.idx]
sample.sizes.subset <- c(sample.sizes.subset, 487511)
```
and in this case, `L.star = 43`. The users will then need to provide a joint set of `m.star` IVs for the `L.star` exposures: two matrices for the exposure beta and se (both `m.star x L.star`), and two vectors for outcome beta and se (both length `m.star`), as well as a `m.star x L.star` *p*-value matrix from the exposure GWAS (only IVs reaching `cutoff` argument, which default is set to `5e-8` are used for the re-analysis). Again, it cannot be overstated enough that the set of `m.star` IVs should be independent. Below loads the 5 objects required for step 2 if the users were to provide their own data (in addition to `sample.sizes.subset`):

```
beta.exposure.mat <- readRDS("beta.exposure.mat.RDS")
se.exposure.mat <- readRDS("se.exposure.mat.RDS")
beta.outcome.vec <- readRDS("beta.outcome.vec.RDS")
se.outcome.vec <- readRDS("se.outcome.vec.RDS")
pval.exposure.mat <- readRDS("pval.exposure.mat.RDS")
```
which upon running
```
step2.res <- mvmr.cml.susie.step2(sample.sizes.subset = sample.sizes.subset, beta.exposure.mat = beta.exposure.mat, se.exposure.mat = se.exposure.mat, beta.outcome.vec = beta.outcome.vec, se.outcome.vec = se.outcome.vec, pval.exposure.mat = pval.exposure.mat, use.openGWAS = FALSE)
```
should also yield identical results as the OpenGWAS dependent version in the README.

Finally, step 3 only depends on `mvdat`, `invalid.idx` and `theta.vec` (all been put together in step 2) and the genetic correlation matrix `rho.mat` (which should be `(L.star + 1) x (L.star + 1)`):
```
rho.mat <- matrix(0, 250, 250)
rho.mat[1:249,1:249] <- readRDS("metdrho.RDS")
rho.mat[250,250] <- 1

rho.mat <- rho.mat[c(subset.idx,250),c(subset.idx,250)]

step3.res <- mvmr.cml.susie.step3(step2.res$mvdat, step2.res$invalid.idx, step2.res$theta.vec, rho.mat)
```
Please refer back to the above README on `mvmr.cml.susie.step3` for the intepretation of the results (this is getting long!).

### References

[1] Elsworth, Ben, et al. "The MRC IEU OpenGWAS data infrastructure." BioRxiv (2020): 2020-08.

[2] Borges, Maria Carolina, et al. "Role of circulating polyunsaturated fatty acids on cardiovascular diseases risk: analysis using Mendelian randomization and fatty acid genetic association data from over 114,000 UK Biobank participants." BMC medicine 20.1 (2022): 1-14.

[3] Bellenguez, Céline, et al. "New insights into the genetic etiology of Alzheimer’s disease and related dementias." Nature genetics 54.4 (2022): 412-436.

[4] Bulik-Sullivan, Brendan, et al. "An atlas of genetic correlations across human diseases and traits." Nature genetics 47.11 (2015): 1236-1241.

[5] Wang, Gao, et al. "A simple new approach to variable selection in regression, with application to genetic fine mapping." Journal of the Royal Statistical Society Series B: Statistical Methodology 82.5 (2020): 1273-1300.

[6] Lin, Zhaotong, Haoran Xue, and Wei Pan. "Robust multivariable Mendelian randomization based on constrained maximum likelihood." The American Journal of Human Genetics 110.4 (2023): 592-605.

[7] Xue, Haoran, Xiaotong Shen, and Wei Pan. "Constrained maximum likelihood-based Mendelian randomization robust to both correlated and uncorrelated pleiotropic effects." The American Journal of Human Genetics 108.7 (2021): 1251-1269.
