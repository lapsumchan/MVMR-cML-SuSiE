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

We will illustrate our function via the same UK Biobank (UKB) metabolite example used in our manuscript, and the outcome of interest is Alzheimer's disease (AD). The summary statistics of the 249 UKB metabolites by Borges et al.<sup>[1]</sup> are available from the OpenGWAS database<sup>[2]</sup> with `met-d` prefix:

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
[1] 0.382700852 0.175206015 0.961084990 0.023165548 0.014269879 0.001754046
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

Similar to step 1, we need the exposure IDs (`exposure.ids.subset`), a vector containing 43 `met-d` IDs, `outcome.id` and `sample.sizes.subset` (which should be same length as `exposure.ids.subset` corresponding to their sample sizes). Nonetheless for simplicity, we use

```
sample.sizes.subset <- rep(85934, 43)
```
again, as it contains the minimum sample size. Then we can run step 2 of MVMR-cML-SuSiE using

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
[1,] 0.001298515         0.001179678 0.949818624 0.004922675 0.0007830701  0.000784514
[2,] 0.005830479         0.048772676 0.002128803 0.002302367 0.0028554353  0.002940393
[3,] 0.023255814         0.023255814 0.023255814 0.023255814 0.0232558140  0.023255814
[4,] 0.023255814         0.023255814 0.023255814 0.023255814 0.0232558140  0.023255814
[5,] 0.023255814         0.023255814 0.023255814 0.023255814 0.0232558140  0.023255814
[6,] 0.023255814         0.023255814 0.023255814 0.023255814 0.0232558140  0.023255814
```
Note that Gln (glutamine) has 0.949818624 (95.0%) PIP for the first signal cluster. We noticed when there are no clear signal (likely due to overspecification on the number of clusters) from an assumed cluster (i.e., cluster 3 to 10), SuSiE assigns a uniform probability to all exposures, in this case, 0.023255814 = 1/43. Thus, we can utilize the following code:
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
and re-run MVMR-cML for each of the 1 x 18 = 18 probable models. Using 1 of the 18 models with `met-d-Gln` and `met-d-ApoB_by_ApoA1` as an example:

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

n <- 85934

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
### References

[1] Borges, Maria Carolina, et al. "Role of circulating polyunsaturated fatty acids on cardiovascular diseases risk: analysis using Mendelian randomization and fatty acid genetic association data from over 114,000 UK Biobank participants." BMC medicine 20.1 (2022): 1-14.

[2] Elsworth, Ben, et al. "The MRC IEU OpenGWAS data infrastructure." BioRxiv (2020): 2020-08.

[3] Bellenguez, Céline, et al. "New insights into the genetic etiology of Alzheimer’s disease and related dementias." Nature genetics 54.4 (2022): 412-436.

[4] Bulik-Sullivan, Brendan, et al. "An atlas of genetic correlations across human diseases and traits." Nature genetics 47.11 (2015): 1236-1241.

[5] Wang, Gao, et al. "A simple new approach to variable selection in regression, with application to genetic fine mapping." Journal of the Royal Statistical Society Series B: Statistical Methodology 82.5 (2020): 1273-1300.
