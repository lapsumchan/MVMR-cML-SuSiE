mvmr.cml.susie.step1 <- function(exposure.ids = NULL, 
                                 outcome.id = NULL, 
                                 sample.sizes, 
                                 beta.exposure.ls = NULL, 
                                 se.exposure.ls = NULL, 
                                 beta.outcome.ls = NULL, 
                                 se.outcome.ls = NULL, 
                                 use.openGWAS = TRUE) {
  
  # Check if OpenGWAS data extraction is required
  if (use.openGWAS) {
    
    # Ensure exposure.ids and outcome.id are provided
    if (is.null(exposure.ids) || is.null(outcome.id)) {
      stop("Please provide 'exposure.ids' and 'outcome.id' when GWAS data is not provided.")
    }
    
    L <- length(exposure.ids)
    
    # Check length of sample.sizes
    if (length(sample.sizes) != L) {
      stop("Length of 'sample.sizes' must be equal to the number of exposures (L) when GWAS data is not provided.")
    }
    
    pval.vec <- rep(NA, L)
    
    for (i in 1:L) {
      # Get instruments
      exposure.dat <- extract_instruments(exposure.ids[i])
      
      # Get effects of instruments on outcome
      outcome.dat <- extract_outcome_data(snps = exposure.dat$SNP, outcomes = outcome.id)
      
      # Harmonise the exposure and outcome data
      dat <- harmonise_data(exposure.dat, outcome.dat)
      
      n <- min(sample.sizes[i], outcome.dat$samplesize.outcome)
      
      # Perform UVMR-cML
      cML.result <- mr_cML(dat$beta.exposure,
                           dat$beta.outcome,
                           dat$se.exposure,
                           dat$se.outcome,
                           n = n,
                           random_start = 100,
                           random_seed = 1)
      pval.vec[i] <- cML.result$MA_BIC_p
    }
    
  } else {
    # Ensure user-provided data is available
    if (is.null(beta.exposure.ls) || is.null(se.exposure.ls) || is.null(beta.outcome.ls) || is.null(se.outcome.ls)) {
      stop("Please provide 'beta.exposure.ls', 'se.exposure.ls', 'beta.outcome.ls', 'se.outcome.ls' when GWAS data is provided.")
    }
    
    # Check that all inputs are lists and of the same length
    if (!is.list(beta.exposure.ls) || !is.list(se.exposure.ls) || !is.list(beta.outcome.ls) || !is.list(se.outcome.ls)) {
      stop("All beta and se inputs must be lists when GWAS data is provided.")
    }
    
    L <- length(beta.exposure.ls)
    
    if (length(se.exposure.ls) != L || length(beta.outcome.ls) != L || length(se.outcome.ls) != L) {
      stop("All lists (beta.exposure.ls, se.exposure.ls, beta.outcome.ls, se.outcome.ls) must have the same length.")
    }
    
    # Check length of sample.sizes
    if (length(sample.sizes) != (L + 1)) {
      stop("Length of 'sample.sizes' must be equal to the number of exposures (L) + 1 (outcome) when GWAS data is provided.")
    }
    
    pval.vec <- rep(NA, L)
    
    for (i in 1:L) {
      print(i)
      # Extract vectors for the ith exposure
      beta.exposure <- beta.exposure.ls[[i]]
      se.exposure <- se.exposure.ls[[i]]
      beta.outcome <- beta.outcome.ls[[i]]
      se.outcome <- se.outcome.ls[[i]]
      
      # Calculate sample size for cML
      n <- min(sample.sizes[i], sample.sizes[L + 1])
      
      # Perform UVMR-cML using provided data
      cML.result <- mr_cML(beta.exposure,
                           beta.outcome,
                           se.exposure,
                           se.outcome,
                           n = n,
                           random_start = 100,
                           random_seed = 1)
      pval.vec[i] <- cML.result$MA_BIC_p
    }
  }
  
  return(pval.vec)
}

mvmr.cml.susie.step2 <- function(exposure.ids.subset = NULL, 
                                 outcome.id = NULL, 
                                 sample.sizes.subset, 
                                 beta.exposure.mat = NULL, 
                                 se.exposure.mat = NULL, 
                                 beta.outcome.vec = NULL, 
                                 se.outcome.vec = NULL, 
                                 pval.exposure.mat = NULL,
                                 use.openGWAS = TRUE,
                                 cutoff = 5e-8) {
  
  # Check if OpenGWAS data extraction is required
  if (use.openGWAS) {
    
    # Ensure exposure.ids.subset and outcome.id are provided
    if (is.null(exposure.ids.subset) || is.null(outcome.id)) {
      stop("Please provide 'exposure.ids.subset' and 'outcome.id' when GWAS data is not provided.")
    }
    
    # Number of exposures in the subset
    L.star <- length(exposure.ids.subset)
    
    # Check length of sample.sizes.subset
    if (length(sample.sizes.subset) != L.star) {
      stop("Length of 'sample.sizes.subset' must be equal to the number of exposures (L.star) when GWAS data is not provided.")
    }
    
    # Get instruments jointly for the subset of exposures
    exposure.dat <- mv_extract_exposures(exposure.ids.subset)
    
    # Get effects of instruments on outcome
    outcome.dat <- extract_outcome_data(snps = exposure.dat$SNP, outcomes = outcome.id)
    
    # Harmonise the exposure and outcome data
    mvdat <- mv_harmonise_data(exposure.dat, outcome.dat)
    
  } else {
    # Ensure user-provided data is available
    if (is.null(beta.exposure.mat) || is.null(se.exposure.mat) || 
        is.null(beta.outcome.vec) || is.null(se.outcome.vec) || 
        is.null(pval.exposure.mat)) {  # Check for pval.exposure.mat
      stop("Please provide 'beta.exposure.mat', 'se.exposure.mat', 'beta.outcome.vec', 'se.outcome.vec', and 'pval.exposure.mat' when GWAS data is provided.")
    }
    
    # Check dimensions of input matrices and vectors
    L.star <- ncol(beta.exposure.mat)
    m <- nrow(beta.exposure.mat)
    
    if (ncol(se.exposure.mat) != L.star || nrow(se.exposure.mat) != m) {
      stop("Dimensions of 'se.exposure.mat' must match 'beta.exposure.mat'.")
    }
    
    if (length(beta.outcome.vec) != m || length(se.outcome.vec) != m) {
      stop("Length of 'beta.outcome.vec' and 'se.outcome.vec' must match the number of rows in 'beta.exposure.mat'.")
    }
    
    if (ncol(pval.exposure.mat) != L.star || nrow(pval.exposure.mat) != m) {
      stop("Dimensions of 'pval.exposure.mat' must match 'beta.exposure.mat'.")
    }
    
    # Check length of sample.sizes.subset
    if (length(sample.sizes.subset) != (L.star + 1)) {
      stop("Length of 'sample.sizes.subset' must be equal to the number of exposures (L.star) + 1 (outcome) when GWAS data is provided.")
    }
    
    mvdat <- list(exposure_beta = beta.exposure.mat,
                  exposure_se = se.exposure.mat,
                  outcome_beta = beta.outcome.vec,
                  outcome_se = se.outcome.vec,
                  exposure_pval = pval.exposure.mat)
  }
  
  uvmr.ls <- list()
  
  # Focus only on QTL for each exposure based on p-values
  for (i in 1:L.star) {
    uvmr.ls[[i]] <- which(mvdat$exposure_pval[, i] < cutoff)
  }
  
  # Initialize result vectors and lists
  pval.vec <- rep(NA, L.star)
  theta.vec <- rep(NA, L.star)
  uvmr.invalid.ls <- list()
  
  for (i in 1:L.star) {
    print(i)
    n <- min(sample.sizes.subset[i], 
             if (use.openGWAS) outcome.dat$samplesize.outcome else sample.sizes.subset[L.star + 1])
    
    # Perform UVMR-cML
    cML.result <- mr_cML(mvdat$exposure_beta[uvmr.ls[[i]], i],
                         mvdat$outcome_beta[uvmr.ls[[i]]],
                         mvdat$exposure_se[uvmr.ls[[i]], i],
                         mvdat$outcome_se[uvmr.ls[[i]]],
                         n = n,
                         random_start = 100,
                         random_seed = 1)
    
    pval.vec[i] <- cML.result$MA_BIC_p
    theta.vec[i] <- cML.result$MA_BIC_theta
    
    if (length(cML.result$BIC_invalid) != 0) {
      uvmr.invalid.ls[[i]] <- uvmr.ls[[i]][cML.result$BIC_invalid]
    } else {
      uvmr.invalid.ls[[i]] <- cML.result$BIC_invalid
    }
  }
  
  # IV is considered as invalid if invalid in any one of the exposures
  rm.idx <- sort(Reduce(union, uvmr.invalid.ls))
  return(list(mvdat = mvdat, invalid.idx = rm.idx, theta.vec = theta.vec))
}

# Iterative SuSiE
mvmr.cml.susie.step3 <- function(mvdat, invalid.idx, theta.vec, rho.mat, 
                                 L = 10, max.iter = 200, tol = 1e-10) {
  
  m.star <- dim(mvdat$exposure_beta)[1]
  L.star <- dim(mvdat$exposure_beta)[2]
  
  # Check if the dimension of rho.mat is (L + 1) x (L + 1)
  if (ncol(rho.mat) != (L.star + 1) || nrow(rho.mat) != (L.star + 1)) {
    stop(paste("The dimension of 'rho.mat' should be", (L.star + 1), "x", (L.star + 1), 
               "but it is", nrow(rho.mat), "x", ncol(rho.mat), "."))
  }
  
  X <- mvdat$exposure_beta
  Y <- mvdat$outcome_beta
  SigmaX <- mvdat$exposure_se
  seY <- mvdat$outcome_se
  
  # Create a list of genetic correlation corresponding to individual i
  sei.list <- list()
  
  for (i in 1:m.star) {
    sei.list[[i]] <- c(SigmaX[i,], seY[i])
  }
  
  # Create a list of variance corresponding to individual i
  Sigmai.list <- list()
  
  for (i in 1:m.star) {
    Sigmai.list[[i]] <- diag(sei.list[[i]]) %*% rho.mat %*% diag(sei.list[[i]])
  }
  
  # Create a list of profile likelihood variance (denominator) corresponding to individual i
  new.sigma.vec <- rep(NA, m.star)
  
  gamma.vec <- c(theta.vec, -1)
  for (i in 1:m.star) {
    new.sigma.vec[i] <- sqrt(t(gamma.vec) %*% Sigmai.list[[i]] %*% gamma.vec)
  }
  
  # Transform exposure and outcome using profile likehood weights
  Xnew <- X / new.sigma.vec
  Ynew <- Y / new.sigma.vec
  
  k <- length(invalid.idx)
  
  # Remove invalid IVs previously identified in step 2 if present
  if (k > 0) {
    X.sub <- Xnew[-invalid.idx,]
    Y.sub <- Ynew[-invalid.idx]
  } else {
    X.sub <- Xnew
    Y.sub <- Ynew
  }
  
  res <- susie(X.sub, Y.sub, L = 10, max = max.iter, intercept = FALSE)
  
  theta.vec <- unname(coef(res))[-1]
  prev.theta <- theta.vec
  iter <- 0
  diff <- 1e300
  thres <- tol
  while (diff > thres) {
    iter <- iter + 1
    
    gamma.vec <- c(unname(coef(res)[-1]), -1)
    
    for (i in 1:m.star) {
      new.sigma.vec[i] <- sqrt(t(gamma.vec) %*% Sigmai.list[[i]] %*% gamma.vec)
    }
    
    Xnew <- X / new.sigma.vec
    Ynew <- Y / new.sigma.vec
    
    if (k > 0) {
      X.sub <- Xnew[-invalid.idx,]
      Y.sub <- Ynew[-invalid.idx]
    } else {
      X.sub <- Xnew
      Y.sub <- Ynew
    }
    
    res <- susie(X.sub, Y.sub, L = L, max = max.iter, intercept = FALSE)
    theta.vec <- unname(coef(res))[-1]
    diff <- sum((prev.theta - theta.vec)^2)
    prev.norm <- diff
    if (diff > 0) {
      prev.theta <- theta.vec
    } else if (diff < 0) {
      iter <- iter - 1
      break
    }
    prev.norm <- diff
  }
  
  return(res)
}
