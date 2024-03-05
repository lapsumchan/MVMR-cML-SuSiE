mvmr.cml.susie.step1 <- function(exposure.ids, outcome.id, sample.sizes) {
  
  L <- length(exposure.ids)
  pval.vec <- rep(NA, L)
  
  for (i in 1:L) {
    print(i)
    # Get instruments
    exposure.dat <- extract_instruments(exposure.ids[i])
    
    # Get effects of instruments on outcome
    outcome.dat <- extract_outcome_data(snps = exposure.dat$SNP, outcomes = outcome.id)
    
    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure.dat, outcome.dat)
    
    n <- min(sample_sizes[i], outcome.dat$samplesize.outcome)
    
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
  return(pval.vec)
}

mvmr.cml.susie.step2 <- function(exposure.ids.subset, outcome.id, sample.sizes.subset) {
  
  # Get instruments jointly this time
  exposure.dat <- mv_extract_exposures(exposure.ids.subset)
  
  # Get effects of instruments on outcome
  outcome.dat <- extract_outcome_data(snps = exposure.dat$SNP, outcomes = outcome.id)
  
  # Harmonise the exposure and outcome data
  mvdat <- mv_harmonise_data(exposure.dat, outcome.dat)
  
  L.star <- length(exposure.ids.subset)
  
  uvmr.ls <- list()
  
  # Focus only on QTL for each exposure
  for (i in 1:L.star) {
    uvmr.ls[[i]] <- which(mvdat$exposure_pval[,i] < 5e-8)
  }
  
  pval.vec <- rep(NA, L.star)
  theta.vec <- rep(NA, L.star)
  for (i in 1:L.star) {
    
    n <- min(sample.sizes.subset[i], outcome.dat$samplesize.outcome)
    
    # Perform UVMR-cML
    cML.result <- mr_cML(mvdat$exposure_beta[uvmr.ls[[i]],i],
                         mvdat$outcome_beta[uvmr.ls[[i]]],
                         mvdat$exposure_se[uvmr.ls[[i]],i],
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
  
  X <- mvdat$exposure_beta
  Y <- mvdat$outcome_beta
  SigmaX <- mvdat$exposure_se
  seY <- mvdat$outcome_se
  
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
    }
    
    res <- susie(X.sub, Y.sub, L = L, max = max.iter, intercept = FALSE)
    theta.vec <- unname(coef(res))[-1]
    diff <- sum((prev.theta - theta.vec)^2)
    prev.norm <- diff
    if (diff > 0) {
      prev.theta <- theta.vec2
    } else if (diff < 0) {
      iter <- iter - 1
      break
    }
    prev.norm <- diff
  }
  return(res)
}