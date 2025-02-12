#' Step 1 of MVMR-cML-SuSiE
#'
#' This function performs the first step of the MVMR-cML-SuSiE approach to obtain each exposure's univariable
#' Mendelian randomization (UVMR) p-value, used for initial filtering of exposures. You can choose to pull
#' GWAS data using \code{OpenGWAS} (by setting \code{use.openGWAS = TRUE}) or provide your own summary
#' statistics (by setting \code{use.openGWAS = FALSE})
#'
#' @param exposure.ids A length \eqn{L} character vector of openGWAS IDs corresponding to the exposures,
#'   where \eqn{L} is the number of intended exposures
#' @param outcome.id A character corresponding to the outcome openGWAS ID
#' @param sample.sizes A numeric vector of sample sizes corresponding to each exposure and outcome
#'   If \code{use.openGWAS = TRUE}, the length should be \eqn{L} (exposures only).
#'   Otherwise, the length should be \eqn{L + 1}, with the last element for the outcome sample size
#' @param beta.exposure.ls A list (length \eqn{L}) of numeric vectors with effect sizes (betas) for each SNP on exposure
#' @param se.exposure.ls A list (length \eqn{L}) of numeric vectors with standard errors for each SNP on exposure
#' @param beta.outcome.ls A list (length \eqn{L}) of numeric vectors with effect sizes (betas) for each SNP on outcome
#' @param se.outcome.ls A list (length \eqn{L}) of numeric vectors with standard errors for each SNP on outcome
#' @param use.openGWAS A logical indicating whether to use \code{OpenGWAS} to extract GWAS data based on
#'   \code{exposure.ids} and \code{outcome.id}. If \code{FALSE}, requires user-provided summary statistics
#'   (lists for \code{beta.exposure.ls}, \code{se.exposure.ls}, \code{beta.outcome.ls}, and \code{se.outcome.ls})
#'   Default is \code{TRUE}
#'
#' @return A numeric vector (length \eqn{L}) of p-values, each corresponding to the UVMR p-value for one exposure
#'
#' @export
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

#' Step 2 of MVMR-cML-SuSiE
#'
#' This function performs the second step of the MVMR-cML-SuSiE approach. After identifying
#' exposures of interest in Step 1 (e.g., by filtering UVMR p-values), step 2 extracts the
#' relevant instruments for those exposures, identifies invalid instruments (if any),
#' and obtains an initial estimate from univariable Mendelian randomization (UVMR) for each exposure
#'
#' @param exposure.ids.subset A length \eqn{L^{*}} character vector of exposure IDs (a subset from step 1)
#' @param outcome.id A character specifying the outcome ID
#' @param sample.sizes.subset A length \eqn{L^{*}} numeric vector of sample sizes corresponding to
#'   \code{exposure.ids.subset} (and optionally the outcome as the last element if user-provided data)
#' @param beta.exposure.mat A size \eqn{m^{*} \times L^{*}} matrix of exposure effect sizes, where
#'   \eqn{m^{*}} is the number of SNPs
#' @param se.exposure.mat A size \eqn{m^{*} \times L^{*}} matrix of exposure standard errors
#' @param beta.outcome.vec A length \eqn{m^{*}} numeric vector of outcome effect sizes for the same SNPs
#' @param se.outcome.vec A length \eqn{m^{*}} numeric vector of standard errors for \code{beta.outcome.vec}
#' @param pval.exposure.mat A \eqn{m^{*} \times L^{*}} matrix of exposure p-values
#' @param use.openGWAS A logical indicating whether to extract data via \code{OpenGWAS}. If \code{FALSE},
#'   you must supply the \code{*.mat} and \code{*.vec} arguments. Default is \code{TRUE}
#' @param cutoff A numeric threshold for instrument selection. Default is \eqn{5 \times 10^{-8}}
#'
#' @return A list containing:
#' \item{mvdat}{A list containing the relevant data (exposure/outcome betas, SEs, etc)}
#' \item{invalid.idx}{A vector of indices of invalid instruments identified based off the UVMR procedure}
#' \item{theta.vec}{A numeric vector of initial values for each exposure}
#'
#' @export
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

#' Step 3 of MVMR-cML-SuSiE (Iterative SuSiE)
#'
#' This function performs the iterative SuSiE algorithm to identify potential exposure signal clusters
#'
#' @param mvdat A list returned by \code{mvmr.cml.susie.step2}, containing exposure and outcome summary statistics
#' @param invalid.idx An integer vector of invalid IV indices as identified in step 2
#' @param theta.vec A numeric vector of initial exposure effect estimates from step 2
#' @param rho.mat A genetic correlation matrix of size \eqn{(L^{*} + 1) \times (L^{*} + 1)}, where \eqn{L^{*}} is
#'   the number of exposures. The last row/column corresponds to the outcome
#' @param S The number of single effects (signal clusters) in SuSiE. Default is 10.
#' @param max.iter The maximum number of iterations for SuSiE. Default is 200.
#' @param tol A numeric value indicating the convergence threshold for the iterative procedure. Default is \eqn{1 \times 10^{-10}}.
#' @return A fitted \code{susie} object. The key SuSiE output is the posterior inclusion probability (PIP) matrix \code{alpha}
#'   used for determining signal clusters
#'
#' @export
mvmr.cml.susie.step3 <- function(mvdat, invalid.idx, theta.vec, rho.mat,
                                 S = 10, max.iter = 200, tol = 1e-10) {

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

  res <- susie(X.sub, Y.sub, L = S, max_iter = max.iter, intercept = FALSE)

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

    res <- susie(X.sub, Y.sub, L = S, max_iter = max.iter, intercept = FALSE)
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
