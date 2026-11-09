# Run against an installed package with Rscript --vanilla.
# These are small execution/regression checks, not a scientific validation study.
library(MVMRcMLSuSiE)
stopifnot(!"package:MRcML" %in% search(),
          !"package:TwoSampleMR" %in% search())

expect_error <- function(expr, pattern = NULL) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(error, "error")) stop("Expected an error, but the call succeeded.")
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(error), fixed = TRUE)) {
    stop("Unexpected error: ", conditionMessage(error))
  }
  invisible(error)
}

# Deterministic, small summary-statistic inputs with two distinct exposures.
i <- seq_len(12L)
bx <- cbind(exposure_A = 0.08 + 0.018 * i,
            exposure_B = 0.11 + 0.013 * i + 0.025 * cos(i))
sx <- cbind(0.008 + 0.0002 * i, 0.011 + 0.0001 * i)
colnames(sx) <- colnames(bx)
by <- 0.4 * bx[, 1L] + 0.005 * sin(2 * i)
sy <- rep(0.02, length(i))
n <- c(20000, 22000, 25000)
px <- matrix(1e-10, nrow(bx), ncol(bx), dimnames = dimnames(bx))

reference <- lapply(seq_len(ncol(bx)), function(j) {
  MRcML::mr_cML(bx[, j], by, sx[, j], sy,
                n = min(n[j], n[3L]), random_start = 100, random_seed = 1)
})
stopifnot(all(vapply(reference, function(x) {
  is.finite(x$MA_BIC_p) && is.finite(x$MA_BIC_theta)
}, logical(1))))

p <- MVMRcMLSuSiE::mvmr.cml.susie.step1(
  sample.sizes = n,
  beta.exposure.ls = lapply(seq_len(ncol(bx)), function(j) bx[, j]),
  se.exposure.ls = lapply(seq_len(ncol(bx)), function(j) sx[, j]),
  beta.outcome.ls = rep(list(by), ncol(bx)),
  se.outcome.ls = rep(list(sy), ncol(bx)),
  use.openGWAS = FALSE
)
stopifnot(isTRUE(all.equal(as.numeric(p),
  vapply(reference, function(x) x$MA_BIC_p, numeric(1)), tolerance = 1e-10)))

stage2 <- MVMRcMLSuSiE::mvmr.cml.susie.step2(
  sample.sizes.subset = n,
  beta.exposure.mat = bx, se.exposure.mat = sx,
  beta.outcome.vec = by, se.outcome.vec = sy,
  pval.exposure.mat = px, use.openGWAS = FALSE
)
expected_invalid <- sort(unique(unlist(lapply(reference, function(x) x$BIC_invalid))))
stopifnot(
  identical(stage2$mvdat$exposure_beta, bx),
  isTRUE(all.equal(as.numeric(stage2$theta.vec),
    vapply(reference, function(x) x$MA_BIC_theta, numeric(1)), tolerance = 1e-10)),
  identical(as.integer(stage2$invalid.idx), as.integer(expected_invalid)),
  !"package:MRcML" %in% search()
)

# The row exclusion must retain a matrix even with just one exposure.
single <- list(exposure_beta = bx[, 1L, drop = FALSE],
               exposure_se = sx[, 1L, drop = FALSE],
               outcome_beta = by, outcome_se = sy)
fit <- MVMRcMLSuSiE::mvmr.cml.susie.step3(
  single, invalid.idx = 2L, theta.vec = 0.4,
  rho.mat = diag(2L), S = 1L, max.outer.iter = 200L
)
stopifnot(inherits(fit, "susie"), ncol(fit$alpha) == 1L,
          all(is.finite(fit$alpha)), all(is.finite(stats::coef(fit))))

# Force nonconvergence in an isolated copy of the function, without changing
# either package namespace. The real outer loop must stop after its first
# update; the alternating mock makes its squared coefficient change nonzero.
bounded_step3 <- MVMRcMLSuSiE::mvmr.cml.susie.step3
mock_env <- new.env(parent = environment(bounded_step3))
mock_env$calls <- 0L
mock_env$.test_susie <- function(X, Y, ...) {
  stopifnot(is.matrix(X))
  mock_env$calls <- mock_env$calls + 1L
  list(test_coef = c(0, rep(if (mock_env$calls %% 2L) 0.1 else 0.8, ncol(X))),
       converged = TRUE)
}
mock_env$.test_coef <- function(object, ...) object$test_coef
replace_calls <- function(x) {
  if (identical(x, quote(susieR::susie))) return(as.name(".test_susie"))
  if (identical(x, quote(stats::coef))) return(as.name(".test_coef"))
  if (is.call(x)) return(as.call(lapply(as.list(x), replace_calls)))
  x
}
body(bounded_step3) <- replace_calls(body(bounded_step3))
environment(bounded_step3) <- mock_env
expect_error(bounded_step3(single, integer(0), 0.4, diag(2L),
                          S = 1L, tol = 1e-10, max.outer.iter = 1L),
             "Outer iteration limit reached")
stopifnot(mock_env$calls == 2L)

# Bad indices, covariance dimensions and nonpositive sampling SEs must fail
# before numerical fitting. Empty instrument sets must fail explicitly.
expect_error(MVMRcMLSuSiE::mvmr.cml.susie.step3(single, 13L, 0.4, diag(2L)))
expect_error(MVMRcMLSuSiE::mvmr.cml.susie.step3(single, integer(0), 0.4, diag(3L)))
bad <- single
bad$outcome_se[1L] <- 0
expect_error(MVMRcMLSuSiE::mvmr.cml.susie.step3(bad, integer(0), 0.4, diag(2L)))
expect_error(MVMRcMLSuSiE::mvmr.cml.susie.step2(
  sample.sizes.subset = n, beta.exposure.mat = bx, se.exposure.mat = sx,
  beta.outcome.vec = by, se.outcome.vec = sy,
  pval.exposure.mat = matrix(1, nrow(bx), ncol(bx)), use.openGWAS = FALSE
))

cat("PASS: offline steps 1-2, one-exposure step 3, outer cap, and input errors.\n")
