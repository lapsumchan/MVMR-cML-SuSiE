# Internal checks for the existing estimators; no data are dropped or repaired.
.check_mr_vectors <- function(bx, by, sx, sy, n, label) {
  values <- list(bx, by, sx, sy)
  if (!all(vapply(values, is.numeric, logical(1))) ||
      length(bx) < 2L || any(lengths(values) != length(bx)) ||
      !all(vapply(values, function(x) all(is.finite(x)), logical(1))) ||
      any(sx <= 0) || any(sy <= 0)) {
    stop(paste0(label, ": supply at least two aligned instruments, finite effects, and positive finite SEs."),
         call. = FALSE)
  }
  if (length(n) != 1L || !is.finite(n) || n <= 1) {
    stop(paste0(label, ": sample size must be finite and greater than one."), call. = FALSE)
  }
}

.check_positive_integer <- function(x, label) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != floor(x)) {
    stop(paste0(label, " must be a positive integer."), call. = FALSE)
  }
}

.check_profile_sd <- function(x) {
  if (any(!is.finite(x)) || any(x <= 0)) {
    stop("Profile variances must be positive and finite; check rho.mat, SEs, and theta.vec.",
         call. = FALSE)
  }
}
