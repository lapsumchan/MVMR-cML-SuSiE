# Run from your existing MVMR-cML-SuSiE folder:
# Rscript --vanilla offline-full.R
# Or, in an updated checkout: Rscript --vanilla inst/scripts/offline-full.R
# Uses the already installed package; no installation or online retrieval.
if (dir.exists("local-library")) {
  .libPaths(c(normalizePath("local-library"), .libPaths()))
}
if (!requireNamespace("MVMRcMLSuSiE", quietly = TRUE)) {
  stop("MVMRcMLSuSiE is not installed. Run from the folder containing local-library, or install the development package first.",
       call. = FALSE)
}
library(MVMRcMLSuSiE)

run_offline_full <- function(output_root = "local-results") {
  output_dir <- file.path(output_root, paste0("MVMRcMLSuSiE-full-",
    format(Sys.time(), "%Y%m%d-%H%M%S"), "-", Sys.getpid()))
  if (!dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create a new output directory: ", output_dir, call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  connection <- file(file.path(output_dir, "run.log"), open = "wt")
  old_output <- sink.number()
  old_message <- sink.number(type = "message")
  sink(connection, split = TRUE)
  sink(connection, type = "message")
  on.exit({
    if (old_message == 2L) sink(type = "message") else
      sink(getConnection(old_message), type = "message")
    while (sink.number() > old_output) sink()
    close(connection)
  }, add = TRUE)
  assert <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
  }
  write_csv <- function(x, filename) {
    utils::write.csv(x, file.path(output_dir, filename), row.names = FALSE)
  }
  package_record <- function(package) {
    d <- utils::packageDescription(package)
    d[intersect(c("Package", "Version", "RemoteUsername", "RemoteRepo",
                  "RemoteSha", "GithubSHA1"), names(d))]
  }
  started <- Sys.time()
  script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  code_path <- system.file("R", "MVMRcMLSuSiE.rdb", package = "MVMRcMLSuSiE")
  results <- list(
    status = "RUNNING", completed_stage = "none", started = started,
    output_dir = output_dir, runner_revision = "full-offline-v1",
    scope = paste("Recompute all 249 step-1 screens; select using p < 0.05/27;",
                  "recompute all selected step-2 fits from bundled joint inputs;",
                  "run step 3 using only these newly computed step-2 results."),
    host = Sys.info(), session = utils::sessionInfo(),
    packages = setNames(lapply(c("MVMRcMLSuSiE", "MRcML", "susieR"), package_record),
                        c("MVMRcMLSuSiE", "MRcML", "susieR")),
    settings = list(screening_cutoff = 0.05 / 27, instrument_cutoff = 5e-8,
                    outcome_sample_size = 487511, seed = 1L,
                    MRcML_random_start = 100L, MRcML_random_seed = 1L,
                    S = 10L, max.iter = 200L, tol = 1e-10, max.outer.iter = 200L),
    installed_code_md5 = if (file.exists(code_path)) tools::md5sum(code_path) else NA,
    runner_md5 = if (length(script_path) == 1L && file.exists(script_path))
      tools::md5sum(script_path) else NA
  )
  checkpoint <- function() {
    saveRDS(results, file.path(output_dir, "results.rds"))
    flush(connection)
  }
  writeLines(capture.output(print(results$session)), file.path(output_dir, "sessionInfo.txt"))
  cat("MVMRcMLSuSiE complete offline pipeline\n")
  cat("Started: ", format(started, tz = "UTC", usetz = TRUE), "\n", sep = "")
  cat("Host: ", Sys.info()[["sysname"]], " / ", Sys.info()[["machine"]], "\n", sep = "")
  cat("R: ", R.version.string, "\n", sep = "")
  cat("Scope: ", results$scope, "\n", sep = "")
  cat("Results: ", output_dir, "\n", sep = "")
  print(results$packages)
  checkpoint()

  results <- tryCatch({
    assert("max.outer.iter" %in% names(formals(MVMRcMLSuSiE::mvmr.cml.susie.step3)),
           "Install the development package with the bounded outer loop first.")
    data_env <- new.env(parent = emptyenv())
    utils::data("example.dat", package = "MVMRcMLSuSiE", envir = data_env)
    dat <- data_env$example.dat
    n <- dat$step1$sample.sizes
    ids <- names(n)
    assert(length(ids) == 249L && !anyNA(ids) && !anyDuplicated(ids),
           "Expected 249 unique exposure IDs in the bundled sample sizes.")
    input_fields <- c("beta.exposure.ls", "se.exposure.ls", "beta.outcome.ls", "se.outcome.ls")
    assert(all(vapply(dat$step1[input_fields], length, integer(1)) == length(ids)),
           "The bundled step-1 lists do not match the 249-exposure axis.")
    assert(length(dat$step1$step1res) == length(ids),
           "Historical step-1 results do not match the exposure axis.")
    # The step-1 lists and outcome vectors are unnamed: retain their bundled order.
    raw <- dat$step2
    joint_ids <- colnames(raw$beta.exposure.mat)
    snp_ids <- rownames(raw$beta.exposure.mat)
    assert(length(joint_ids) == 43L && !anyNA(joint_ids) && !anyDuplicated(joint_ids) &&
             length(snp_ids) == 187L && !anyNA(snp_ids) && !anyDuplicated(snp_ids),
           "Expected the bundled joint input to contain 187 unique SNPs and 43 unique exposures.")
    assert(identical(dimnames(raw$beta.exposure.mat), dimnames(raw$se.exposure.mat)) &&
             identical(dimnames(raw$beta.exposure.mat), dimnames(raw$pval.exposure.mat)),
           "The joint exposure matrices have different axes.")
    assert(length(raw$beta.outcome.vec) == length(snp_ids) &&
             length(raw$se.outcome.vec) == length(snp_ids),
           "Outcome vectors do not match the joint SNP rows.")
    saved2 <- raw$step2res
    assert(identical(joint_ids, colnames(saved2$mvdat$exposure_beta)) &&
             identical(snp_ids, rownames(saved2$mvdat$exposure_beta)),
           "Historical step-2 results have different exposure/SNP axes.")
    for (field in c("beta", "se")) {
      assert(isTRUE(all.equal(raw[[paste0(field, ".exposure.mat")]],
                              saved2$mvdat[[paste0("exposure_", field)]], tolerance = 0)) &&
               isTRUE(all.equal(as.numeric(raw[[paste0(field, ".outcome.vec")]]),
                                as.numeric(saved2$mvdat[[paste0("outcome_", field)]]), tolerance = 0)),
             paste("Raw and historical step-2", field, "inputs differ."))
    }
    cat("PASS: bundled input dimensions and joint axes agree.\n")
    set.seed(results$settings$seed)
    results$rng_kind <- RNGkind()

    cat("\nSTEP 1: recomputing all 249 screens; this can take tens of minutes.\n")
    cat("Each completed exposure is checkpointed in results.rds.\n")
    tick <- proc.time()[["elapsed"]]
    p <- rep(NA_real_, length(ids))
    for (j in seq_along(ids)) {
      cat("Screen ", j, "/", length(ids), ": ", ids[j], " ...\n", sep = "")
      flush(connection)
      # Each package call already uses random_seed=1, so splitting the loop
      # preserves its calculations while allowing per-exposure checkpoints.
      p[j] <- MVMRcMLSuSiE::mvmr.cml.susie.step1(
        sample.sizes = c(n[j], results$settings$outcome_sample_size),
        beta.exposure.ls = dat$step1$beta.exposure.ls[j],
        se.exposure.ls = dat$step1$se.exposure.ls[j],
        beta.outcome.ls = dat$step1$beta.outcome.ls[j],
        se.outcome.ls = dat$step1$se.outcome.ls[j], use.openGWAS = FALSE)
      assert(length(p[j]) == 1L && is.finite(p[j]) && p[j] >= 0 && p[j] <= 1,
             paste("Step 1 returned an invalid p-value for", ids[j]))
      results$step1 <- list(p_values = setNames(p, ids), completed_exposures = j,
                            elapsed_seconds = proc.time()[["elapsed"]] - tick)
      checkpoint()
    }
    assert(length(p) == length(ids) && all(is.finite(p)) && all(p >= 0 & p <= 1),
           "Step 1 did not return 249 finite probabilities.")
    selected <- which(p < results$settings$screening_cutoff)
    selected_ids <- ids[selected]
    historical_selected <- dat$step1$step1res < results$settings$screening_cutoff
    results$step1 <- list(p_values = setNames(p, ids), completed_exposures = length(ids),
                          selected_ids = selected_ids,
                          elapsed_seconds = proc.time()[["elapsed"]] - tick)
    comparison1 <- data.frame(exposure = ids, p_value = p,
      historical_p_value = dat$step1$step1res, difference = p - dat$step1$step1res,
      selected = p < results$settings$screening_cutoff, historical_selected = historical_selected)
    results$comparison <- list(
      step1_maximum_absolute_p_difference = max(abs(comparison1$difference)),
      step1_selection_identical = identical(comparison1$selected, historical_selected),
      step1_changed_exposures = ids[comparison1$selected != historical_selected])
    write_csv(comparison1, "step1-comparison.csv")
    results$completed_stage <- "step1"
    checkpoint()
    cat("STEP 1 complete: ", length(selected_ids), " exposures selected in ",
        results$step1$elapsed_seconds, " seconds.\n", sep = "")
    # The bundle contains joint data for this selection only. Never substitute
    # saved screening decisions if the newly computed selection changes.
    assert(identical(selected_ids, joint_ids), paste(
      "New selection differs from the bundled 43-exposure joint input.",
      "Step-1 results were saved; matching joint data are needed to continue.",
      "No historical selection or fit has been substituted."))
    assert(all(selected_ids %in% rownames(dat$step3$metdrho)) &&
             all(selected_ids %in% colnames(dat$step3$metdrho)),
           "The covariance matrix is missing selected exposure IDs.")

    cat("\nSTEP 2: recomputing all ", length(selected_ids), " selected fits on all ",
        length(snp_ids), " joint-input rows.\n", sep = "")
    tick <- proc.time()[["elapsed"]]
    step2 <- MVMRcMLSuSiE::mvmr.cml.susie.step2(
      sample.sizes.subset = c(n[selected_ids], results$settings$outcome_sample_size),
      beta.exposure.mat = raw$beta.exposure.mat[, selected_ids, drop = FALSE],
      se.exposure.mat = raw$se.exposure.mat[, selected_ids, drop = FALSE],
      beta.outcome.vec = raw$beta.outcome.vec, se.outcome.vec = raw$se.outcome.vec,
      pval.exposure.mat = raw$pval.exposure.mat[, selected_ids, drop = FALSE],
      use.openGWAS = FALSE, cutoff = results$settings$instrument_cutoff
    )
    assert(length(step2$theta.vec) == length(selected_ids) && all(is.finite(step2$theta.vec)),
           "Step 2 returned malformed initial effects.")
    results$step2 <- list(result = step2, elapsed_seconds = proc.time()[["elapsed"]] - tick)
    results$comparison$step2_maximum_absolute_theta_difference <-
      max(abs(step2$theta.vec - saved2$theta.vec))
    results$comparison$step2_invalid_set_identical <- setequal(step2$invalid.idx, saved2$invalid.idx)
    results$comparison$step2_new_invalid_SNPs <- snp_ids[setdiff(step2$invalid.idx, saved2$invalid.idx)]
    results$comparison$step2_removed_invalid_SNPs <- snp_ids[setdiff(saved2$invalid.idx, step2$invalid.idx)]
    write_csv(data.frame(exposure = selected_ids, theta = step2$theta.vec,
                        historical_theta = saved2$theta.vec,
                        difference = step2$theta.vec - saved2$theta.vec), "step2-comparison.csv")
    write_csv(data.frame(row = seq_along(snp_ids), SNP = snp_ids,
      invalid = seq_along(snp_ids) %in% step2$invalid.idx,
      historical_invalid = seq_along(snp_ids) %in% saved2$invalid.idx), "instrument-comparison.csv")
    results$completed_stage <- "step2"
    checkpoint()
    cat("STEP 2 complete: ", length(step2$invalid.idx), " invalid instruments in ",
        results$step2$elapsed_seconds, " seconds.\n", sep = "")

    cat("\nSTEP 3: fitting from the NEW step-2 estimates and invalid-IV set.\n")
    rho <- diag(length(selected_ids) + 1L)
    dimnames(rho) <- list(c(selected_ids, "outcome_AD"), c(selected_ids, "outcome_AD"))
    rho[selected_ids, selected_ids] <- dat$step3$metdrho[selected_ids, selected_ids, drop = FALSE]
    tick <- proc.time()[["elapsed"]]
    fit <- MVMRcMLSuSiE::mvmr.cml.susie.step3(
      step2$mvdat, step2$invalid.idx, step2$theta.vec, rho,
      S = results$settings$S, max.iter = results$settings$max.iter,
      tol = results$settings$tol, max.outer.iter = results$settings$max.outer.iter)
    convergence <- attr(fit, "cml.convergence")
    assert(isTRUE(convergence$outer.converged) && isTRUE(convergence$inner.converged),
           "Step 3 did not confirm both outer and final inner convergence.")
    assert(identical(colnames(fit$alpha), selected_ids) && all(is.finite(fit$alpha)) &&
             all(is.finite(stats::coef(fit))) && length(fit$pip) == length(selected_ids) &&
             all(is.finite(fit$pip)), "Step 3 returned malformed results.")
    saved3 <- dat$step3$step3res
    assert(identical(dim(fit$alpha), dim(saved3$alpha)) &&
             identical(colnames(fit$alpha), colnames(saved3$alpha)) &&
             length(stats::coef(fit)) == length(stats::coef(saved3)),
           "Historical and newly fitted step-3 result axes differ.")
    results$step3 <- list(result = fit, rho = rho, convergence = convergence,
                          elapsed_seconds = proc.time()[["elapsed"]] - tick)
    results$comparison$step3_maximum_absolute_alpha_difference <- max(abs(fit$alpha - saved3$alpha))
    results$comparison$step3_maximum_absolute_coefficient_difference <-
      max(abs(stats::coef(fit) - stats::coef(saved3)))
    results$comparison$step3_maximum_absolute_PIP_difference <- max(abs(fit$pip - saved3$pip))
    write_csv(data.frame(exposure = selected_ids, PIP = fit$pip,
      historical_PIP = saved3$pip, difference = fit$pip - saved3$pip), "step3-PIP-comparison.csv")
    write_csv(data.frame(term = c("(Intercept)", selected_ids),
      coefficient = as.numeric(stats::coef(fit)), historical_coefficient = as.numeric(stats::coef(saved3)),
      difference = as.numeric(stats::coef(fit) - stats::coef(saved3))), "step3-coefficient-comparison.csv")
    alpha_axes <- expand.grid(component = seq_len(nrow(fit$alpha)), exposure = selected_ids,
                              stringsAsFactors = FALSE)
    write_csv(data.frame(alpha_axes, alpha = as.vector(fit$alpha),
      historical_alpha = as.vector(saved3$alpha), difference = as.vector(fit$alpha - saved3$alpha)),
      "step3-alpha-comparison.csv")
    results$completed_stage <- "step3"
    results$status <- "PASS"
    cat("STEP 3 complete.\n")
    print(convergence)
    cat("\nHistorical comparisons (numeric differences are descriptive):\n")
    print(results$comparison)
    cat("PASS means the complete bundled offline computation finished with finite outputs\n",
        "and confirmed final convergence. It is not a scientific calibration test.\n",
        "This run uses pre-harmonized bundled inputs and excludes online retrieval\n",
        "and the separate post-selection MVMR-cML analysis.\n", sep = "")
    results
  }, error = function(e) {
    results$status <- "FAIL"
    results$error <- conditionMessage(e)
    cat("FAIL: ", conditionMessage(e), "\n", sep = "")
    results
  })
  results$finished <- Sys.time()
  results$total_elapsed_seconds <- as.numeric(difftime(results$finished, started, units = "secs"))
  results$session <- utils::sessionInfo()
  writeLines(capture.output(print(results$session)), file.path(output_dir, "sessionInfo.txt"))
  checkpoint()
  summary <- c(paste("Overall status:", results$status),
               paste("Last completed stage:", results$completed_stage),
               paste("Total elapsed seconds:", round(results$total_elapsed_seconds, 3)),
               capture.output(print(results$comparison)),
               if (!is.null(results$error)) paste("Error:", results$error),
               paste("Saved logs, results and comparison CSV files in:", output_dir))
  writeLines(summary, file.path(output_dir, "summary.txt"))
  cat("\n", paste(summary, collapse = "\n"), "\n", sep = "")
  invisible(results)
}

.cml_full_result <- run_offline_full()
if (!identical(.cml_full_result$status, "PASS")) {
  stop("Complete offline run stopped. See ", file.path(.cml_full_result$output_dir, "run.log"),
       call. = FALSE)
}
invisible(.cml_full_result)
