# Small offline execution check for the installed MVMRcMLSuSiE package.
# Run in a fresh R session. This does not repeat all 249 screening analyses,
# perform downstream MVMR-cML inference, or establish new scientific results.
library(MVMRcMLSuSiE)

run_offline_smoke <- function(output_root = "local-results") {
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  output_dir <- file.path(output_root, paste0("MVMRcMLSuSiE-", stamp, "-", Sys.getpid()))
  if (!dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create a new output directory: ", output_dir, call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  log_path <- file.path(output_dir, "run.log")
  log_connection <- file(log_path, open = "wt")
  previous_output_sinks <- sink.number()
  previous_message_sink <- sink.number(type = "message")
  sink(log_connection, split = TRUE)
  sink(log_connection, type = "message")
  on.exit({
    if (previous_message_sink == 2L) sink(type = "message") else
      sink(getConnection(previous_message_sink), type = "message")
    while (sink.number() > previous_output_sinks) sink()
    close(log_connection)
  }, add = TRUE)

  assert <- function(condition, message) {
    if (!isTRUE(condition)) stop(message, call. = FALSE)
  }
  package_record <- function(package) {
    description <- utils::packageDescription(package)
    description[intersect(c("Package", "Version", "RemoteUsername", "RemoteRepo",
                            "RemoteSha", "GithubSHA1"), names(description))]
  }
  started <- Sys.time()
  results <- list(
    status = "RUNNING", started = started, output_dir = output_dir,
    scope = paste("One bundled exposure through step 1; one bundled exposure through step 2;",
                  "full 43-exposure step 3 from saved upstream stage-2 results."),
    source_baseline = list(repository = "lapsumchan/MVMR-cML-SuSiE",
                           upstream_commit = "da71c5cd62bda8b0e89f8732ca6013585fb9584b"),
    host = list(system = Sys.info(), R_version = R.version.string,
                R_platform = R.version$platform, R_architecture = R.version$arch),
    packages = lapply(c("MVMRcMLSuSiE", "MRcML", "susieR"), package_record)
  )
  names(results$packages) <- c("MVMRcMLSuSiE", "MRcML", "susieR")
  cat("MVMRcMLSuSiE offline smoke run\n")
  cat("Started: ", format(started, tz = "UTC", usetz = TRUE), "\n", sep = "")
  cat("Host: ", Sys.info()[["sysname"]], " / ", Sys.info()[["machine"]], "\n", sep = "")
  cat("R: ", R.version.string, " / ", R.version$platform, "\n", sep = "")
  cat("Scope: ", results$scope, "\n", sep = "")
  cat("Results: ", output_dir, "\n", sep = "")
  print(results$packages)

  results <- tryCatch({
    assert(!"package:MRcML" %in% search() && !"package:TwoSampleMR" %in% search(),
           "Restart R without attaching MRcML or TwoSampleMR, then run this script again.")
    assert("max.outer.iter" %in% names(formals(MVMRcMLSuSiE::mvmr.cml.susie.step3)),
           "This smoke script requires the local-runtime patch with max.outer.iter.")
    data_env <- new.env(parent = emptyenv())
    utils::data("example.dat", package = "MVMRcMLSuSiE", envir = data_env)
    dat <- data_env$example.dat
    n_exposure <- dat$step1$sample.sizes
    exposure_ids <- names(n_exposure)
    assert(length(exposure_ids) == 249L && !anyDuplicated(exposure_ids),
           "The bundled 249-exposure identity vector is missing or duplicated.")
    assert(length(dat$step1$step1res) == length(exposure_ids),
           "Saved screening results and exposure identities have different lengths.")
    selected <- which(dat$step1$step1res < 0.05 / 27)
    selected_ids <- exposure_ids[selected]
    saved_stage2 <- dat$step2$step2res
    assert(identical(selected_ids, colnames(saved_stage2$mvdat$exposure_beta)),
           "Saved screening identities/order do not match the bundled stage-2 columns.")
    assert(identical(selected_ids, colnames(dat$step2$beta.exposure.mat)) &&
             identical(dimnames(dat$step2$beta.exposure.mat),
                       dimnames(dat$step2$se.exposure.mat)) &&
             identical(dimnames(dat$step2$beta.exposure.mat),
                       dimnames(dat$step2$pval.exposure.mat)),
           "Bundled stage-2 input matrices have inconsistent instrument/exposure axes.")
    for (field in c("beta", "se")) {
      assert(isTRUE(all.equal(dat$step2[[paste0(field, ".exposure.mat")]],
                              saved_stage2$mvdat[[paste0("exposure_", field)]],
                              tolerance = 0)),
             paste("Raw and saved stage-2 exposure", field, "inputs differ."))
      assert(isTRUE(all.equal(as.numeric(dat$step2[[paste0(field, ".outcome.vec")]]),
                              as.numeric(saved_stage2$mvdat[[paste0("outcome_", field)]]),
                              tolerance = 0)),
             paste("Raw and saved stage-2 outcome", field, "inputs differ."))
    }
    assert(all(selected_ids %in% rownames(dat$step3$metdrho)) &&
             all(selected_ids %in% colnames(dat$step3$metdrho)),
           "The covariance input is missing one or more selected exposure identities.")
    results$inputs <- list(selected_exposures = selected_ids,
                            stage2_instruments = rownames(saved_stage2$mvdat$exposure_beta),
                            selection_threshold = 0.05 / 27,
                            outcome_sample_size = 487511L)
    cat("PASS: bundled exposure identities/order and saved/raw stage-2 inputs agree.\n")

    # Exercise step 1 on the smallest complete bundled exposure with >=3 IVs.
    iv_counts <- lengths(dat$step1$beta.exposure.ls)
    eligible <- which(iv_counts >= 3L)
    assert(length(eligible) > 0L, "No bundled step-1 exposure has at least three IVs.")
    j <- eligible[order(iv_counts[eligible], exposure_ids[eligible])[1L]]
    tick <- proc.time()[["elapsed"]]
    p <- MVMRcMLSuSiE::mvmr.cml.susie.step1(
      sample.sizes = c(unname(n_exposure[j]), 487511L),
      beta.exposure.ls = dat$step1$beta.exposure.ls[j],
      se.exposure.ls = dat$step1$se.exposure.ls[j],
      beta.outcome.ls = dat$step1$beta.outcome.ls[j],
      se.outcome.ls = dat$step1$se.outcome.ls[j], use.openGWAS = FALSE
    )
    assert(length(p) == 1L && is.finite(p) && p >= 0 && p <= 1,
           "Step 1 did not return a finite probability.")
    results$step1 <- list(exposure = exposure_ids[j], index = j,
                          instruments = iv_counts[j], p_value = as.numeric(p),
                          saved_p_value = dat$step1$step1res[j],
                          difference_from_saved = as.numeric(p) - dat$step1$step1res[j],
                          elapsed_seconds = proc.time()[["elapsed"]] - tick)
    cat("PASS: step 1, one complete bundled exposure.\n")
    print(results$step1)

    # Exercise step 2 on one exposure, keeping every row of its joint input.
    eligible_counts <- colSums(dat$step2$pval.exposure.mat < 5e-8)
    eligible <- which(eligible_counts >= 3L)
    assert(length(eligible) > 0L, "No bundled step-2 exposure has three qualifying IVs.")
    k <- eligible[order(eligible_counts[eligible], selected_ids[eligible])[1L]]
    tick <- proc.time()[["elapsed"]]
    small_stage2 <- MVMRcMLSuSiE::mvmr.cml.susie.step2(
      sample.sizes.subset = c(unname(n_exposure[selected_ids[k]]), 487511L),
      beta.exposure.mat = dat$step2$beta.exposure.mat[, k, drop = FALSE],
      se.exposure.mat = dat$step2$se.exposure.mat[, k, drop = FALSE],
      beta.outcome.vec = dat$step2$beta.outcome.vec,
      se.outcome.vec = dat$step2$se.outcome.vec,
      pval.exposure.mat = dat$step2$pval.exposure.mat[, k, drop = FALSE],
      use.openGWAS = FALSE
    )
    assert(length(small_stage2$theta.vec) == 1L && is.finite(small_stage2$theta.vec),
           "Step 2 did not return a finite initial effect.")
    qualifying_rows <- which(dat$step2$pval.exposure.mat[, k] < 5e-8)
    assert(all(small_stage2$invalid.idx %in% qualifying_rows),
           "Step-2 invalid-IV indices are not mapped to the original input rows.")
    results$step2 <- list(exposure = selected_ids[k], index_within_selected = k,
                          input_rows = nrow(dat$step2$beta.exposure.mat),
                          qualifying_instruments = eligible_counts[k],
                          theta = as.numeric(small_stage2$theta.vec),
                          saved_theta = as.numeric(saved_stage2$theta.vec[k]),
                          difference_from_saved = as.numeric(small_stage2$theta.vec) -
                            as.numeric(saved_stage2$theta.vec[k]),
                          invalid_rows = small_stage2$invalid.idx,
                          elapsed_seconds = proc.time()[["elapsed"]] - tick)
    cat("PASS: step 2, one exposure with all joint-input rows retained.\n")
    print(results$step2)

    # Recompute full step 3 from saved stage 2. This is deliberately separate
    # from the two small checks above; it uses the original 43-exposure fit.
    rho <- diag(length(selected_ids) + 1L)
    dimnames(rho) <- list(c(selected_ids, "outcome_AD"), c(selected_ids, "outcome_AD"))
    rho[seq_along(selected_ids), seq_along(selected_ids)] <-
      dat$step3$metdrho[selected_ids, selected_ids, drop = FALSE]
    tick <- proc.time()[["elapsed"]]
    fit <- MVMRcMLSuSiE::mvmr.cml.susie.step3(
      saved_stage2$mvdat, saved_stage2$invalid.idx, saved_stage2$theta.vec, rho
    )
    convergence <- attr(fit, "cml.convergence")
    assert(isTRUE(convergence$outer.converged) && isTRUE(convergence$inner.converged),
           "Step 3 did not confirm both outer and final inner convergence.")
    assert(identical(dim(fit$alpha), c(10L, length(selected_ids))) &&
             identical(colnames(fit$alpha), selected_ids) && all(is.finite(fit$alpha)) &&
             all(is.finite(stats::coef(fit))),
           "Step 3 returned unexpected dimensions, identities, or non-finite results.")
    saved_fit <- dat$step3$step3res
    assert(identical(dim(fit$alpha), dim(saved_fit$alpha)) &&
             identical(colnames(fit$alpha), colnames(saved_fit$alpha)),
           "Saved and recomputed step-3 alpha arrays do not have the same axes.")
    results$step3 <- list(exposures = length(selected_ids),
                          input_instruments = nrow(saved_stage2$mvdat$exposure_beta),
                          excluded_instruments = length(saved_stage2$invalid.idx),
                          convergence = convergence,
                          maximum_absolute_alpha_difference = max(abs(fit$alpha - saved_fit$alpha)),
                          maximum_absolute_coefficient_difference =
                            max(abs(stats::coef(fit) - stats::coef(saved_fit))),
                          elapsed_seconds = proc.time()[["elapsed"]] - tick)
    cat("PASS: full step 3 from saved stage 2; inner and outer convergence confirmed.\n")
    print(results$step3)
    cat("Differences from historical saved values are descriptive; no exact cross-version\n",
        "reproduction or scientific calibration claim is made by these execution checks.\n", sep = "")
    results$recomputed_step2 <- small_stage2
    results$recomputed_step3 <- fit
    results$status <- "PASS"
    results
  }, error = function(error) {
    results$status <- "FAIL"
    results$error <- conditionMessage(error)
    cat("FAIL: ", conditionMessage(error), "\n", sep = "")
    results
  })
  results$finished <- Sys.time()
  results$total_elapsed_seconds <- as.numeric(difftime(results$finished, started, units = "secs"))
  results$session <- utils::sessionInfo()
  writeLines(capture.output(print(results$session)), file.path(output_dir, "sessionInfo.txt"))
  saveRDS(results, file.path(output_dir, "results.rds"))
  cat("\nOverall status: ", results$status, "\n", sep = "")
  cat("Total elapsed seconds: ", results$total_elapsed_seconds, "\n", sep = "")
  cat("Saved run.log, sessionInfo.txt, and results.rds in:\n", output_dir, "\n", sep = "")
  invisible(results)
}

.cml_smoke_result <- run_offline_smoke()
if (!identical(.cml_smoke_result$status, "PASS")) {
  stop("Offline smoke run failed. See ", file.path(.cml_smoke_result$output_dir, "run.log"),
       call. = FALSE)
}
invisible(.cml_smoke_result)
