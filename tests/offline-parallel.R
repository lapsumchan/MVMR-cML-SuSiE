# Focused serial/PSOCK equivalence and warning/error transport checks.
script <- system.file("scripts", "offline-full.R", package = "MVMRcMLSuSiE")
stopifnot(nzchar(script))
runner <- new.env(parent = globalenv())
for (expr in parse(script)) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      identical(expr[[2L]], as.name(".cml_args"))) break
  eval(expr, runner)
}

check_parallel <- function() {
  data_env <- new.env()
  utils::data("example.dat", package = "MVMRcMLSuSiE", envir = data_env)
  dat <- data_env$example.dat
  # Different real exposures with small instrument sets keep this check short.
  selected <- order(lengths(dat$step1$beta.exposure.ls))[seq_len(6L)]
  jobs <- lapply(selected, function(j) list(index = j,
    exposure = names(dat$step1$sample.sizes)[j],
    sample_sizes = c(dat$step1$sample.sizes[j], 487511),
    bx = dat$step1$beta.exposure.ls[[j]], sx = dat$step1$se.exposure.ls[[j]],
    by = dat$step1$beta.outcome.ls[[j]], sy = dat$step1$se.outcome.ls[[j]]))
  serial <- lapply(jobs, runner$.cml_capture_task, task_fun = runner$.cml_step1_task)
  cluster <- parallel::makePSOCKcluster(2L, rscript_args = "--vanilla", outfile = tempfile())
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, runner$.cml_init_worker, .libPaths(), RNGkind())
  parallel_results <- parallel::parLapply(cluster, jobs, runner$.cml_capture_task,
                                          task_fun = runner$.cml_step1_task)
  for (i in seq_along(jobs)) {
    a <- serial[[i]]
    b <- parallel_results[[i]]
    stopifnot(is.null(a$error), is.null(b$error), identical(a$index, b$index),
              identical(a$exposure, b$exposure), identical(a$value, b$value),
              identical(a$warnings, b$warnings))
  }
  cat("PASS: six real exposure fits agree exactly between serial and two PSOCK workers.\n")

  # Deliberately emit more than R's default 50-warning limit, then fail one job.
  noisy_task <- function(job) {
    for (i in seq_len(55L)) warning("intentional transport-test warning", call. = FALSE)
    if (job$index == 2L) stop("intentional transport-test error", call. = FALSE)
    job$index / 10
  }
  environment(noisy_task) <- baseenv()
  noisy_jobs <- list(list(index = 1L, exposure = "fixture_a"),
                     list(index = 2L, exposure = "fixture_b"))
  noisy <- parallel::parLapply(cluster, noisy_jobs, runner$.cml_capture_task,
                               task_fun = noisy_task)
  stopifnot(length(noisy) == 2L, is.null(noisy[[1L]]$error),
            identical(noisy[[1L]]$value, 0.1),
            grepl("intentional transport-test error", noisy[[2L]]$error),
            nrow(noisy[[1L]]$warnings) == 55L, nrow(noisy[[2L]]$warnings) == 55L,
            all(noisy[[1L]]$warnings$exposure == "fixture_a"),
            all(noisy[[2L]]$warnings$exposure == "fixture_b"))
  cat("PASS: all 110 warnings survived worker transport, including warnings before a fitting error.\n")
}

check_parallel()
