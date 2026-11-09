# Run from the downloaded/cloned repository:
# Rscript --vanilla scripts/install-local.R
# Installs to a project-local library, then runs the bundled offline check.
if (getRversion() < "4.3.0") stop("Install R 4.3.0 or newer before running this script.")
if (!file.exists("DESCRIPTION") ||
    read.dcf("DESCRIPTION")[1L, "Package"] != "MVMRcMLSuSiE") {
  stop("Run this script from the MVMR-cML-SuSiE repository directory.")
}
local_library <- file.path(getwd(), "local-library")
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_library, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 300)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = local_library)
}
remotes::install_local(".", lib = local_library, dependencies = NA,
                       upgrade = "never", force = TRUE, build_vignettes = FALSE)
source(system.file("scripts", "offline-smoke.R", package = "MVMRcMLSuSiE"))
