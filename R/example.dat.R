#' Example data for MVMR-cML-SuSiE
#'
#' This dataset, \code{example.dat}, is a list of lists containing a real data example based on 249
#' metabolomics GWAS from the UK Biobank (UKB) and an Alzheimer's disease (AD) GWAS by Bellenguez et al. (2022).
#' It demonstrates how MVMR-cML-SuSiE can be run entirely from user-supplied data (without \code{OpenGWAS})
#' and also contains precomputed outputs from each of the three steps of MVMR-cML-SuSiE
#'
#' \strong{Structure:}
#'
#' \itemize{
#'   \item \code{step1} : A list relevant to step 1
#'     \itemize{
#'       \item \code{sample.sizes}    - A numeric vector of sample sizes for \eqn{L = 249} exposures.
#'       \item \code{beta.exposure.ls} - A list of length \eqn{L}, each an exposure effect-size vector
#'         (used if \code{use.openGWAS = FALSE})
#'       \item \code{se.exposure.ls}   - A list of length \eqn{L}, each an exposure standard-error vector
#'         (used if \code{use.openGWAS = FALSE})
#'       \item \code{beta.outcome.ls}  - A list of length \eqn{L}, each an outcome effect-size vector
#'         (used if \code{use.openGWAS = FALSE})
#'       \item \code{se.outcome.ls}    - A list of length \eqn{L}, each an outcome standard-error vector
#'         (used if \code{use.openGWAS = FALSE})
#'       \item \code{step1res}         - The resulting p-value vector (length \eqn{L}) from step 1
#'     }
#'
#'   \item \code{step2} : A list relevant to step 2
#'     \itemize{
#'       \item \code{beta.exposure.mat} - A size \eqn{m^{*} \times L^{*}} matrix of exposure effect sizes,
#'         where \eqn{m^{*} = 187} and \eqn{L^{*} = 43} (used if \code{use.openGWAS = FALSE})
#'       \item \code{se.exposure.mat}   - A size \eqn{m^{*} \times L^{*}} matrix of exposure standard errors
#'         (used if \code{use.openGWAS = FALSE})
#'       \item \code{beta.outcome.vec}  - A length-\eqn{m^{*}} vector of outcome effect sizes (used if \code{use.openGWAS = FALSE})
#'       \item \code{se.outcome.vec}    - A length-\eqn{m^{*}} vector of outcome standard errors (used if \code{use.openGWAS = FALSE})
#'       \item \code{pval.exposure.mat} - A \eqn{m^{*} \times L^{*}} matrix of exposure p-values
#'         (used if \code{use.openGWAS = FALSE})
#'       \item \code{step2res}         - The result object from step 2, including \code{mvdat},
#'         \code{invalid.idx}, and \code{theta.vec}.
#'     }
#'
#'   \item \code{step3} : A list relevant to step 3
#'     \itemize{
#'       \item \code{metdrho}  - A genetic correlation matrix of size \eqn{L \times L} (249 x 249). Subsetting is needed
#'         to obtain \eqn{L^{*} \times L^{*}} matrix (a subset of the 44 x 44 matrix, which includes the outcome) used in step 3
#'       \item \code{step3res} - The final \code{susie} model output from step 3
#'     }
#' }
#'
#' @format A named list with three components: \code{step1}, \code{step2}, and \code{step3}
#'
#' @examples
#' \dontrun{
#' # Load the example data
#' data("example.dat")
#'
#' # Inspect its structure:
#' str(example.dat)
#'
#' # Examine step1 results
#' head(example.dat$step1$step1res)
#'
#' # Examine step 3 SuSiE object
#' str(example.dat$step3$step3res)
#' }
