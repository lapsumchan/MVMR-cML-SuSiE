\name{example.dat}
\alias{example.dat}
\docType{data}
\title{Bundled summary statistics and saved MVMR-cML-SuSiE results}
\format{A named list with three components:
\describe{
\item{step1}{Sample sizes, exposure/outcome effect and SE lists for 249 exposures,
and the saved univariable screening p-values.}
\item{step2}{Exposure effect, SE and p-value matrices for 187 instruments and
43 exposures, outcome effect and SE vectors, and the saved step-2 result.}
\item{step3}{A named 249 by 249 correlation matrix and the saved SuSiE fit.}
}}
\description{The packaged UK Biobank metabolite and Alzheimer disease GWAS
example permits offline execution of the three implemented steps.
The saved results are historical references, not results recomputed on loading.}
\usage{data("example.dat")}
\examples{
data("example.dat", package = "MVMRcMLSuSiE")
dim(example.dat$step2$beta.exposure.mat)
}
\keyword{datasets}
