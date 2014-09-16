\name{missMethyl-package}
\alias{missMethyl-package}
\alias{missMethyl}
\docType{package}
\title{
Introduction to the missMethyl package
}
\description{
missMethyl is a library for the analysis of Illumina's 450K human methylation BeadChip. Specifically, functions for SWAN normalisation and differential variability analysis are provided. SWAN normalisation uses probe specific information, and the differential variability procedure uses linear models which can handle any designed experiment.
}
\details{
\tabular{ll}{
Package: \tab missMethyl\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2014-06-30\cr
License: \tab GPL-2\cr
}
Normalisation of the 450K arrays can be performed using the function \code{SWAN}. 

Differential variability analysis can be performed by calling \code{varFit} followed by \code{topVar} for a list of the top ranked differentially variable CpGs between conditions.

More detailed help documentation is provided in each function's help page.
}
\author{
Belinda Phipson and Jovana Maksimovic

Maintainer: Belinda Phipson <belinda.phipson@mcri.edu.au>, Jovana Maksimovic <jovana.maksimovic@mcri.edu.au>
}
\references{
Maksimovic, J., Gordon, L., Oshlack, A. (2012).  SWAN: Subset-quantile within array normalization for illumina infinium HumanMethylation450 BeadChips. \emph{Genome Biology}, 13:R44.

Phipson, B., and Oshlack, A. (2014). DiffVar: A new method for detecting differential variability with application to methylation in cancer and aging. \emph{Genome Biology, Accepted 10 September 2014}.
}
