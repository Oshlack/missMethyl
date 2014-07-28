\name{varFit}
\alias{varFit}
\alias{varFit.default}
\alias{varFit.MethylSet}
\title{
Testing for differential variability
}
\description{
Fit linear model on mean absolute deviations for each CpG given a series of methylation arrays
}
\usage{
varFit(data, design = NULL, coef = NULL, type = NULL, trend = TRUE, robust = TRUE)
}
\arguments{
  \item{data}{
	object of class \code{MethylSet} or \code{matrix} with rows corresponding to the features of interest such as CpG sites and columns corresponding to samples or arrays
}
  \item{design}{
	the design matrix of the experiment, with rows corresponding to arrays/samples and columns to coefficients to be estimated. Defaults to the unit vector.
}
  \item{coef}{
	scalar indicating which column in the design matrix corresponds to the comparison of interest to test for differential variability. Defaults to last column in the design matrix.
}
  \item{type}{
	character string, \code{"AD"} for absolute residuals or \code{"SQ"} for squared residuals. Default is absolute.
}
  \item{trend}{
	logical, if true fits a mean variance trend on the absolute or squared deviations
}
  \item{robust}{
	logical, if true performs robust empirical Bayes shrinkage of the variances for the moderated t statistics
}
}
\details{
This function depends on the \code{limma} package and is used to rank features such as CpG sites or genes in order of evidence of differential variability between two groups. Additional covariates can be taken into account with the design matrix. The analysis is performed on the M-values, defined as the log base 2 ratio of the methylated signal to the unmethylated signal. If a \code{MethylSet} object is supplied, M-values are extracted with an offset of 100 added to the numerator and denominator. A measure of variability is calculated for each CpG in each sample by subtracting out the group mean and taking the absolute or squared deviation. A linear model is then fitted to the absolute or squared deviations. The residuals of the linear model fit are subjected to empirical Bayes shrinkage and moderated t statistics (Smyth, 2004) calculated. False discovery rates are calculated using the method of Benjamini and Hochberg (1995).

Since the output is of class \code{MArrayLM}, any functions that can be applied to fit objects from \code{lmFit} and \code{eBayes} can be applied, for example, \code{topTable} and \code{decideTests}. 
}
\value{
	produces an object of class \code{MArrayLM} (see \code{\link{MArrayLM-class}}) containing everything found in a fitted model object produced by \code{lmFit} and \code{eBayes}
}
\references{
Phipson, B., and Oshlack, A. (2013). A method for detecting differential variability in methylation data shows CpG islands are highly variably methylated in cancers. \emph{Manuscript submitted}.

Smyth, G.K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. \emph{Statistical Applications in Genetics and Molecular Biology}, Volume \bold{3}, Article 3.

Smyth, G. K. (2005). Limma: linear models for microarray data. In: \emph{Bioinformatics and Computational Biology Solutions using R and Bioconductor}. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, 2005. 

Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. \emph{Journal of the Royal Statistical Society Series}, B, \bold{57}, 289-300.
}
\seealso{
\code{\link{topVar}}, \code{\link{getLeveneResiduals}}, \code{\link{lmFit}}, \code{\link{eBayes}}, \code{\link{topTable}}, \code{\link{decideTests}}
}
\author{
	Belinda Phipson
}
\examples{
# Randomly generate data for a 2 group problem with 100 CpG sites and 5 arrays in each group. 

y<-matrix(rnorm(1000),ncol=10)

group<-factor(rep(c(1,2),each=5))
design<-model.matrix(~group)

# Fit linear model for differential variability
vfit<-varFit(y,design,coef=2)

# Look at top table of results

topVar(vfit,coef=2)
}

