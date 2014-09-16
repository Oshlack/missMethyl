\name{getLeveneResiduals}
\alias{getLeveneResiduals}
\title{
Obtain Levene residuals
}
\description{
Obtain absolute or squared Levene residuals for each CpG given a series of methylation arrays
}
\usage{
getLeveneResiduals(data, design = NULL, type = NULL)
}
\arguments{
  \item{data}{
	object of class \code{matrix} of M values, with rows corresponding to features of interest such as CpG sites and columns corresponding to samples or arrays
}
  \item{design}{
	the design matrix of the experiment, with rows corresponding to arrays/samples and columns to coefficients to be estimated. Defaults to the unit vector.
}
  \item{type}{
	character string, \code{"AD"} for absolute residuals or \code{"SQ"} for squared residuals. Default is \code{"AD"}.
}
}
\details{
	This function will return absolute or squared Levene residuals given a matrix of M values and a design matrix. This can be used for graphing purposes or for downstream analysis such a gene set testing based on differential variability rather than differential methylation. If no design matrix is given, the residuals are determined by treating all samples as coming from one group. 
}
\value{
	Returns a list with three components. \code{data} contains a matrix of absolute or squared residuals, \code{AvgVar} is a vector of sample variances and \code{LogVarRatio} corresponds to the columns of the design matrix and is usually the ratios of the log of the group variances.
}
\references{
Phipson, B., and Oshlack, A. (2014). A method for detecting differential variability in methylation data shows CpG islands are highly variably methylated in cancers. \emph{Genome Biology, Accepted 10 September 2014}.
}
\author{
	Belinda Phipson
}

\seealso{
\code{\link{varFit}}
}
\examples{
# Randomly generate data for a 2 group problem with 100 CpG sites and 5 arrays in each group
y <- matrix(rnorm(1000),ncol=10)

group <- factor(rep(c(1,2),each=5))
design <- model.matrix(~group)

# Get absolute Levene Residuals
resid <- getLeveneResiduals(y,design)

# Plot the first CpG
barplot(resid$data[1,],col=rep(c(2,4),each=5),ylab="Absolute Levene Residuals",names=group)
}
