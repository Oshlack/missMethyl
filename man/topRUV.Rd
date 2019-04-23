\name{topRUV}
\alias{topRUV}
\title{
Table of top-ranked differentially methylated CpGs obatained from a 
differential methylation analysis using RUV
}

\description{
Extract a table of the top-ranked CpGs from a linear model fit after 
performing a differential methylation analysis using \code{RUVfit}.	
}

\usage{
topRUV(fit, number=10, p.value.cut = 1, 
cut.on = c("p.ebayes.BH","p.BH","p.rsvar.BH","p.evar.BH","p.rsvar.ebayes.BH"),
sort.by = c("p.ebayes.BH","p.BH","p.rsvar.BH","p.evar.BH","p.rsvar.ebayes.BH"))
}

\arguments{
\item{fit}{
An object containing a linear model fit produced by \code{RUVfit}, followed by
\code{RUVadj}. The fit object should be of class \code{MArrayLM}.
}

\item{number}{
integer, maximum number of genes to list. Default is 10.
}

\item{p.value.cut}{
numeric, cutoff value for adjusted p-values. Only features with lower p-values 
are listed. Must ve between 0 and 1. Default is 1.
}

\item{cut.on}{
numeric, the type of adjusted p-value that the cutoff should be applied to. 
Default is \code{p.ebayes.BH}. Other options are: \code{p.BH}, \code{p.rsvar.BH},
\code{p.evar.BH} or \code{p.rsvar.ebayes.BH}.
}

\item{sort.by}{
character string, the type of adjusted p-value that should be used for sorting. 
Default is \code{p.ebayes.BH}. Other options are: \code{p.BH}, \code{p.rsvar.BH}, 
\code{p.evar.BH} or \code{p.rsvar.ebayes.BH}. 
}
}

\details{
This function summarises the results of a differential methylation analysis 
performed using \code{RUVfit}, followed by \code{RUVadj}. The top ranked CpGs 
are selected by first ranking the adjusted p-values (Default: \code{p.ebayes.BH}), 
then ranking the raw p-values (Default: \code{p.ebayes}). 
}

\value{
Produces a dataframe with rows corresponding to the top \code{number} CpGs and 
the following columns:
\item{coefficients}{The estimated coefficients of the factor(s) of interest.
}
\item{sigma2}{Estimates of the features' variances.
}
\item{t}{t statistics for the factor(s) of interest.
}
\item{p}{P-values for the factor(s) of interest.
}
\item{multiplier}{The constant by which \code{sigma2} must be multiplied in order 
to get an estimate of the variance of \code{coefficients}
}
\item{df}{The number of residual degrees of freedom.
}
\item{W}{The estimated unwanted factors.
}
\item{alpha}{The estimated coefficients of W.
}
\item{byx}{The coefficients in a regression of Y on X (after both Y and X have 
been "adjusted" for Z).  Useful for projection plots.
}
\item{bwx}{The coefficients in a regression of W on X (after X has been 
"adjusted" for Z). Useful for projection plots.
}
\item{X}{\code{X}. Included for reference.
}
\item{k}{\code{k}. Included for reference.
}
\item{ctl}{\code{ctl}. Included for reference.
}
\item{Z}{\code{Z}. Included for reference.
}
\item{fullW0}{Can be used to speed up future calls of \code{RUVfit}.
}

The following columns may or may not be present depending on the options 
selected when \code{RUVadj} was run:
\item{p.rsvar}{P-values, after applying the method of rescaled variances.
}
\item{p.evar}{P-values, after applying the method of empirical variances.
}
\item{p.ebayes}{P-values, after applying the empirical bayes method of Smyth (2004).
}
\item{p.rsvar.ebayes}{P-values, after applying the empirical bayes method of 
Smyth (2004) and the method of rescaled variances.
}
\item{p.BH}{P-values adjusted for false discovery rate (FDR) using the method
of Benjamini and Hochberg (1995).
}
\item{p.rsvar.BH}{FDR-adjusted p-values, after applying the method of rescaled 
variances.
}
\item{p.evar.BH}{FDR-adjusted p-values, after applying the method of empirical 
variances.
}
\item{p.ebayes.BH}{FDR-adjusted p-values, after applying the empirical bayes 
method of Smyth (2004).
}
\item{p.rsvar.ebayes.BH}{FDR-adjusted p-values, after applying the empirical 
bayes method of Smyth (2004) and the method of rescaled variances.
}
}

\references{
Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: 
a practical and powerful approach to multiple testing. \emph{Journal of the 
Royal Statistical Society Series}, B, \bold{57}, 289-300.
    
Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing 
differential expression in microarray experiments. \emph{Statistical 
Applications in Genetics and Molecular Biology}, Volume 3, Article 3. 
\url{http://www.statsci.org/smyth/pubs/ebayes.pdf}.
}

\author{
Jovana Maksimovic \email{jovana.maksimovic@mcri.edu.au}
}

\seealso{
\code{\link{RUVfit}}, \code{\link{RUVadj}}, \code{\link[limma:marraylm]{MArrayLM}} 
}

\examples{
if(require(minfi) & require(minfiData) & require(limma)){

# Get methylation data for a 2 group comparison
meth <- getMeth(MsetEx)
unmeth <- getUnmeth(MsetEx)
Mval <- log2((meth + 100)/(unmeth + 100))

group<-factor(pData(MsetEx)$Sample_Group)
design<-model.matrix(~group)

# Perform initial analysis to empirically identify negative control features 
# when *not* known a priori
lFit = lmFit(Mval,design)
lFit2 = eBayes(lFit)
lTop = topTable(lFit2,coef=2,num=Inf)

# The negative control features should *not* be associated with factor of interest
# but *should* be affected by unwanted variation 
ctl = rownames(Mval) \%in\% rownames(lTop[lTop$adj.P.Val > 0.5,])

# Perform RUV adjustment and fit
fit = RUVfit(data=Mval, design=design, coef=2, ctl=ctl)
fit2 = RUVadj(fit)

# Look at table of top results
top = topRUV(fit2)
}
}

