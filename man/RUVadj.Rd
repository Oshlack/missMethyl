\name{RUVadj}
\alias{RUVadj}

\title{
Adjust estimated variances
}

\description{
Calculate rescaled variances, empirical variances, etc. For use with RUV model 
fits produced using \code{RUVfit}.
}

\usage{
RUVadj(fit, ebayes = TRUE, evar = FALSE, rsvar = FALSE, ...)
}

\arguments{
\item{fit}{
An object of class \code{MArrayLM} as returned by \code{RUVfit}.
}
\item{ebayes}{
A logical variable.  Should empirical bayes estimates be calculated?
}
\item{evar}{
A logical variable.  Should empirical variance estimates be calculated?
}
\item{rsvar}{
A logical variable.  Should rescaled variance estimates be calculated?
}
\item{...}{
additional arguments that can be passed to \code{\link{variance_adjust}}. See 
linked function documentation for details.
}
}

\details{
Adjust variance. By default only the empirical bayes method of Smyth (2004) is 
performed.
}

\value{
An object of class \code{MArrayLM} containing:
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

The following items may or may not be present depending on the options 
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

Gagnon-Bartsch JA, Speed TP. (2012). Using control genes to correct for unwanted 
variation in microarray data. \emph{Biostatistics}. \bold{13}(3), 539-52. 
Available at: \url{http://biostatistics.oxfordjournals.org/content/13/3/539.full}.

Gagnon-Bartsch, Jacob, and Speed. {2013}. Removing Unwanted Variation from High 
Dimensional Data with Negative Controls. Available at: 
\url{http://statistics.berkeley.edu/tech-reports/820}.

Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing 
differential expression in microarray experiments. \emph{Statistical 
Applications in Genetics and Molecular Biology}, Volume 3, Article 3. 
\url{http://www.statsci.org/smyth/pubs/ebayes.pdf}.
}

\author{
Jovana Maksimovic \email{jovana.maksimovic@mcri.edu.au}
}

\seealso{
\code{\linkS4class{MArrayLM}}, \code{\link{RUV2}}, \code{\link{RUV4}}, 
\code{\link{RUVinv}}, \code{\link{RUVrinv}}, \code{\link{p.adjust}}, 
\code{\link{get_empirical_variances}}, \code{\link{sigmashrink}}
}

\examples{
if(require(minfi) & require(minfiData) & require(limma)) {

# Get methylation data for a 2 group comparison
meth <- getMeth(MsetEx)
unmeth <- getUnmeth(MsetEx)
Mval <- log2((meth + 100)/(unmeth + 100))

group<-factor(pData(MsetEx)$Sample_Group)
design<-model.matrix(~group)

# Perform initial analysis to empirically identify negative control features 
# when not known a priori
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

