\name{RUVadj}
\alias{RUVadj}

\title{
RUV adjust
}

\description{
Post-process and summarize the results of call to \code{\link{RUVfit}}.
}

\usage{
RUVadj(Y, fit, var.type=c("ebayes", "standard", "pooled"),
       p.type=c("standard", "rsvar", "evar"), cpginfo=NULL, ...)
}

\arguments{
\item{Y}{
The original data matrix used in the call to \code{\link{RUVfit}}.
}
\item{fit}{
A RUV model fit (a \code{list}) as returned by \code{RUVfit}.
}
\item{var.type}{Which type of estimate for sigma2 should be used 
from the call to \code{\link{variance_adjust}}? The options are "ebayes", 
"standard", or "pooled." See \code{\link{variance_adjust}} for details.
}
\item{p.type}{Which type of p-values should be used from the call to 
\code{\link{variance_adjust}}? The options are "standard", "rsvar", or "evar".
}
\item{cpginfo}{A matrix or dataframe containing information about the CpGs. 
This information is included in the summary that is returned.
}
\item{...}{Other parameters that can be passed to \code{\link{ruv}} function
\code{\link{ruv_summary}}.
}
}

\details{
This function post-processes the results of a call to \code{\link{RUVfit}} 
and then summarizes the output. The post-processing step primarily 
consists of a call to \code{\link{ruv_summary}} and \code{\link{variance_adjust}}, 
which computes various adjustments to variances, t-statistics, and and p-values. 
See \code{\link{variance_adjust}} for details. The \code{var.type} and 
\code{p.type} options determine which of these adjustments are used. 

After post-processing, the results are summarized into a list containing 
4 objects: 1) the data matrix Y; 2) a dataframe R containing information 
about the rows (samples); 3) a dataframe C containing information about 
the columns (features, e.g. genes), and 4) a list misc of other information 
returned by \code{\link{RUVfit}}.
}

\value{
An \code{list} containing:
\item{Y}{The original data matrix..
}
\item{R}{A dataframe of sample-wise information, including X, Z, and any 
other data passed in with \code{rowinfo}.
}
\item{C}{A dataframe of cpg-wise information, including p-values, 
estimated regression coefficients, estimated variances, column means, 
an index of the negative controls, and any other data passed in with 
\code{cpginfo}.
}
\item{misc}{A list of additional information returned by \code{\link{RUVfit}}.
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
lFit <- lmFit(Mval,design)
lFit2 <- eBayes(lFit)
lTop <- topTable(lFit2,coef=2,num=Inf)

# The negative control features should *not* be associated with factor of interest
# but *should* be affected by unwanted variation 
ctl <- rownames(Mval) \%in\% rownames(lTop[lTop$adj.P.Val > 0.5,])

# Perform RUV adjustment and fit
fit <- RUVfit(Y=Mval, X=group, ctl=ctl)
fit2 <- RUVadj(Y=Mval, fit=fit)

# Look at table of top results
top <- topRUV(fit2)
}
}