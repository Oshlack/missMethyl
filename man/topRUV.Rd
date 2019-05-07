\name{topRUV}
\alias{topRUV}
\title{
Table of top-ranked differentially methylated CpGs obatained from a 
differential methylation analysis using RUV
}

\description{
Extract a table of the top-ranked CpGs from a linear model fit after 
performing a differential methylation analysis using \code{RUVfit} 
and \code{RUVadj}.	
}

\usage{
topRUV(fitsum, number = 10, sort.by = c("p","F.p"), p.BH = 1)
}

\arguments{
\item{fitsum}{
An object containing the summary fit object produced by \code{RUVadj}. 
The object should be a \code{list}.
}

\item{number}{
integer, maximum number of genes to list. Default is 10.
}

\item{sort.by}{
character string, what the results should be sorted by. Default is unadjusted p-value. 
}

\item{p.BH}{
numeric, cutoff value for Benjamini-Hochberg adjusted p-values. Only features with 
lower p-values are listed. Must be between 0 and 1. Default is 1.
}

}

\details{
This function summarises the results of a differential methylation analysis 
performed using \code{RUVfit}, followed by \code{RUVadj}. The top ranked CpGs 
are sorted by p-value. 
}

\value{
Produces a dataframe with rows corresponding to the top \code{number} CpGs and 
the following columns:
F.p	F.p.BH	p_X1	p.BH_X1	b_X1	sigma2	var.b_X1	fit.ctl	mean

\item{F.p}{P-values for testing all of the factors of interest simultaneously.
}
\item{F.p.BH}{Benjamini-Hochberg adjusted p-values for testing all of the 
factors of interest simultaneously.
}
\item{p_X1}{p-values for the factor of interest.
}
\item{p.BH_X1}{Benjamini-Hochberg adjusted p-values for the factor of interest.
}
\item{b_X1}{The estimated coefficients of the factor of interest.
}
\item{sigma2}{Estimate of the methylation variance.
}
\item{var.b_X1}{Variance estimate of \code{betahat}.
}
\item{fit.ctl}{logical, indicating whether CpG was designated as a negative control. 
}
\item{mean}{The mean methylation (M-value).
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

group <- factor(pData(MsetEx)$Sample_Group)
design <- model.matrix(~group)

# Perform initial analysis to empirically identify negative control features 
# when *not* known a priori
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