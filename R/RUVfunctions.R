# Extract data corrected for unwanted variation 
# Only to be used for visualisation!


#' Extract values adjusted for unwanted variation by RUVm
#' 
#' Extract values adjusted for unwanted variation by RUVm.
#' 
#' This function extracts values adjusted for unwanted variation by RUVm.
#' These values are ONLY intendeded to be used for visualisation purposes.  It
#' is NOT recommended that they are used for any further analysis.
#' 
#' @param Y A matrix of M-values.
#' @param fit The list \code{list} object produced by \code{RUVfit}.
#' @return An matrix of M-values.
#' @author Jovana Maksimovic 
#' @seealso \code{\linkS4class{MArrayLM}}
#' @examples
#' 
#' if(require(minfi) & require(minfiData) & require(limma)) {
#' 
#' # Get methylation data for a 2 group comparison
#' meth <- getMeth(MsetEx)
#' unmeth <- getUnmeth(MsetEx)
#' Mval <- log2((meth + 100)/(unmeth + 100))
#' 
#' group <- factor(pData(MsetEx)$Sample_Group, labels=c(0,1))
#' design <- model.matrix(~group)
#' 
#' # Perform initial analysis to empirically identify negative control features 
#' # when not known a priori
#' lFit <- lmFit(Mval,design)
#' lFit2 <- eBayes(lFit)
#' lTop <- topTable(lFit2,coef=2,num=Inf)
#' 
#' # The negative control features should *not* be associated with factor of 
#' # interest but *should* be affected by unwanted variation 
#' ctl <- rownames(Mval) %in% rownames(lTop[lTop$adj.P.Val > 0.5,])
#' 
#' # Perform RUV adjustment and fit
#' fit <- RUVfit(Y=Mval, X=group, ctl=ctl)
#' fit2 <- RUVadj(Y=Mval, fit=fit)
#' 
#' # get adjusted values
#' Madj <- getAdj(Y=Mval,fit=fit)
#' }
#' 
#' @export getAdj
getAdj <- function(Y, fit){
  if (is.data.frame(Y)) 
    Y <- data.matrix(Y)
  
  if (mode(Y) != "numeric") 
    stop("'Y' must be a numeric matrix.")
  
  Y <- t(Y)
  
  W <- fit$W
  alpha <- fit$alpha
  Ya <- Y - W %*% alpha
  
  t(Ya)
}

# Get Illumina negative control data


#' Extract intensity data for Illumina negative controls found on 450k or EPIC
#' arrays.
#' 
#' Extracts the intensity data for the Illumina negative controls found on 450k
#' or EPIC arrays and returns a matrix of M-values (log2 ratio of the green to
#' red intensities).
#' 
#' The getINCs function extracts the intensity data for the INCs from an
#' \linkS4class{RGChannelSet} object. The function retrieves both the green and
#' red channel intensity values and returns the data as the log2 ratio of the
#' green to red intensities. Essentially, the INCs are being treated like Type
#' II probes for which the M-values are also given as the log2 ratio of the
#' green to red intensities.
#' 
#' @param rgSet An object of class \code{RGChannelSet}.
#' @return An matrix of M-values.
#' @author Jovana Maksimovic 
#' @seealso \code{\linkS4class{RGChannelSet}}
#' @examples
#' 
#' if (require(minfi) & require(minfiData)) {
#' 
#'   INCs <- getINCs(RGsetEx)
#'   head(INCs)
#'   dim(INCs)
#' }
#' 
#' 
#' @export getINCs
getINCs <- function(rgSet){
  
  ctrls = minfi::getProbeInfo(rgSet, type = "Control")
  M.Neg = minfi::getGreen(rgSet)[ctrls$Address[ctrls$Type == "NEGATIVE"], ]
  U.Neg = minfi::getRed(rgSet)[ctrls$Address[ctrls$Type == "NEGATIVE"], ]
  
  M.Neg[M.Neg == 0] = min(M.Neg[M.Neg != 0])
  U.Neg[U.Neg == 0] = min(U.Neg[U.Neg != 0])
  
  log2(M.Neg/U.Neg)
}

# # Perform linear model fit using RUV
# RUVfit.old <- function(data, design, coef=ncol(design), ctl, 
#                    method=c("inv", "rinv", "ruv4", "ruv2"), k = NULL, ...){
#     
#   method <- match.arg(method)
#   
#   if ((method %in% c("ruv4", "ruv2")) & is.null(k))
#     stop("'k' cannot be NULL if method is 'ruv4' or 'ruv2'.")
#   
#   design <- as.matrix(design)
#   
#   if (mode(design) != "numeric") 
#     stop("'design' must be a numeric matrix.")
#   
#   if (mode(coef) != "numeric")
#     stop("'coef' must be an integer or numeric vector.")
#   
#   if (coef < 1 | coef > ncol(design))
#     stop("'coef' can only contain values >= 1 or <= ncol 'design'.")
#   
#   if (mode(ctl) != "logical")
#     stop("'ctl' must be a logical vector.")
#   
#   X <- as.matrix(design[, coef])
#   Z <- as.matrix(design[, -coef])
# 
#   Y <- as.matrix(data)
#   
#   if (mode(Y) != "numeric") 
#     stop("'data' must be a numeric matrix.")
#   
#   Y <- t(Y)
#     
#   fit <- switch(method, 
#                 inv = RUVinv(Y = Y, X = X, ctl = ctl, Z = Z, ...),
#                 rinv = RUVrinv(Y = Y, X = X, ctl = ctl, Z = Z, k = k, ...),
#                 ruv4 = RUV4(Y = Y, X = X, ctl = ctl, k = k, Z = Z, ...),
#                 ruv2 = RUV2(Y = Y, X = X, ctl = ctl, k = k, Z = Z, ...))
# 
#   return(.toMArrayLM(fit))
#   
# }

# Perform linear model fit using RUV


#' Remove unwanted variation when testing for differential methylation
#' 
#' Provides an interface similar to \code{\link{lmFit}} from
#' \code{\link{limma}} to the \code{\link{RUV2}}, \code{\link{RUV4}},
#' \code{\link{RUVinv}} and \code{\link{RUVrinv}} functions from the
#' \code{\link{ruv}} package, which facilitates the removal of unwanted
#' variation in a differential methylation analysis. A set of negative control
#' variables, as described in the references, must be specified.
#' 
#' This function depends on the \code{\link{ruv}} package and is used to
#' estimate and adjust for unwanted variation in a differential methylation
#' analysis. Briefly, the unwanted factors \code{W} are estimated using
#' negative control variables. \code{Y} is then regressed on the variables
#' \code{X}, \code{Z}, and \code{W}. For methylation data, the analysis is
#' performed on the M-values, defined as the log base 2 ratio of the methylated
#' signal to the unmethylated signal.
#' 
#' @param Y numeric \code{matrix} with rows corresponding to the features of
#' interest such as CpG sites and columns corresponding to samples or arrays.
#' @param X The factor(s) of interest. A m by p matrix, where m is the number
#' of samples and p is the number of factors of interest. Very often p = 1.
#' Factors and dataframes are also permissible, and converted to a matrix by
#' \code{\link{design.matrix}}.
#' @param ctl logical vector, \code{length == nrow(Y)}. Features that are to be
#' used as negative control variables are indicated as TRUE, all other features
#' are FALSE.
#' @param Z Any additional covariates to include in the model, typically a m by
#' q matrix.  Factors and dataframes are also permissible, and converted to a
#' matrix by \code{\link{design.matrix}}. Alternatively, may simply be 1 (the
#' default) for an intercept term. May also be NULL.
#' @param k integer, required if \code{method} is "ruv2" or "ruv4". Indicates
#' the number of unwanted factors to use. Can be 0.
#' @param method character string, indicates which \code{\link{ruv}} method
#' should be used.
#' @param ...  additional arguments that can be passed to \code{\link{RUV2}},
#' \code{\link{RUV4}}, \code{\link{RUVinv}} and \code{\link{RUVrinv}}. See
#' linked function documentation for details.
#' @return A \code{list} containing: \item{betahat}{The estimated coefficients
#' of the factor(s) of interest. A p by n matrix. } \item{sigma2}{Estimates of
#' the features' variances. A vector of length n.  } \item{t}{t statistics for
#' the factor(s) of interest. A p by n matrix. } \item{p}{P-values for the
#' factor(s) of interest. A p by n matrix. } \item{Fstats}{F statistics for
#' testing all of the factors in X simultaneously.. } \item{Fpvals}{P-values
#' for testing all of the factors in X simultaneously. } \item{multiplier}{The
#' constant by which \code{sigma2} must be multiplied in order to get an
#' estimate of the variance of \code{betahat}. } \item{df}{The number of
#' residual degrees of freedom. } \item{W}{The estimated unwanted factors. }
#' \item{alpha}{The estimated coefficients of W. } \item{byx}{The coefficients
#' in a regression of Y on X (after both Y and X have been "adjusted" for Z).
#' Useful for projection plots. } \item{bwx}{The coefficients in a regression
#' of W on X (after X has been "adjusted" for Z). Useful for projection plots.
#' } \item{X}{\code{X}. Included for reference. } \item{k}{\code{k}. Included
#' for reference. } \item{ctl}{\code{ctl}. Included for reference. }
#' \item{Z}{\code{Z}. Included for reference. } \item{fullW0}{Can be used to
#' speed up future calls of \code{RUVfit}. }
#' \item{include.intercept}{\code{include.intercept}. Included for reference. }
#' \item{method}{Character variable with value indicating which RUV method was
#' used.  Included for reference. }
#' @author Jovana Maksimovic 
#' @seealso \code{\link{RUV2}}, \code{\link{RUV4}}, \code{\link{RUVinv}},
#' \code{\link{RUVrinv}}, \code{\link{topRUV}}
#' @references Gagnon-Bartsch JA, Speed TP. (2012). Using control genes to
#' correct for unwanted variation in microarray data. \emph{Biostatistics}.
#' \bold{13}(3), 539-52.  Available at:
#' \url{http://biostatistics.oxfordjournals.org/content/13/3/539.full}.
#' 
#' Gagnon-Bartsch, Jacob, and Speed. 2013. Removing Unwanted Variation from
#' High Dimensional Data with Negative Controls. Available at:
#' \url{http://statistics.berkeley.edu/tech-reports/820}.
#' @examples
#' 
#' if(require(minfi) & require(minfiData) & require(limma)) {
#' # Get methylation data for a 2 group comparison
#' meth <- getMeth(MsetEx)
#' unmeth <- getUnmeth(MsetEx)
#' Mval <- log2((meth + 100)/(unmeth + 100))
#' group <- factor(pData(MsetEx)$Sample_Group)
#' design <- model.matrix(~group)
#' # Perform initial analysis to empirically identify negative control features 
#' # when not known a priori
#' lFit <- lmFit(Mval,design)
#' lFit2 <- eBayes(lFit)
#' lTop <- topTable(lFit2,coef=2,num=Inf)
#' # The negative control features should *not* be associated with factor of 
#' # interest but *should* be affected by unwanted variation 
#' ctl <- rownames(Mval) %in% rownames(lTop[lTop$adj.P.Val > 0.5,])
#' # Perform RUV adjustment and fit
#' fit <- RUVfit(Y=Mval, X=group, ctl=ctl)
#' fit2 <- RUVadj(Y=Mval, fit=fit)
#' # Look at table of top results
#' top <- topRUV(fit2)
#' }
#' 
#' @export RUVfit
RUVfit <- function(Y, X, ctl, Z = 1, k = NULL,
                   method=c("inv", "rinv", "ruv4", "ruv2"), ...){
  
  method <- match.arg(method)
  
  if ((method %in% c("ruv4", "ruv2")) & is.null(k))
    stop("'k' cannot be NULL if method is 'ruv4' or 'ruv2'.")
  
  if (mode(ctl) != "logical")
    stop("'ctl' must be a logical vector.")
  
  if (is.data.frame(Y)) 
    Y <- data.matrix(Y)
  
  if (mode(Y) != "numeric") 
    stop("'Y' must be a numeric matrix.")
  
  Y <- t(Y)
  
  fit <- switch(method, 
                inv = ruv::RUVinv(Y = Y, X = X, ctl = ctl, Z = Z, ...),
                rinv = ruv::RUVrinv(Y = Y, X = X, ctl = ctl, Z = Z, k = k, ...),
                ruv4 = ruv::RUV4(Y = Y, X = X, ctl = ctl, k = k, Z = Z, ...),
                ruv2 = ruv::RUV2(Y = Y, X = X, ctl = ctl, k = k, Z = Z, ...))
  
  return(fit)
  
}

# Calculate rescaled variances, empirical variances, etc. 
# For use with RUV model fits.


#' RUV adjust
#' 
#' Post-process and summarize the results of call to \code{\link{RUVfit}}.
#' 
#' This function post-processes the results of a call to \code{\link{RUVfit}}
#' and then summarizes the output. The post-processing step primarily consists
#' of a call to \code{\link{ruv_summary}} and \code{\link{variance_adjust}},
#' which computes various adjustments to variances, t-statistics, and and
#' p-values.  See \code{\link{variance_adjust}} for details. The
#' \code{var.type} and \code{p.type} options determine which of these
#' adjustments are used.
#' 
#' After post-processing, the results are summarized into a list containing 4
#' objects: 1) the data matrix Y; 2) a dataframe R containing information about
#' the rows (samples); 3) a dataframe C containing information about the
#' columns (features, e.g. genes), and 4) a list misc of other information
#' returned by \code{\link{RUVfit}}.
#' 
#' @param Y The original data matrix used in the call to \code{\link{RUVfit}}.
#' @param fit A RUV model fit (a \code{list}) as returned by \code{RUVfit}.
#' @param var.type Which type of estimate for sigma2 should be used from the
#' call to \code{\link{variance_adjust}}? The options are "ebayes", "standard",
#' or "pooled." See \code{\link{variance_adjust}} for details.
#' @param p.type Which type of p-values should be used from the call to
#' \code{\link{variance_adjust}}? The options are "standard", "rsvar", or
#' "evar".
#' @param cpginfo A matrix or dataframe containing information about the CpGs.
#' This information is included in the summary that is returned.
#' @param ... Other parameters that can be passed to \code{\link{ruv}} function
#' \code{\link{ruv_summary}}.
#' @return An \code{list} containing: \item{Y}{The original data matrix.. }
#' \item{R}{A dataframe of sample-wise information, including X, Z, and any
#' other data passed in with \code{rowinfo}. } \item{C}{A dataframe of cpg-wise
#' information, including p-values, estimated regression coefficients,
#' estimated variances, column means, an index of the negative controls, and
#' any other data passed in with \code{cpginfo}. } \item{misc}{A list of
#' additional information returned by \code{\link{RUVfit}}. }
#' @author Jovana Maksimovic \email{jovana.maksimovic@@mcri.edu.au}
#' @seealso \code{\linkS4class{MArrayLM}}, \code{\link{RUV2}},
#' \code{\link{RUV4}}, \code{\link{RUVinv}}, \code{\link{RUVrinv}},
#' \code{\link{p.adjust}}, \code{\link{get_empirical_variances}},
#' \code{\link{sigmashrink}}
#' @references Benjamini, Y., and Hochberg, Y. (1995). Controlling the false
#' discovery rate: a practical and powerful approach to multiple testing.
#' \emph{Journal of the Royal Statistical Society Series}, B, \bold{57},
#' 289-300.
#' 
#' Gagnon-Bartsch JA, Speed TP. (2012). Using control genes to correct for
#' unwanted variation in microarray data. \emph{Biostatistics}. \bold{13}(3),
#' 539-52.  Available at:
#' \url{http://biostatistics.oxfordjournals.org/content/13/3/539.full}.
#' 
#' Gagnon-Bartsch, Jacob, and Speed. 2013. Removing Unwanted Variation from
#' High Dimensional Data with Negative Controls. Available at:
#' \url{http://statistics.berkeley.edu/tech-reports/820}.
#' 
#' Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing
#' differential expression in microarray experiments. \emph{Statistical
#' Applications in Genetics and Molecular Biology}, Volume 3, Article 3.
#' \url{http://www.statsci.org/smyth/pubs/ebayes.pdf}.
#' @examples
#' 
#' if(require(minfi) & require(minfiData) & require(limma)) {
#' 
#' # Get methylation data for a 2 group comparison
#' meth <- getMeth(MsetEx)
#' unmeth <- getUnmeth(MsetEx)
#' Mval <- log2((meth + 100)/(unmeth + 100))
#' 
#' group<-factor(pData(MsetEx)$Sample_Group)
#' design<-model.matrix(~group)
#' 
#' # Perform initial analysis to empirically identify negative control features 
#' # when not known a priori
#' lFit <- lmFit(Mval,design)
#' lFit2 <- eBayes(lFit)
#' lTop <- topTable(lFit2,coef=2,num=Inf)
#' 
#' # The negative control features should *not* be associated with factor of 
#' # interest but *should* be affected by unwanted variation 
#' ctl <- rownames(Mval) %in% rownames(lTop[lTop$adj.P.Val > 0.5,])
#' 
#' # Perform RUV adjustment and fit
#' fit <- RUVfit(Y=Mval, X=group, ctl=ctl)
#' fit2 <- RUVadj(Y=Mval, fit=fit)
#' 
#' # Look at table of top results
#' top <- topRUV(fit2)
#' }
#' 
#' @export RUVadj
RUVadj <- function(Y, fit, var.type=c("ebayes", "standard", "pooled"),
                   p.type=c("standard", "rsvar", "evar"), cpginfo=NULL, ...){
  
  if (is.data.frame(Y)) 
    Y <- data.matrix(Y)
  
  if (mode(Y) != "numeric") 
    stop("'Y' must be a numeric matrix.")
  
  Y <- t(Y)
  
  var.type <- match.arg(var.type)
  p.type <- match.arg(p.type)
  
  fitsum <- ruv::ruv_summary(Y, fit, colinfo=cpginfo, var.type=var.type, 
                        p.type=p.type, ...)
  
  return(fitsum)
}

# .toMArrayLM <- function(fit){
#     
#     obj <- new("MArrayLM")
#     
#     obj$coefficients <- t(fit$betahat)
#     obj$sigma2 <- fit$sigma2
#     obj$tvals <- fit$tvals
#     obj$pvals <- fit$pvals 
#     obj$multiplier <- fit$multiplier
#     obj$df <- fit$df 
#     obj$W <- t(fit$W) 
#     obj$alpha <- t(fit$alpha) 
#     obj$byx <- t(fit$byx) 
#     obj$bwx <- t(fit$bwx) 
#     obj$X <- fit$X 
#     obj$k <- fit$k
#     obj$ctl <- fit$ctl
#     obj$Z <- fit$Z 
#     obj$fullW0 <- t(fit$fullW0) 
#     obj$lambda <- fit$lambda
#     obj$t <- t(fit$t)
#     obj$p <- t(fit$p)
#     
#     slots <- attributes(fit)$names
#     
#     if ("p.rsvar" %in% slots) 
#         obj$p.rsvar <- t(fit$p.rsvar)
#     
#     if ("p.evar" %in% slots) 
#         obj$p.evar <- t(fit$p.evar)
#     
#     if ("p.ebayes" %in% slots) 
#         obj$p.ebayes <- t(fit$p.ebayes)
#     
#     if ("p.rsvar.ebayes" %in% slots) 
#         obj$p.rsvar.ebayes <- t(fit$p.rsvar.ebayes)
#     
#     if ("p.BH" %in% slots) 
#         obj$p.BH <- t(fit$p.BH)
#     
#     if ("p.rsvar.BH" %in% slots) 
#         obj$p.rsvar.BH <- t(fit$p.rsvar.BH)
#     
#     if ("p.evar.BH" %in% slots) 
#         obj$p.evar.BH <- t(fit$p.evar.BH)
#     
#     if ("p.ebayes.BH" %in% slots) 
#         obj$p.ebayes.BH <- t(fit$p.ebayes.BH)
#     
#     if ("p.rsvar.ebayes.BH" %in% slots) 
#         obj$p.rsvar.ebayes.BH <- t(fit$p.rsvar.ebayes.BH)
#     
#     obj
# }

# .toList <- function(fit){
# 
#     obj <- list()
#     
#     obj$betahat <- t(fit$coefficients)
#     obj$sigma2 <- fit$sigma2
#     obj$tvals <- fit$tvals
#     obj$pvals <- fit$pvals 
#     obj$multiplier <- fit$multiplier
#     obj$df <- fit$df 
#     obj$W <- t(fit$W) 
#     obj$alpha <- t(fit$alpha) 
#     obj$byx <- t(fit$byx) 
#     obj$bwx <- t(fit$bwx) 
#     obj$X <- fit$X 
#     obj$k <- fit$k
#     obj$ctl <- fit$ctl
#     obj$Z <- fit$Z 
#     obj$fullW0 <- t(fit$fullW0) 
#     obj$lambda <- fit$lambda
#     obj$t <- t(fit$t)
#     obj$p <- t(fit$p)
#     
#     obj
# }

# Calculate rescaled variances, empirical variances, etc. 
# For use with RUV model fits.
# RUVadj.old <- function(fit, ebayes = TRUE, evar = FALSE, rsvar = FALSE, ...){
#     
#     fit <- variance_adjust(fit = .toList(fit), ebayes = ebayes, evar = evar, 
#                            rsvar = rsvar, ...)
#     
#     return(.toMArrayLM(fit))
# }



#' Table of top-ranked differentially methylated CpGs obatained from a
#' differential methylation analysis using RUV
#' 
#' Extract a table of the top-ranked CpGs from a linear model fit after
#' performing a differential methylation analysis using \code{RUVfit} and
#' \code{RUVadj}.
#' 
#' This function summarises the results of a differential methylation analysis
#' performed using \code{RUVfit}, followed by \code{RUVadj}. The top ranked
#' CpGs are sorted by p-value.
#' 
#' @param fitsum An object containing the summary fit object produced by
#' \code{RUVadj}.  The object should be a \code{list}.
#' @param number integer, maximum number of genes to list. Default is 10.
#' @param sort.by character string, what the results should be sorted by.
#' Default is unadjusted p-value.
#' @param p.BH numeric, cutoff value for Benjamini-Hochberg adjusted p-values.
#' Only features with lower p-values are listed. Must be between 0 and 1.
#' Default is 1.
#' @return Produces a dataframe with rows corresponding to the top
#' \code{number} CpGs and the following columns: F.p F.p.BH p_X1 p.BH_X1 b_X1
#' sigma2 var.b_X1 fit.ctl mean
#' 
#' \item{F.p}{P-values for testing all of the factors of interest
#' simultaneously. } \item{F.p.BH}{Benjamini-Hochberg adjusted p-values for
#' testing all of the factors of interest simultaneously. }
#' \item{p_X1}{p-values for the factor of interest. }
#' \item{p.BH_X1}{Benjamini-Hochberg adjusted p-values for the factor of
#' interest. } \item{b_X1}{The estimated coefficients of the factor of
#' interest. } \item{sigma2}{Estimate of the methylation variance. }
#' \item{var.b_X1}{Variance estimate of \code{betahat}. }
#' \item{fit.ctl}{logical, indicating whether CpG was designated as a negative
#' control.  } \item{mean}{The mean methylation (M-value). }
#' @author Jovana Maksimovic \email{jovana.maksimovic@@mcri.edu.au}
#' @seealso \code{\link{RUVfit}}, \code{\link{RUVadj}},
#' \code{\link[limma:marraylm]{MArrayLM}}
#' @references Benjamini, Y., and Hochberg, Y. (1995). Controlling the false
#' discovery rate: a practical and powerful approach to multiple testing.
#' \emph{Journal of the Royal Statistical Society Series}, B, \bold{57},
#' 289-300.
#' 
#' Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing
#' differential expression in microarray experiments. \emph{Statistical
#' Applications in Genetics and Molecular Biology}, Volume 3, Article 3.
#' \url{http://www.statsci.org/smyth/pubs/ebayes.pdf}.
#' @examples
#' 
#' if(require(minfi) & require(minfiData) & require(limma)){
#' 
#' # Get methylation data for a 2 group comparison
#' meth <- getMeth(MsetEx)
#' unmeth <- getUnmeth(MsetEx)
#' Mval <- log2((meth + 100)/(unmeth + 100))
#' 
#' group <- factor(pData(MsetEx)$Sample_Group)
#' design <- model.matrix(~group)
#' 
#' # Perform initial analysis to empirically identify negative control features 
#' # when *not* known a priori
#' lFit <- lmFit(Mval,design)
#' lFit2 <- eBayes(lFit)
#' lTop <- topTable(lFit2,coef=2,num=Inf)
#' 
#' # The negative control features should *not* be associated with factor of 
#' # interest but *should* be affected by unwanted variation 
#' ctl <- rownames(Mval) %in% rownames(lTop[lTop$adj.P.Val > 0.5,])
#' 
#' # Perform RUV adjustment and fit
#' fit <- RUVfit(Y=Mval, X=group, ctl=ctl)
#' fit2 <- RUVadj(Y=Mval, fit=fit)
#' 
#' # Look at table of top results
#' top <- topRUV(fit2)
#' }
#' 
#' @export topRUV
topRUV <- function (fitsum, number = 10, sort.by = c("p","F.p"), p.BH = 1){
  
  tab <- fitsum$C
  
  if (p.BH < 1) {
    
    sig <- (tab[,grepl("p.BH_", colnames(tab))] < p.BH)
    
    if (any(is.na(sig)))
      sig[is.na(sig)] <- FALSE
    
    if (all(!sig))
      return(data.frame())
    
    tab <- tab[sig,]
  }
  
  sort.by <- match.arg(sort.by)
  
  ord <- switch(sort.by, p = order(tab[,grepl("p_", colnames(tab))], 
                                   decreasing = FALSE), 
                F.p = order(tab$F.p, decreasing=FALSE))
  
  if (nrow(tab) < number)
    number <- nrow(tab)
  
  if (number < 1)
    return(data.frame())
  
  top <- ord[1:number]
  
  tab[top,]
}

# topRUV.old <- function (fit, number = 10, p.value.cut = 1,
#                          cut.on = c("p.ebayes.BH","p.BH","p.rsvar.BH",
#                                     "p.evar.BH","p.rsvar.ebayes.BH"),
#                          sort.by = c("p.ebayes.BH","p.BH","p.rsvar.BH",
#                                      "p.evar.BH","p.rsvar.ebayes.BH")){
#   
#   coefficients <- fit$coefficients
#   t <- fit$t
#   p <- fit$p
#   ID <- rownames(coefficients)
# 
#   p.rsvar <- NULL
#   p.evar <- NULL
#   p.ebayes <- NULL
#   p.rsvar.ebayes <- NULL
#   p.BH <- NULL
#   p.rsvar.BH <- NULL
#   p.evar.BH <- NULL
#   p.ebayes.BH <- NULL
#   p.rsvar.ebayes.BH <- NULL
#   
#   slots <- attributes(fit)$names
#   
#   if ("p.rsvar" %in% slots) 
#     p.rsvar <- fit$p.rsvar
#   
#   if ("p.evar" %in% slots) 
#     p.evar <- fit$p.evar
#   
#   if ("p.ebayes" %in% slots) 
#     p.ebayes <- fit$p.ebayes
#   
#   if ("p.rsvar.ebayes" %in% slots) 
#     p.rsvar.ebayes <- fit$p.rsvar.ebayes
#   
#   if ("p.BH" %in% slots) 
#     p.BH <- fit$p.BH
#   
#   if ("p.rsvar.BH" %in% slots) 
#     p.rsvar.BH <- fit$p.rsvar.BH
#   
#   if ("p.evar.BH" %in% slots) 
#     p.evar.BH <- fit$p.evar.BH
#   
#   if ("p.ebayes.BH" %in% slots) 
#     p.ebayes.BH <- fit$p.ebayes.BH
#   
#   if ("p.rsvar.ebayes.BH" %in% slots) 
#     p.rsvar.ebayes.BH <- fit$p.rsvar.ebayes.BH
#   
#   tab <- data.frame(row.names = ID, coefficients = coefficients, t = t, p = p, 
#                     stringsAsFactors = FALSE, check.names=FALSE)
#   
#   if (!is.null(p.BH)) 
#     tab <- data.frame(tab, p.BH = p.BH, stringsAsFactors = FALSE, 
#                       check.names=FALSE)
#   
#   if (!is.null(p.rsvar)) 
#     tab <- data.frame(tab, p.rsvar = p.rsvar, stringsAsFactors = FALSE, 
#                       check.names=FALSE)
#   
#   if (!is.null(p.rsvar.BH)) 
#     tab <- data.frame(tab, p.rsvar.BH = p.rsvar.BH, stringsAsFactors = FALSE, 
#                       check.names=FALSE)
#   
#   if (!is.null(p.evar)) 
#     tab <- data.frame(tab, p.evar = p.evar, stringsAsFactors = FALSE, 
#                       check.names=FALSE)
#   
#   if (!is.null(p.evar.BH)) 
#     tab <- data.frame(tab, p.evar.BH = p.evar.BH, stringsAsFactors = FALSE, 
#                       check.names=FALSE)
#   
#   if (!is.null(p.rsvar.ebayes)) 
#     tab <- data.frame(tab, p.rsvar.ebayes = p.rsvar.ebayes, 
#                       stringsAsFactors = FALSE, check.names=FALSE)
#   
#   if (!is.null(p.rsvar.ebayes.BH)) 
#     tab <- data.frame(tab, p.rsvar.ebayes.BH = p.rsvar.ebayes.BH, 
#                       stringsAsFactors = FALSE, check.names=FALSE)
#   
#   if (!is.null(p.ebayes)) 
#     tab <- data.frame(tab, p.ebayes = p.ebayes, stringsAsFactors = FALSE, 
#                       check.names=FALSE)
#   
#   if (!is.null(p.ebayes.BH)) 
#     tab <- data.frame(tab, p.ebayes.BH = p.ebayes.BH, stringsAsFactors = FALSE, 
#                       check.names=FALSE)
#   
#   if (p.value.cut < 1) {
#     
#     cut.on <- match.arg(cut.on)
#     
#     if (!(cut.on %in% slots))
#       cat(sprintf("Cannot threshold on '%s' because it is not in the 'fit' object.", 
#                   cut.on))
#     
#     adj.p.value <- get(cut.on, fit)
#     
#     sig <- (adj.p.value < p.value.cut)
#     
#     if (any(is.na(sig)))
#       sig[is.na(sig)] <- FALSE
#     
#     if (all(!sig))
#       return(data.frame())
#     
#     tab <- tab[sig,]
#   }
#                  
#   sort.by <- match.arg(sort.by)
#   
#   if (!(sort.by %in% slots))
#     stop(sprintf("Cannot sort by '%s' because it is not in the 'fit' object.", 
#                  sort.by))
#                  
#   ord <- switch(sort.by, p.ebayes.BH = order(tab$p.ebayes.BH, tab$p.ebayes, 
#                                              decreasing = FALSE), 
#                 p.BH = order(tab$p.BH, tab$p, decreasing=FALSE), 
#                 p.rsvar.BH = order(tab$p.rsvar.BH, tab$p.rsvar, 
#                                    decreasing=FALSE), 
#                 p.evar.BH = order(tab$p.evar.BH, tab$p.evar, 
#                                   decreasing=FALSE), 
#                 p.rsvar.ebayes.BH = order(tab$p.rsvar.ebayes.BH, 
#                                           tab$p.rsvar.ebayes, decreasing=FALSE))
#   
#   if (nrow(tab) < number)
#     number <- nrow(tab)
#                  
#   if (number < 1)
#     return(data.frame())
#   
#   top <- ord[1:number]
#   
#   tab[top,]
# }
