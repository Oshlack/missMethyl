#' Testing for differential variability
#' 
#' Fit linear model on mean absolute or squared deviations for each CpG given a
#' series of methylation arrays
#' 
#' This function depends on the \code{limma} package and is used to rank
#' features such as CpG sites or genes in order of evidence of differential
#' variability between different comparisons corresponding to the columns of
#' the design matrix.  A measure of variability is calculated for each CpG in
#' each sample by subtracting out the group mean and taking the absolute or
#' squared deviation. A linear model is then fitted to the absolute or squared
#' deviations. The residuals of the linear model fit are subjected to empirical
#' Bayes shrinkage and moderated t statistics (Smyth, 2004) calculated. False
#' discovery rates are calculated using the method of Benjamini and Hochberg
#' (1995).
#' 
#' Please always specify the \code{coef} parameter in the call to \code{varFit}, 
#' which indicates which groups are to be tested for differential variability.
#' If \code{coef} is not specified, then group means are estimated based on all
#' the columns of the design matrix and subtracted out before testing for
#' differential variability. If the design matrix contains nuisance parameters,
#' then subsetting the design matrix columns by \code{coef} should remove these
#' columns from the design matrix. If the design matrix includes an intercept
#' term, this should be included in \code{coef}. The nuisance parameters are
#' included in the linear model fit to the absolute or squared deviations, but
#' should not be considered when subtracting group means to obtain the
#' deviations. Note that design matrices without an intercept term are
#' permitted, and specific contrasts tested using the function
#' \code{contrasts.varFit}.
#' 
#' For methylation data, the analysis is performed on the M-values, defined as
#' the log base 2 ratio of the methylated signal to the unmethylated signal. If
#' a \code{MethylSet} object is supplied, M-values are extracted with an offset
#' of 100 added to the numerator and denominator.
#' 
#' For testing differential variability on RNA-Seq data, a \code{DGEList}
#' object can be supplied directly to the function. A \code{voom}
#' transformation is applied before testing for differential variability. The
#' weights calculated in \code{voom} are used in the linear model fit.
#' 
#' Since the output is of class \code{MArrayLM}, any functions that can be
#' applied to fit objects from \code{lmFit} and \code{eBayes} can be applied,
#' for example, \code{topTable} and \code{decideTests}.
#' 
#' @aliases varFit varFit.default varFit.DGEList varFit.MethylSet
#' @param data Object of class \code{MethylSet} or \code{matrix} of M-values 
#' with rows corresponding to the features of interest such as CpG sites and 
#' columns corresponding to samples or arrays.
#' @param design The design matrix of the experiment, with rows corresponding
#' to arrays/samples and columns to coefficients to be estimated. Defaults to
#' the unit vector.
#' @param coef The columns of the design matrix containing the comparisons to
#' test for differential variability. Defaults to all columns of design matrix.
#' @param type Character string, \code{"AD"} for absolute residuals or
#' \code{"SQ"} for squared residuals. Default is absolute.
#' @param trend Logical, if true fits a mean variance trend on the absolute or
#' squared deviations.
#' @param robust Logical, if true performs robust empirical Bayes shrinkage of
#' the variances for the moderated t statistics.
#' @param weights Non-negative observation weights. Can be a numeric matrix of
#' individual weights, of same size as the object matrix, or a numeric vector
#' of array weights, or a numeric vector of gene/feature weights.
#' @return Produces an object of class \code{MArrayLM} (see
#' \code{\link{MArrayLM-class}}) containing everything found in a fitted model
#' object produced by \code{lmFit} and \code{eBayes} as well as a vector
#' containing the sample CpG-wise variances and a matrix of LogVarRatios
#' corresponding to the differential variability analysis.
#' @author Belinda Phipson
#' @seealso \code{\link{contrasts.varFit}}, \code{\link{topVar}},
#' \code{\link{getLeveneResiduals}}, \code{\link{lmFit}}, \code{\link{eBayes}},
#' \code{\link{topTable}}, \code{\link{decideTests}}, \code{\link{voom}}
#' @references Phipson, B., and Oshlack, A. (2014). A method for detecting
#' differential variability in methylation data shows CpG islands are highly
#' variably methylated in cancers. \emph{Genome Biology}, \bold{15}:465.
#' 
#' Smyth, G.K. (2004). Linear models and empirical Bayes methods for assessing
#' differential expression in microarray experiments. \emph{Statistical
#' Applications in Genetics and Molecular Biology}, Volume \bold{3}, Article 3.
#' 
#' Smyth, G. K. (2005). Limma: linear models for microarray data. In:
#' \emph{Bioinformatics and Computational Biology Solutions using R and
#' Bioconductor}. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber
#' (eds), Springer, New York, 2005.
#' 
#' Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery
#' rate: a practical and powerful approach to multiple testing. \emph{Journal
#' of the Royal Statistical Society Series}, B, \bold{57}, 289-300.
#' @examples
#' 
#' # Randomly generate data for a 2 group problem with 100 CpG sites and 5 
#' # arrays in each # group. 
#' 
#' y<-matrix(rnorm(1000),ncol=10)
#' 
#' group<-factor(rep(c(1,2),each=5))
#' design<-model.matrix(~group)
#' 
#' # Fit linear model for differential variability
#' vfit<-varFit(y,design,coef=c(1,2))
#' 
#' # Look at top table of results
#' topVar(vfit,coef=2)
#' 
#' @rdname varFit
#' @export varFit
varFit <- function(data,design=NULL,coef=NULL,type=NULL,trend=TRUE,robust=TRUE,
                   weights=NULL)
UseMethod("varFit")

#' @return \code{NULL}
#'
#' @rdname varFit
#' @method varFit MethylSet
#' @export
varFit.MethylSet <- function(data,design=NULL,coef=NULL,type=NULL,trend=TRUE,
                             robust=TRUE,weights=NULL)
{
    message("Extracting M values from MethylSet object.")
    meth<-minfi::getMeth(data)
    unmeth<-minfi::getUnmeth(data)
    Mval<-log2((meth+100)/(unmeth+100))
    varFit.default(data=Mval,design=design,coef=coef,type=type,trend=trend,
                   robust=robust,weights=weights)
}

#' @return \code{NULL}
#'
#' @rdname varFit
#' @method varFit DGEList
#' @export
varFit.DGEList <- function(data,design=NULL,coef=NULL,type=NULL,trend=TRUE,
                           robust=TRUE,weights=NULL)
{
    message("Converting counts to log counts-per-million using voom.")
    v <- limma::voom(data,design)
    varFit.default(data=v$E,design=design,coef=coef,type=type,trend=trend,
                   robust=robust,weights=v$weights)
}

#' @return \code{NULL}
#'
#' @rdname varFit
#' @method varFit default
#' @export
varFit.default <- function(data,design=NULL,coef=NULL,type=NULL,trend=TRUE,
                           robust=TRUE,weights=NULL)
#   Test for differential variability using linear modelling
#   Belinda Phipson
#   15 July 2014. Updated 2 February 2015.
{
#   Check data
    data <- as.matrix(data)
    
#   Check for beta values
    if(min(data) >= 0 & max(data) <= 1){
#   Replace zeros with minimum beta observed > 0    
        zeros <- data==0
        sort.data <- c(sort(data[!zeros]))
        data[zeros] <- min(sort.data)
#   Replace ones with maximum beta observed < 1
        ones <- data==1
        sort.data <- c(sort(data[!ones]))
        data[ones] <- max(sort.data)
#   Perform logit transformation
        data<-log(data/(1-data))
    }

#   Check design, coef and get Levene residuals
    if(is.null(design)){
        z <- getLeveneResiduals(data,design=design,type=type)    
        design <- matrix(1,ncol(data),1)
        coef <- 1
    }    
    else{        
        design <- as.matrix(design)
        if(is.null(coef)){
            message("coef not specified. Using all columns of design matrix.")
            coef <- c(1:ncol(design))
        }
        z <- getLeveneResiduals(data,design=design[,coef],type=type)
    }
        
    fit <- limma::lmFit(z$data,design,weights=weights)
    if(is.null(fit$Amean)) 
        fit$Amean <- rowMeans(z$data, na.rm = TRUE)

    fit <- limma::eBayes(fit,trend=trend,robust=robust)
     
    fit$AvgVar <- z$AvgVar
    fit$LogVarRatio <- matrix(0,ncol=ncol(design),nrow=nrow(fit))
    colnames(fit$LogVarRatio) <- colnames(design)
    fit$LogVarRatio[,coef] <- z$LogVarRatio

    fit
}



#' Compute contrasts for a varFit object.
#' 
#' Compute estimated coefficients, standard errors and LogVarRatios for a given
#' set of contrasts.
#' 
#' This function calls the \code{contrasts.fit} function in \code{limma} to
#' compute coefficients and standard errors for the specified contrasts
#' corresponding to a linear model fit obtained from the \code{varFit}
#' function. LogVarRatios are also computed in terms of the contrasts. A
#' contrasts matrix can be computed using the \code{makeContrasts} function.
#' 
#' @param fit List containing a linear model fit produced by \code{varFit}. The
#' fit object should be of class \code{MArrayLM}.
#' @param contrasts Numeric matrix with rows corresponding to coefficients in
#' \code{fit} and columns containing contrasts.
#' @return A list object of the same class as \code{fit}.
#' @author Belinda Phipson
#' @seealso \code{varFit}, \code{contrasts.fit}, \code{makeContrasts}
#' @examples
#' 
#' # Randomly generate data for a 3 group problem with 100 CpG sites and 4 
#' # arrays in each group. 
#' 
#' library(limma)
#' 
#' y<-matrix(rnorm(1200),ncol=12)
#' 
#' group<-factor(rep(c(1,2,3),each=4))
#' design<-model.matrix(~0+group)
#' colnames(design)<-c("grp1","grp2","grp3")
#' design
#' 
#' # Fit linear model for differential variability
#' # Please always specify the coef parameter in the call to varFit
#' vfit<-varFit(y,design,coef=c(1,2,3))
#' 
#' # Specify contrasts
#' contr<-makeContrasts(grp2-grp1,grp3-grp1,grp3-grp2,levels=colnames(design))
#' 
#' # Compute contrasts from fit object
#' vfit.contr<-contrasts.varFit(vfit,contrasts=contr)
#' 
#' summary(decideTests(vfit.contr))
#' 
#' # Look at top table of results for first contrast
#' topVar(vfit.contr,coef=1)
#' 
#' @export contrasts.varFit
contrasts.varFit <- function(fit, contrasts=NULL)
{
    if(!is.null(fit$Amean))
        trend <- TRUE
    else
        trend <- FALSE
    if(length(fit$df.prior)>1)
        robust <- TRUE
    else
    	robust <- FALSE
    fit.contr <- limma::contrasts.fit(fit, contrasts=contrasts)
    fit.contr <- limma::eBayes(fit.contr, trend=trend, robust=robust)
    
    fit.contr$AvgVar <- fit$AvgVar
    fit.contr$LogVarRatio <- fit$LogVarRatio %*% contrasts
    
    fit.contr
}



#' Table of top-ranked differentially variable CpGs
#' 
#' Extract a table of the top-ranked CpGs from a linear model fit after a
#' differential variability analysis.
#' 
#' This function summarises the results of a differential variability analysis
#' performed with \code{varFit}. The p-values from the comparison of interest
#' are adjusted using Benjamini and Hochberg's false discovery rate with the
#' function \code{p.adjust}. The top ranked CpGs are selected by first ranking
#' the adjusted p-values, then ranking the raw p-values. At this time no other
#' sorting option is catered for.
#' 
#' @param fit List containing a linear model fit produced by \code{varFit}. The
#' fit object should be of class \code{MArrayLM}.
#' @param coef Column number or column name specifying which coefficient of the
#' linear model fit is of interest. It should be the same coefficient that the
#' differential variability testing was performed on. Default is last column of
#' fit object.
#' @param number Maximum number of genes to list. Default is 10.
#' @param sort Logical, default is TRUE. Sorts output according the P-value.
#' FALSE will return results in same order as fit object.
#' @return Produces a dataframe with rows corresponding to the top CpGs and the
#' following columns: \item{genelist }{one or more columns of annotation for
#' each CpG, if the gene information is available in \code{fit}} \item{AvgVar
#' }{average of the absolute or squared Levene residuals across all samples}
#' \item{DiffVar}{estimate of the difference in the Levene residuals
#' corresponding to the comparison of interest} \item{t }{moderated
#' t-statistic} \item{P.Value }{raw p-value} \item{Adj.P.Value }{adjusted
#' p-value}
#' @author Belinda Phipson
#' @seealso \code{varFit}, \code{p.adjust}
#' @references Phipson, B., and Oshlack, A. (2014). A method for detecting
#' differential variability in methylation data shows CpG islands are highly
#' variably methylated in cancers. \emph{Genome Biology}, \bold{15}:465.
#' 
#' Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery
#' rate: a practical and powerful approach to multiple testing. \emph{Journal
#' of the Royal Statistical Society Series}, B, \bold{57}, 289-300.
#' @examples
#' 
#' # Randomly generate data for a 2 group problem with 100 CpG sites and 5 
#' # arrays in each group. 
#' 
#' y<-matrix(rnorm(1000),ncol=10)
#' 
#' group<-factor(rep(c(1,2),each=5))
#' design<-model.matrix(~group)
#' 
#' # Fit linear model for differential variability
#' vfit<-varFit(y,design)
#' 
#' # Look at top table of results
#' topVar(vfit,coef=2)
#' 
#' @export topVar
topVar <- function(fit,coef = NULL,number=10,sort=TRUE)
{
    if(is.null(coef)) coef = ncol(fit)
    p.adj <- stats::p.adjust(fit$p.value[,coef],method="BH")
    if(!is.null(fit$genes)){
        if(is.null(fit$LogVarRatio))
            out <- data.frame(fit$genes, SampleVar = fit$AvgVar, 
                              DiffLevene =  fit$coeff[,coef], t = fit$t[,coef], 
                              P.Value = fit$p.value[,coef], Adj.P.Value = p.adj)
        else    
            out <- data.frame(fit$genes, SampleVar = fit$AvgVar, 
                              LogVarRatio = fit$LogVarRatio[,coef], 
                              DiffLevene =  fit$coeff[,coef], t = fit$t[,coef], 
                              P.Value = fit$p.value[,coef], Adj.P.Value = p.adj)
    }
    else{
        if(is.null(fit$LogVarRatio))
            out <- data.frame(SampleVar = fit$AvgVar, 
                              DiffLevene =  fit$coeff[,coef], t = fit$t[,coef], 
                              P.Value = fit$p.value[,coef], Adj.P.Value = p.adj)
        else
            out <- data.frame(SampleVar = fit$AvgVar, 
                              LogVarRatio = fit$LogVarRatio[,coef], 
                              DiffLevene =  fit$coeff[,coef], t = fit$t[,coef], 
                              P.Value = fit$p.value[,coef], Adj.P.Value = p.adj)
    }
    if(sort){
        o <- order(p.adj,out$P.Value) 
        out[o,][1:number,]
    }
    else 
        out[1:number,]
}



#' Obtain Levene residuals
#' 
#' Obtain absolute or squared Levene residuals for each CpG given a series of
#' methylation arrays
#' 
#' This function will return absolute or squared Levene residuals given a
#' matrix of M values and a design matrix. This can be used for graphing
#' purposes or for downstream analysis such a gene set testing based on
#' differential variability rather than differential methylation. If no design
#' matrix is given, the residuals are determined by treating all samples as
#' coming from one group.
#' 
#' @param data Object of class \code{matrix} of M values, with rows
#' corresponding to features of interest such as CpG sites and columns
#' corresponding to samples or arrays.
#' @param design The design matrix of the experiment, with rows corresponding
#' to arrays/samples and columns to coefficients to be estimated. Defaults to
#' the unit vector.
#' @param type Character string, \code{"AD"} for absolute residuals or
#' \code{"SQ"} for squared residuals. Default is \code{"AD"}.
#' @return Returns a list with three components. \code{data} contains a matrix
#' of absolute or squared residuals, \code{AvgVar} is a vector of sample
#' variances and \code{LogVarRatio} corresponds to the columns of the design
#' matrix and is usually the ratios of the log of the group variances.
#' @author Belinda Phipson
#' @seealso \code{\link{varFit}}
#' @references Phipson, B., and Oshlack, A. (2014). A method for detecting
#' differential variability in methylation data shows CpG islands are highly
#' variably methylated in cancers. \emph{Genome Biology}, \bold{15}:465.
#' @examples
#' 
#' # Randomly generate data for a 2 group problem with 100 CpG sites and 5 
#' # arrays in each group
#' y <- matrix(rnorm(1000),ncol=10)
#' 
#' group <- factor(rep(c(1,2),each=5))
#' design <- model.matrix(~group)
#' 
#' # Get absolute Levene Residuals
#' resid <- getLeveneResiduals(y,design)
#' 
#' # Plot the first CpG
#' barplot(resid$data[1,],col=rep(c(2,4),each=5),
#' ylab="Absolute Levene Residuals",names=group)
#' 
#' @export getLeveneResiduals
getLeveneResiduals<-function(data,design=NULL,type=NULL)
{
#   Check data
    data <- as.matrix(data)

#   Calculate variance across all samples
    y <- rowMeans(data)
    AvgVar <- rowSums((data-y)^2)/(ncol(data)-1)
    
#   Absolute levine residuals or squared levine residuals? 
#   Default is absolute levine residuals.
    if(is.null(type)) type = "AD"

#   Check design
    if(is.null(design)){
        message("Design matrix not specified. Treating all samples as one group.")
        if(type=="AD")
            z <- abs(data-y)
        if(type=="SQ")
            z <- (data-y)^2
        list(data=z,AvgVar=AvgVar,LogVarRatio=NULL)
    }    

    else{ 
#   Calculating means based on design matrix    
        QR <- qr(design)
        V.inv <- chol2inv(QR$qr,size=QR$rank)
       	des.y <- data %*% design
	beta.hat <- des.y %*% t(V.inv)
	ybar <- beta.hat %*% t(design)
		
#   Getting sample sizes to calculate leverage factors
        n.inv <- diag(design %*% V.inv %*% t(design))
        n <- 1/n.inv
        lvg <- sqrt(n/(n-1))
        lvg <- matrix(rep(lvg,nrow(data)),byrow=TRUE,nrow=nrow(data))
        
#   Calculating group variances to get LogVarRatios        
        z <- (data-ybar)^2
        a <- z %*% design
        b <- t(V.inv) %*% t(design)
	var.mat <- (a %*% b)*n/(n-1)
	c <- design %*% t(V.inv)
	LogVarRatio <- log(var.mat) %*% c
	colnames(LogVarRatio) <- colnames(design)
	        
#   Calculating Levene residuals        
        if(type=="AD"){
            z <- abs(data-ybar)
            z <- z*lvg
        }
        if(type=="SQ"){
            z <- z*lvg
        }
        list(data=z,AvgVar=AvgVar,LogVarRatio=LogVarRatio)
    }
}

