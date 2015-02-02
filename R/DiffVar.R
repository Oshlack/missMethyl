varFit <- function(data,design=NULL,coef=NULL,type=NULL,trend=TRUE,robust=TRUE,weights=NULL)
UseMethod("varFit")

varFit.MethylSet <- function(data,design=NULL,coef=NULL,type=NULL,trend=TRUE,robust=TRUE,weights=NULL)
{
    message("Extracting M values from MethylSet object.")
    meth<-getMeth(data)
    unmeth<-getUnmeth(data)
    Mval<-log2((meth+100)/(unmeth+100))
    varFit.default(data=Mval,design=design,coef=coef,type=type,trend=trend,robust=robust,weights=weights)
}

varFit.DGEList <- function(data,design=NULL,coef=NULL,type=NULL,trend=TRUE,robust=TRUE,weights=NULL)
{
    message("Converting counts to log counts-per-million using voom.")
    v <- voom(data,design)
    varFit.default(data=v$E,design=design,coef=coef,type=type,trend=trend,robust=robust,weights=v$weights)
}

varFit.default <- function(data,design=NULL,coef=NULL,type=NULL,trend=TRUE,robust=TRUE,weights=NULL)
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
        if(is.null(coef)) 
            coef <- c(1,ncol(design))
        z <- getLeveneResiduals(data,design=design[,coef],type=type)
    }
        
    fit <- lmFit(z$data,design,weights=weights)
    if(is.null(fit$Amean)) 
        fit$Amean <- rowMeans(z$data, na.rm = TRUE)

    fit <- eBayes(fit,trend=trend,robust=robust)
     
    fit$AvgVar <- z$AvgVar
    fit$LogVarRatio <- matrix(0,ncol=ncol(design),nrow=nrow(fit))
    colnames(fit$LogVarRatio) <- colnames(design)
    fit$LogVarRatio[,coef] <- z$LogVarRatio

    fit
}

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
    fit.contr <- contrasts.fit(fit, contrasts=contrasts)
    fit.contr <- eBayes(fit.contr, trend=trend, robust=robust)
    
    fit.contr$AvgVar <- fit$AvgVar
    fit.contr$LogVarRatio <- fit$LogVarRatio %*% contrasts
    
    fit.contr
}

topVar <- function(fit,coef = NULL,number=10,sort=TRUE)
{
    if(is.null(coef)) coef = ncol(fit)
    p.adj <- p.adjust(fit$p.value[,coef],method="BH")
    if(!is.null(fit$genes)){
        if(is.null(fit$LogVarRatio))
            out <- data.frame(fit$genes, SampleVar = fit$AvgVar, DiffLevene =  fit$coeff[,coef], t = fit$t[,coef], P.Value = fit$p.value[,coef], Adj.P.Value = p.adj)
        else    
            out <- data.frame(fit$genes, SampleVar = fit$AvgVar, LogVarRatio = fit$LogVarRatio[,coef], DiffLevene =  fit$coeff[,coef], t = fit$t[,coef], P.Value = fit$p.value[,coef], Adj.P.Value = p.adj)
    }
    else{
        if(is.null(fit$LogVarRatio))
            out <- data.frame(SampleVar = fit$AvgVar, DiffLevene =  fit$coeff[,coef], t = fit$t[,coef], P.Value = fit$p.value[,coef], Adj.P.Value = p.adj)
        else
            out <- data.frame(SampleVar = fit$AvgVar, LogVarRatio = fit$LogVarRatio[,coef], DiffLevene =  fit$coeff[,coef], t = fit$t[,coef], P.Value = fit$p.value[,coef], Adj.P.Value = p.adj)
    }
    if(sort){
        o <- order(p.adj,out$P.Value) 
        out[o,][1:number,]
    }
    else 
        out[1:number,]
}

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

