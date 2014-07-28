varFit <- function(data,design=NULL,coef=NULL,type=NULL,trend=TRUE,robust=TRUE)
UseMethod("varFit")

varFit.MethylSet <- function(data,design=NULL,coef=NULL,type=NULL,trend=TRUE,robust=TRUE)
{
    message("Extracting M values from MethylSet object.")
    meth<-getMeth(data)
    unmeth<-getUnmeth(data)
    Mval<-log2((meth+100)/(unmeth+100))
    varFit.default(data=Mval,design=design,coef=coef,type=type,trend=trend,robust=robust)
}

varFit.default <- function(data,design=NULL,coef=NULL,type=NULL,trend=TRUE,robust=TRUE)
#   Test for differential variability using linear modelling
#   Belinda Phipson
#   15 July 2014.
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
    
#   Check design
    if(is.null(design)){
        z <- getLeveneResiduals(data,design=design,coef=coef,type=type)    
        design <- matrix(1,ncol(data),1)
    }    
    else{        
        design <- as.matrix(design)
        if(is.null(coef))
            coef <- ncol(design)    
        z <- getLeveneResiduals(data,design=design,coef=coef,type=type)
    }
        
    fit <- lmFit(z$data,design)
    if(is.null(fit$Amean)) 
        fit$Amean <- rowMeans(z$data, na.rm = TRUE)

    fit <- eBayes(fit,trend=trend,robust=robust)
     
    fit$AvgVar <- z$AvgVar
    fit$LogVarRatio <- z$LogVarRatio

    fit
}

topVar <- function(fit,coef = NULL,number=10,sort=TRUE)
{
    if(is.null(coef)) coef = ncol(fit)
    p.adj <- p.adjust(fit$p.value[,coef],method="BH")
    if(!is.null(fit$genes)){
        if(is.null(fit$LogVarRatio))
            out <- data.frame(fit$genes, SampleVar = fit$AvgVar, DiffLevene =  fit$coeff[,coef], t = fit$t[,coef], P.Value = fit$p.value[,coef], Adj.P.Value = p.adj)
        else    
            out <- data.frame(fit$genes, SampleVar = fit$AvgVar, LogVarRatio = fit$LogVarRatio, DiffLevene =  fit$coeff[,coef], t = fit$t[,coef], P.Value = fit$p.value[,coef], Adj.P.Value = p.adj)
    }
    else{
        if(is.null(fit$LogVarRatio))
            out <- data.frame(SampleVar = fit$AvgVar, DiffLevene =  fit$coeff[,coef], t = fit$t[,coef], P.Value = fit$p.value[,coef], Adj.P.Value = p.adj)
        else
            out <- data.frame(SampleVar = fit$AvgVar, LogVarRatio = fit$LogVarRatio, DiffLevene =  fit$coeff[,coef], t = fit$t[,coef], P.Value = fit$p.value[,coef], Adj.P.Value = p.adj)
    }
    if(sort){
        o <- order(p.adj,out$P.Value) 
        out[o,][1:number,]
    }
    else 
        out[1:number,]
}

getLeveneResiduals<-function(data,design=NULL,coef=NULL,type=NULL)
{
#   Check data
    data <- as.matrix(data)

#   Calculate variance across all samples
    y <- rowMeans(data)
    AvgVar <- rowSums((data-y)^2)/(ncol(data)-1)
    
#   Do you want absolute levine residuals or squared levine residuals? Default is absolute.
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
#   Check coef
        if(is.null(coef))
            coef <- ncol(design)

        o <- order(design[,coef])
        o.design <- design[o,]
        o.data <- data[,o]
        group <- o.design[,coef]
        i.group <- names(table(group))
    
#   Check i.group for length 2
        if(length(i.group)!=2) 
            stop("Design matrix not parameterised correctly. Comparison must be two group.")    
  
#   table always orders i.group "0","1"
  
        i.data1 <- o.data[,group==i.group[1]]
        i.data2 <- o.data[,group==i.group[2]]
  
        if(ncol(i.data1)==1 | ncol(i.data2)==1) stop("Each group must have sample size greater than 1.")
  
        y1 <- rowMeans(i.data1)
        y2 <- rowMeans(i.data2)
                      
#   Calculate leverage according to sample size
        n1 <- ncol(i.data1)
        n2 <- ncol(i.data2)
        l1 <- sqrt(n1/(n1-1))
        l2 <- sqrt(n2/(n2-1))
  
#   Calculate absolute Levene residuals
        if(type=="AD"){
            z1 <- abs(i.data1 - y1)
            var1 <- rowSums(z1^2)/(n1-1)
            z1 <- z1*l1
            z2 <- abs(i.data2 - y2)
            var2 <- rowSums(z2^2)/(n2-1)
            z2 <- z2*l2
        }
        if(type=="SQ"){
            z1 <- (i.data1 - y1)^2
            var1 <- rowSums(z1)/(n1-1)
            z1 <- z1*l1
            z2 <- (i.data2 - y2)^2
            var2 <- rowSums(z2)/(n2-1)
            z2 <- z2*l2
        }
  
#     Get data ready for mod t test and put it back into the original order
        z <- cbind(z1,z2)
        z <- z[,order(o)]
          
        LogVarRatio = log(var2/var1)
          
        list(data=z,AvgVar=AvgVar,LogVarRatio=LogVarRatio)
    }
}
