RUVfit <- function(data, design, coef=ncol(design), ctl, method=c("inv", "rinv", "ruv4", "ruv2"), k = NULL, ...){
    
  method <- match.arg(method)
  
  if ((method %in% c("ruv4", "ruv2")) & is.null(k))
    stop("'k' cannot be NULL if method is 'ruv4' or 'ruv2'.")
  
  design <- as.matrix(design)
  
  if (mode(design) != "numeric") 
    stop("'design' must be a numeric matrix.")
  
  if (mode(coef) != "numeric")
    stop("'coef' must be an integer or numeric vector.")
  
  if (coef < 1 | coef > ncol(design))
    stop("'coef' can only contain values >= 1 or <= ncol 'design'.")
  
  if (mode(ctl) != "logical")
    stop("'ctl' must be a logical vector.")
  
  X <- as.matrix(design[, coef])
  Z <- as.matrix(design[, -coef])

  Y <- as.matrix(data)
  
  if (mode(Y) != "numeric") 
    stop("'data' must be a numeric matrix.")
  
  Y <- t(Y)
    
  fit <- switch(method, 
                inv = RUVinv(Y = Y, X = X, ctl = ctl, Z = Z, ...),
                rinv = RUVrinv(Y = Y, X = X, ctl = ctl, Z = Z, k = k, ...),
                ruv4 = RUV4(Y = Y, X = X, ctl = ctl, k = k, Z = Z, ...),
                ruv2 = RUV2(Y = Y, X = X, ctl = ctl, k = k, Z = Z, ...))

  return(.toMArrayLM(fit))
  
}

.toMArrayLM <- function(fit){
    
    obj <- new("MArrayLM")
    
    obj$coefficients <- t(fit$betahat)
    obj$sigma2 <- fit$sigma2
    obj$tvals <- fit$tvals
    obj$pvals <- fit$pvals 
    obj$multiplier <- fit$multiplier
    obj$df <- fit$df 
    obj$W <- t(fit$W) 
    obj$alpha <- t(fit$alpha) 
    obj$byx <- t(fit$byx) 
    obj$bwx <- t(fit$bwx) 
    obj$X <- fit$X 
    obj$k <- fit$k
    obj$ctl <- fit$ctl
    obj$Z <- fit$Z 
    obj$fullW0 <- t(fit$fullW0) 
    obj$lambda <- fit$lambda
    obj$t <- t(fit$t)
    obj$p <- t(fit$p)
    
    slots <- attributes(fit)$names
    
    if ("p.rsvar" %in% slots) 
        obj$p.rsvar <- t(fit$p.rsvar)
    
    if ("p.evar" %in% slots) 
        obj$p.evar <- t(fit$p.evar)
    
    if ("p.ebayes" %in% slots) 
        obj$p.ebayes <- t(fit$p.ebayes)
    
    if ("p.rsvar.ebayes" %in% slots) 
        obj$p.rsvar.ebayes <- t(fit$p.rsvar.ebayes)
    
    if ("p.BH" %in% slots) 
        obj$p.BH <- t(fit$p.BH)
    
    if ("p.rsvar.BH" %in% slots) 
        obj$p.rsvar.BH <- t(fit$p.rsvar.BH)
    
    if ("p.evar.BH" %in% slots) 
        obj$p.evar.BH <- t(fit$p.evar.BH)
    
    if ("p.ebayes.BH" %in% slots) 
        obj$p.ebayes.BH <- t(fit$p.ebayes.BH)
    
    if ("p.rsvar.ebayes.BH" %in% slots) 
        obj$p.rsvar.ebayes.BH <- t(fit$p.rsvar.ebayes.BH)
    
    obj
}

.toList <- function(fit){

    obj <- list()
    
    obj$betahat <- t(fit$coefficients)
    obj$sigma2 <- fit$sigma2
    obj$tvals <- fit$tvals
    obj$pvals <- fit$pvals 
    obj$multiplier <- fit$multiplier
    obj$df <- fit$df 
    obj$W <- t(fit$W) 
    obj$alpha <- t(fit$alpha) 
    obj$byx <- t(fit$byx) 
    obj$bwx <- t(fit$bwx) 
    obj$X <- fit$X 
    obj$k <- fit$k
    obj$ctl <- fit$ctl
    obj$Z <- fit$Z 
    obj$fullW0 <- t(fit$fullW0) 
    obj$lambda <- fit$lambda
    obj$t <- t(fit$t)
    obj$p <- t(fit$p)
    
    obj
}

RUVadj <- function(fit, ebayes = TRUE, evar = FALSE, rsvar = FALSE, ...){
    
    fit <- variance_adjust(fit = .toList(fit), ebayes = ebayes, evar = evar, rsvar = rsvar, ...)
    
    return(.toMArrayLM(fit))
}


topRUV <- function (fit, number = 10, p.value.cut = 1,
                         cut.on = c("p.ebayes.BH","p.BH","p.rsvar.BH","p.evar.BH","p.rsvar.ebayes.BH"),
                         sort.by = c("p.ebayes.BH","p.BH","p.rsvar.BH","p.evar.BH","p.rsvar.ebayes.BH")){
  
  coefficients <- fit$coefficients
  t <- fit$t
  p <- fit$p
  ID <- rownames(coefficients)

  p.rsvar <- NULL
  p.evar <- NULL
  p.ebayes <- NULL
  p.rsvar.ebayes <- NULL
  p.BH <- NULL
  p.rsvar.BH <- NULL
  p.evar.BH <- NULL
  p.ebayes.BH <- NULL
  p.rsvar.ebayes.BH <- NULL
  
  slots <- attributes(fit)$names
  
  if ("p.rsvar" %in% slots) 
    p.rsvar <- fit$p.rsvar
  
  if ("p.evar" %in% slots) 
    p.evar <- fit$p.evar
  
  if ("p.ebayes" %in% slots) 
    p.ebayes <- fit$p.ebayes
  
  if ("p.rsvar.ebayes" %in% slots) 
    p.rsvar.ebayes <- fit$p.rsvar.ebayes
  
  if ("p.BH" %in% slots) 
    p.BH <- fit$p.BH
  
  if ("p.rsvar.BH" %in% slots) 
    p.rsvar.BH <- fit$p.rsvar.BH
  
  if ("p.evar.BH" %in% slots) 
    p.evar.BH <- fit$p.evar.BH
  
  if ("p.ebayes.BH" %in% slots) 
    p.ebayes.BH <- fit$p.ebayes.BH
  
  if ("p.rsvar.ebayes.BH" %in% slots) 
    p.rsvar.ebayes.BH <- fit$p.rsvar.ebayes.BH
  
  tab <- data.frame(row.names = ID, coefficients = coefficients, t = t, p = p, stringsAsFactors = FALSE, check.names=FALSE)
  
  if (!is.null(p.BH)) 
    tab <- data.frame(tab, p.BH = p.BH, stringsAsFactors = FALSE, check.names=FALSE)
  
  if (!is.null(p.rsvar)) 
    tab <- data.frame(tab, p.rsvar = p.rsvar, stringsAsFactors = FALSE, check.names=FALSE)
  
  if (!is.null(p.rsvar.BH)) 
    tab <- data.frame(tab, p.rsvar.BH = p.rsvar.BH, stringsAsFactors = FALSE, check.names=FALSE)
  
  if (!is.null(p.evar)) 
    tab <- data.frame(tab, p.evar = p.evar, stringsAsFactors = FALSE, check.names=FALSE)
  
  if (!is.null(p.evar.BH)) 
    tab <- data.frame(tab, p.evar.BH = p.evar.BH, stringsAsFactors = FALSE, check.names=FALSE)
  
  if (!is.null(p.rsvar.ebayes)) 
    tab <- data.frame(tab, p.rsvar.ebayes = p.rsvar.ebayes, stringsAsFactors = FALSE, check.names=FALSE)
  
  if (!is.null(p.rsvar.ebayes.BH)) 
    tab <- data.frame(tab, p.rsvar.ebayes.BH = p.rsvar.ebayes.BH, stringsAsFactors = FALSE, check.names=FALSE)
  
  if (!is.null(p.ebayes)) 
    tab <- data.frame(tab, p.ebayes = p.ebayes, stringsAsFactors = FALSE, check.names=FALSE)
  
  if (!is.null(p.ebayes.BH)) 
    tab <- data.frame(tab, p.ebayes.BH = p.ebayes.BH, stringsAsFactors = FALSE, check.names=FALSE)
  
  if (p.value.cut < 1) {
    
    cut.on <- match.arg(cut.on)
    
    if (!(cut.on %in% slots))
      cat(sprintf("Cannot threshold on '%s' because it is not in the 'fit' object.", cut.on))
    
    adj.p.value <- get(cut.on, fit)
    
    sig <- (adj.p.value < p.value.cut)
    
    if (any(is.na(sig)))
      sig[is.na(sig)] <- FALSE
    
    if (all(!sig))
      return(data.frame())
    
    tab <- tab[sig,]
  }
                 
  sort.by <- match.arg(sort.by)
  
  if (!(sort.by %in% slots))
    stop(sprintf("Cannot sort by '%s' because it is not in the 'fit' object.", sort.by))
                 
  ord <- switch(sort.by, p.ebayes.BH = order(tab$p.ebayes.BH, tab$p.ebayes, decreasing = FALSE), 
                p.BH = order(tab$p.BH, tab$p, decreasing=FALSE), 
                p.rsvar.BH = order(tab$p.rsvar.BH, tab$p.rsvar, decreasing=FALSE), 
                p.evar.BH = order(tab$p.evar.BH, tab$p.evar, decreasing=FALSE), 
                p.rsvar.ebayes.BH = order(tab$p.rsvar.ebayes.BH, tab$p.rsvar.ebayes, decreasing=FALSE))
  
  if (nrow(tab) < number)
    number <- nrow(tab)
                 
  if (number < 1)
    return(data.frame())
  
  top <- ord[1:number]
  
  tab[top,]
}


