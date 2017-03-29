SWAN <- function(data, verbose = FALSE)
UseMethod("SWAN")

SWAN.MethyLumiSet <- function(data, verbose = FALSE){
    
    cat("[SWAN] MethyLumiSet -> MethylSet\n")    
    M <- methylated(data)
    U <- unmethylated(data)
    
    mSet <- MethylSet(Meth = M, Unmeth = U, colData = DataFrame(phenoData(data)@data),
                 annotation = annotation(data))
    mSet@preprocessMethod <- c(rg.norm = "Raw (no normalization or bg correction)",
                                minfi = as.character(packageVersion("minfi")), 
                                manifest = as.character(packageVersion("IlluminaHumanMethylation450kmanifest")))
    SWAN.default(data = mSet, verbose = verbose)
}

SWAN.RGChannelSet <- function(data, verbose = FALSE){
    
    cat("[SWAN] RGChannelSet -> MethylSet\n")
    mSet <- preprocessRaw(data)
    SWAN.default(data = mSet, verbose = verbose)
}

SWAN.default <- function(data, verbose = FALSE){
  
    if(!is(data, "MethylSet")) stop("'data' should be a 'MethylSet'")

    cat("[SWAN] Preparing normalization subset\n")

    if(annotation(data)["array"] == "IlluminaHumanMethylation450k"){
      cat("450k\n")
      manifest <- IlluminaHumanMethylation450kmanifest
      manifestName <- "IlluminaHumanMethylation450kmanifest"
      
    } else if (annotation(data)["array"] == "IlluminaHumanMethylationEPIC") {
      cat("EPIC\n")
      manifest <- IlluminaHumanMethylationEPICmanifest
      manifestName <- "IlluminaHumanMethylationEPICmanifest"
      
    }
    
    CpG.counts <- rbind(data.frame(getProbeInfo(manifest, type = "I")[, c("Name", "nCpG")], Type="I"),
                        data.frame(getProbeInfo(manifest, type = "II")[, c("Name", "nCpG")], Type="II"))
    
    inBoth <- intersect(getManifestInfo(manifest,type="locusNames"), 
                        featureNames(data))
    CpG.counts <- CpG.counts[CpG.counts$Name %in% inBoth,]
                          
    subset <- min(table(CpG.counts$nCpG[CpG.counts$Type == "I" & CpG.counts$nCpG %in% 1:3]),
                    table(CpG.counts$nCpG[CpG.counts$Type == "II" & CpG.counts$nCpG %in% 1:3]))
    xNormSet <- vector("list", 2)
    xNormSet[[1]] <- .getSubset(CpG.counts$nCpG[CpG.counts$Type=="I"], subset)
    xNormSet[[2]] <- .getSubset(CpG.counts$nCpG[CpG.counts$Type=="II"], subset)
    
    data <- data[inBoth, ]
    methData <- getMeth(data)
    unmethData <- getUnmeth(data)
    normMethData <- matrix(NA_real_, ncol = ncol(methData), nrow = nrow(methData))
    normUnmethData <- normMethData
    
    cat("[SWAN] Normalizing methylated channel\n")
    normMethData <- sapply(1:ncol(data), function(i){
        if(verbose) cat(sprintf("[SWAN] Normalizing array %d of %d\n", i, ncol(data)))
        .normalizeTypes(methData[CpG.counts$Name[CpG.counts$Type=="I"], i], 
                        methData[CpG.counts$Name[CpG.counts$Type=="II"], i], xNormSet)})
    
    cat("[SWAN] Normalizing unmethylated channel\n")
    normUnmethData <- sapply(1:ncol(data), function(i){
        if(verbose) cat(sprintf("[SWAN] Normalizing array %d of %d\n", i, ncol(data)))
        .normalizeTypes(unmethData[CpG.counts$Name[CpG.counts$Type=="I"], i],
                        unmethData[CpG.counts$Name[CpG.counts$Type=="II"], i], xNormSet)})

    normSet <- MethylSet(Meth = normMethData, Unmeth = normUnmethData, colData = colData(data),
                  annotation = annotation(data))
 
    featureNames(normSet) <- featureNames(data)
    sampleNames(normSet) <- sampleNames(data)
    normSet@preprocessMethod <- c(rg.norm = sprintf("SWAN (based on a MethylSet preprocessed as '%s')",
                                                      preprocessMethod(data)[1]),
                                    minfi = as.character(packageVersion("minfi")),
                                    manifest = as.character(packageVersion(manifestName)))
    normSet
}

.getSubset <- function(counts, subset){
  
  x <- unlist(lapply(1:3, function(i) sample(which(counts == i), subset)))
  return(1:length(counts) %in% x)
}

.normalizeTypes <- function(intensityI, intensityII, xNormSet) {
    xTarget <- .aveQuantile(list(intensityI[xNormSet[[1]]], intensityII[xNormSet[[2]]]))
    xNorm <- unlist(.subsetQuantileNorm(list(intensityI, intensityII), xNormSet, xTarget))
    xNorm
}

.aveQuantile <- function(X) {
    nbrOfChannels <- length(X)
    if (nbrOfChannels == 1) {
          return(X)
    }
    nbrOfObservations <- unlist(lapply(X, FUN = length), use.names = FALSE)
    maxNbrOfObservations <- max(nbrOfObservations)
    if (maxNbrOfObservations == 1) {
          return(X)
    }
    ## nbrOfFiniteObservations <- rep(maxNbrOfObservations, times = nbrOfChannels)
    quantiles <- (0:(maxNbrOfObservations - 1))/(maxNbrOfObservations - 1)
    xTarget <- vector("double", maxNbrOfObservations)
    for (cc in 1:nbrOfChannels) {
          Xcc <- X[[cc]]
          Scc <- sort(Xcc)
          nobs <- length(Scc)
          if (nobs < maxNbrOfObservations) {
              ## tt <- !is.na(Xcc)
              bins <- (0:(nobs - 1))/(nobs - 1)
              Scc <- approx(x = bins, y = Scc, xout = quantiles,ties = "ordered")$y
          }
          xTarget <- xTarget + Scc
    }
    rm(Scc, Xcc)
    xTarget <- xTarget/nbrOfChannels
    xTarget
}

.subsetQuantileNorm <- function(x, xNormSet, xTarget) {
  
  bgMin <- min(unlist(x)[unlist(x) != 0])
  if(bgMin < 100) bg <- bgMin else bg <- 100
  
  for(i in 1:length(x)){
    
     n <- length(x[[i]])
     nTarget <- length(xTarget)
     nNormSet <- sum(xNormSet[[i]])
          
     if(nNormSet != nTarget){
         
         targetQuantiles <- (0:(nTarget - 1))/(nTarget - 1)
         r <- rank(x[xNormSet[,i], i])
         xNew <-(r - 1)/(nNormSet - 1)
         xNew <- xNew[order(xNew)]
         xNorm <- approx(x = targetQuantiles, y = xTarget, xout = xNew, ties = "ordered", rule = 2)$y
     } else {
         
         xNorm <- xTarget
     }
          
     r <- rank(x[[i]])
     xNew <- (r - 1)/(n - 1)
     quantiles <- xNew[xNormSet[[i]]]
     quantiles <- quantiles[order(quantiles)]
     xmin <- min(x[[i]][xNormSet[[i]]]) #get min value from subset
     xmax <- max(x[[i]][xNormSet[[i]]]) #get max value from subset
     kmax <- which(xNew > max(quantiles)) 
     kmin<- which(xNew < min(quantiles))
     offsets.max <- x[[i]][kmax] - xmax
     offsets.min <- x[[i]][kmin] - xmin
     x[[i]] <- approx(x = quantiles, y = xNorm, xout = xNew, ties = "ordered")$y #interpolate
     x[[i]][kmax] <- max(xNorm) + offsets.max
     x[[i]][kmin] <- min(xNorm) + offsets.min
     x[[i]] <- ifelse(x[[i]] <= 0, bg, x[[i]])    
  }
  
  x
}


