#' Subset-quantile Within Array Normalisation for Illumina Infinium
#' HumanMethylation450 BeadChips
#' 
#' Subset-quantile Within Array Normalisation (SWAN) is a within array
#' normalisation method for the Illumina Infinium HumanMethylation450 platform.
#' It allows Infinium I and II type probes on a single array to be normalized
#' together.
#' 
#' The SWAN method has two parts. First, an average quantile distribution is
#' created using a subset of probes defined to be biologically similar based on
#' the number of CpGs underlying the probe body. This is achieved by randomly
#' selecting N Infinium I and II probes that have 1, 2 and 3 underlying CpGs,
#' where N is the minimum number of probes in the 6 sets of Infinium I and II
#' probes with 1, 2 or 3 probe body CpGs. If no probes have previously been
#' filtered out e.g. sex chromosome probes, etc. N=11,303. This results in a
#' pool of 3N Infinium I and 3N Infinium II probes. The subset for each probe
#' type is then sorted by increasing intensity. The value of each of the 3N
#' pairs of observations is subsequently assigned to be the mean intensity of
#' the two probe types for that row or 'quantile'. This is the standard
#' quantile procedure. The intensities of the remaining probes are then
#' separately adjusted for each probe type using linear interpolation between
#' the subset probes.
#' 
#' @aliases SWAN SWAN.default SWAN.MethyLumiSet SWAN.RGChannelSet
#' @param data An object of class either \code{MethylSet}, \code{RGChannelSet}
#' or \code{MethyLumiSet}.
#' @param verbose Should the function be verbose?
#' @return An object of class \code{MethylSet}.
#' @note SWAN uses a random subset of probes to perform the within-array
#' normalization. In order to achive reproducible results, the seed needs to be
#' set using \code{set.seed}.
#' @author Jovana Maksimovic 
#' @seealso \code{\linkS4class{RGChannelSet}} and
#' \code{\linkS4class{MethylSet}} as well as \code{\linkS4class{MethyLumiSet}}
#' and \code{\linkS4class{IlluminaMethylationManifest}}.
#' @references J Maksimovic, L Gordon and A Oshlack (2012).  \emph{SWAN: Subset
#' quantile Within-Array Normalization for Illumina Infinium
#' HumanMethylation450 BeadChips}.  Genome Biology 13, R44.
#' @examples
#' 
#' if (require(minfi) & require(minfiData)) {
#' 
#'   set.seed(100)
#'   datSwan1 <- SWAN(RGsetEx)
#'   
#'   dat <- preprocessRaw(RGsetEx)
#'   set.seed(100)
#'   datSwan2 <- SWAN(dat)
#'   
#'   head(getMeth(datSwan2)) == head(getMeth(datSwan1))
#' }
#' 
#' @rdname SWAN
#' @export SWAN
SWAN <- function(data, verbose = FALSE)
UseMethod("SWAN")

#' @return \code{NULL}
#'
#' @rdname SWAN
#' @method SWAN MethyLumiSet
#' @export
SWAN.MethyLumiSet <- function(data, verbose = FALSE){
    
    cat("[SWAN] MethyLumiSet -> MethylSet\n")    
    M <- methylumi::methylated(data)
    U <- methylumi::unmethylated(data)
    
    mSet <- minfi::MethylSet(Meth = M, Unmeth = U, 
                      colData = DataFrame(phenoData(data)@data),
                      annotation = annotation(data))
    mSet@preprocessMethod <- c(rg.norm = "Raw (no normalization or bg 
                               correction)",
                                minfi = as.character(packageVersion("minfi")), 
manifest = as.character(packageVersion("IlluminaHumanMethylation450kmanifest")))
    SWAN.default(data = mSet, verbose = verbose)
}

#' @return \code{NULL}
#'
#' @rdname SWAN
#' @method SWAN RGChannelSet
#' @export
SWAN.RGChannelSet <- function(data, verbose = FALSE){
    
    cat("[SWAN] RGChannelSet -> MethylSet\n")
    mSet <- minfi::preprocessRaw(data)
    SWAN.default(data = mSet, verbose = verbose)
}

#' @return \code{NULL}
#'
#' @rdname SWAN
#' @method SWAN default
#' @export
SWAN.default <- function(data, verbose = FALSE){
  
    if(!is(data, "MethylSet")) stop("'data' should be a 'MethylSet'")

    cat("[SWAN] Preparing normalization subset\n")

    if(annotation(data)["array"] == "IlluminaHumanMethylation450k"){
      cat("450k\n")
      manifest <- minfi::IlluminaHumanMethylation450kmanifest
      manifestName <- "IlluminaHumanMethylation450kmanifest"
      
    } else if (annotation(data)["array"] == "IlluminaHumanMethylationEPIC") {
      cat("EPIC\n")
      manifest <- minfi::IlluminaHumanMethylationEPICmanifest
      manifestName <- "IlluminaHumanMethylationEPICmanifest"
      
    }
    
    CpG.counts <- rbind(data.frame(minfi::getProbeInfo(manifest, 
                                                type = "I")[, c("Name", 
                                                                "nCpG")], 
                                   Type="I"),
                        data.frame(minfi::getProbeInfo(manifest, 
                                                type = "II")[, c("Name", 
                                                                 "nCpG")], 
                                   Type="II"))
    
    inBoth <- intersect(minfi::getManifestInfo(manifest,type="locusNames"), 
                        featureNames(data))
    CpG.counts <- CpG.counts[CpG.counts$Name %in% inBoth,]
                          
    subset <- min(table(CpG.counts$nCpG[CpG.counts$Type == "I" & 
                                            CpG.counts$nCpG %in% 1:3]),
                    table(CpG.counts$nCpG[CpG.counts$Type == "II" & 
                                              CpG.counts$nCpG %in% 1:3]))
    xNormSet <- vector("list", 2)
    xNormSet[[1]] <- .getSubset(CpG.counts$nCpG[CpG.counts$Type=="I"], subset)
    xNormSet[[2]] <- .getSubset(CpG.counts$nCpG[CpG.counts$Type=="II"], subset)
    
    data <- data[inBoth, ]
    methData <- minfi::getMeth(data)
    unmethData <- minfi::getUnmeth(data)
    normMethData <- matrix(NA_real_, ncol = ncol(methData), 
                           nrow = nrow(methData))
    normUnmethData <- normMethData
    
    cat("[SWAN] Normalizing methylated channel\n")
    normMethData <- sapply(1:ncol(data), function(i){
        if(verbose) cat(sprintf("[SWAN] Normalizing array %d of %d\n", i, 
                                ncol(data)))
        .normalizeTypes(methData[CpG.counts$Name[CpG.counts$Type=="I"], i], 
                        methData[CpG.counts$Name[CpG.counts$Type=="II"], i], 
                        xNormSet)})
    
    cat("[SWAN] Normalizing unmethylated channel\n")
    normUnmethData <- sapply(1:ncol(data), function(i){
        if(verbose) cat(sprintf("[SWAN] Normalizing array %d of %d\n", i, 
                                ncol(data)))
        .normalizeTypes(unmethData[CpG.counts$Name[CpG.counts$Type=="I"], i],
                        unmethData[CpG.counts$Name[CpG.counts$Type=="II"], i], 
                        xNormSet)})

    normSet <- minfi::MethylSet(Meth = normMethData, Unmeth = normUnmethData, 
                         colData = colData(data),
                  annotation = annotation(data))
 
    featureNames(normSet) <- featureNames(data)
    sampleNames(normSet) <- sampleNames(data)
    normSet@preprocessMethod <- c(rg.norm = sprintf("SWAN (based on a MethylSet 
                                                    preprocessed as '%s')",
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
    xTarget <- .aveQuantile(list(intensityI[xNormSet[[1]]], 
                                 intensityII[xNormSet[[2]]]))
    xNorm <- unlist(.subsetQuantileNorm(list(intensityI, intensityII), 
                                        xNormSet, xTarget))
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
    ## nbrOfFiniteObservations <- rep(maxNbrOfObservations,
    ## times = nbrOfChannels)
    quantiles <- (0:(maxNbrOfObservations - 1))/(maxNbrOfObservations - 1)
    xTarget <- vector("double", maxNbrOfObservations)
    for (cc in 1:nbrOfChannels) {
          Xcc <- X[[cc]]
          Scc <- sort(Xcc)
          nobs <- length(Scc)
          if (nobs < maxNbrOfObservations) {
              ## tt <- !is.na(Xcc)
              bins <- (0:(nobs - 1))/(nobs - 1)
              Scc <- approx(x = bins, y = Scc, xout = quantiles,
                            ties = "ordered")$y
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
         xNorm <- approx(x = targetQuantiles, y = xTarget, xout = xNew, 
                         ties = "ordered", rule = 2)$y
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
     x[[i]] <- approx(x = quantiles, y = xNorm, xout = xNew, 
                      ties = "ordered")$y #interpolate
     x[[i]][kmax] <- max(xNorm) + offsets.max
     x[[i]][kmin] <- min(xNorm) + offsets.min
     x[[i]] <- ifelse(x[[i]] <= 0, bg, x[[i]])    
  }
  
  x
}


