#' Plot the beta value distributions of the Infinium I and II probe types
#' relative to the overall beta value distribution.
#' 
#' Plot the overall density distribution of beta values and the density
#' distributions of the Infinium I and II probe types.
#' 
#' The density distribution of the beta values for a single sample is plotted.
#' The density distributions of the Infinium I and II probes are then plotted
#' individually, showing how they contribute to the overall distribution. This
#' is useful for visualising how using \code{\link{SWAN}} affects the data.
#' 
#' @param data A \code{MethylSet} or a \code{matrix} or a \code{vector}.  We
#' either use the \code{getBeta} function to get Beta values (in the first
#' case) or we assume the matrix or vector contains Beta values.
#' @param legendPos The x and y co-ordinates to be used to position the legend.
#' They can be specified by keyword or in any way which is accepted by
#' \code{\link[grDevices]{xy.coords}}. See \code{\link[graphics]{legend}} for
#' details.
#' @param colors Colors to be used for the different beta value density
#' distributions. Must be a vector of length 3.
#' @param main Plot title.
#' @param lwd The line width to be used for the different beta value density
#' distributions.
#' @param cex.legend The character expansion factor for the legend text.
#' @author Jovana Maksimovic
#' @seealso \code{\link{densityPlot}}, \code{\link{densityBeanPlot}},
#' \code{\link[graphics]{par}}, \code{\link[graphics]{legend}}
#' @return No return value. Plot is produced as a side-effect.
#' @examples
#' 
#' if (require(minfi) & require(minfiData)) {
#'   dat <- preprocessRaw(RGsetEx)
#'   datSwan <- SWAN(dat)
#'   par(mfrow=c(1,2))
#'   densityByProbeType(dat[,1], main="Raw")
#'   densityByProbeType(datSwan[,1], main="SWAN")
#' }
#' 
#' @export densityByProbeType
densityByProbeType <- function(data, legendPos = "top", 
                               colors = c("black", "red", "blue"),
                               main = "", lwd = 3, cex.legend = 1){
    if(is(data, "MethylSet")) {
        if(ncol(data) > 1){
            stop("'data' must only contain one sample")
            
        } else {
            betas <- matrix(minfi::getBeta(data), ncol=1, 
                            dimnames=list(featureNames(data), 
                                          sampleNames(data)))
            
        }

    } else {
        if(is.vector(data) | is.matrix(data)){
            r = range(data, na.rm=TRUE)
            
            if(!(r[1] >= 0 & r[2] <= 1)){
                stop("'data' needs to be a 'MethylSet' or a matrix or vector of 
                     beta values")
                
            } else {
                if (is.vector(data)) {
                    betas <- matrix(data, ncol=1,dimnames=list(names(data),NULL))
                    
                } else {
                    if(ncol(data) > 1){
                        stop("'data' must only contain only one sample")
                        
                    } else {
                        betas <- matrix(data,ncol=1,
                                        dimnames=list(rownames(data),
                                                      colnames(data)))
                        
                    }
                }
            }          
        } else {
            stop("'data' needs to be either a 'MethylSet', 'vector', 'matrix' or 
                 'data.frame'")
            
        } 
    }
    
  if(annotation(data)["array"] == "IlluminaHumanMethylation450k"){
    manifest <- IlluminaHumanMethylation450kmanifest::IlluminaHumanMethylation450kmanifest
  } else if (annotation(data)["array"] == "IlluminaHumanMethylationEPIC") {
    manifest <- IlluminaHumanMethylationEPICmanifest::IlluminaHumanMethylationEPICmanifest
  }
    probeTypes <- rbind(data.frame(minfi::getProbeInfo(manifest, 
                                                type = "I")[, c("Name", 
                                                                "nCpG")], 
                                   Type="I"),
                        data.frame(minfi::getProbeInfo(manifest, 
                                                type = "II")[, c("Name", 
                                                                 "nCpG")], 
                                   Type="II"))
    
    ymax <- max(sapply(1:ncol(betas), function(x) max(stats::density(betas[,x],
                                                              na.rm=TRUE)$y)))
    betas <- matrix(betas[!is.na(betas),],ncol=1,
                    dimnames=list(rownames(betas)[!is.na(betas)],
                                  colnames(betas)))
    total <- nrow(betas)
    type1 <- sum(probeTypes$Type=="I" & probeTypes$Name %in% rownames(betas))
    type2 <- sum(probeTypes$Type=="II" & probeTypes$Name %in% rownames(betas))

    graphics::plot(stats::density(betas), main=main, xlab="Beta values", 
                   ylim=c(0,ymax), lwd=lwd, col=colors[1])
    graphics::lines(suppressWarnings(stats::density(betas[rownames(betas) %in% 
                                    probeTypes$Name[probeTypes$Type=="I"]],
                                   weights=rep(1/total,type1))), col=colors[2], 
          lty=2, lwd=lwd)
    graphics::lines(suppressWarnings(stats::density(betas[rownames(betas) %in% 
                                    probeTypes$Name[probeTypes$Type=="II"]],
                                   weights=rep(1/total,type2))), col=colors[3], 
          lty=2, lwd=lwd)
    graphics::legend(legendPos, c("All probes", "Infinium I", "Infinium II"), 
                     col=colors, lwd=lwd, bg="white", cex=cex.legend, 
                     lty=c(1,2,2))
}
