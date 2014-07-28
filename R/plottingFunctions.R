densityByProbeType <- function(data, legendPos = "top", colors = c("black", "red", "blue"),
                            main = "", lwd = 3, cex.legend = 1){
    if(is(data, "MethylSet")) {
        if(ncol(data) > 1){
            stop("'data' must only contain one sample")
            
        } else {
            betas <- matrix(getBeta(data), ncol=1, dimnames=list(featureNames(data), sampleNames(data)))
            
        }

    } else {
        if(is.vector(data) | is.matrix(data)){
            r = range(data, na.rm=TRUE)
            
            if(!(r[1] >= 0 & r[2] <= 1)){
                stop("'data' needs to be a 'MethylSet' or a matrix or vector of beta values")
                
            } else {
                if (is.vector(data)) {
                    betas <- matrix(data, ncol=1,dimnames=list(names(data),NULL))
                    
                } else {
                    if(ncol(data) > 1){
                        stop("'data' must only contain only one sample")
                        
                    } else {
                        betas <- matrix(data,ncol=1,dimnames=list(rownames(data),colnames(data)))
                        
                    }
                }
            }          
        } else {
            stop("'data' needs to be either a 'MethylSet', 'vector', 'matrix' or 'data.frame'")
            
        } 
    }
    
    rgSet <- new("RGChannelSet", annotation = "IlluminaHumanMethylation450k")
    probeTypes <- rbind(data.frame(getProbeInfo(rgSet, type = "I")[, c("Name", "nCpG")], Type="I"),
                        data.frame(getProbeInfo(rgSet, type = "II")[, c("Name", "nCpG")], Type="II"))
    
    ymax <- max(sapply(1:ncol(betas), function(x) max(density(betas[,x],na.rm=TRUE)$y)))
    betas <- matrix(betas[!is.na(betas),],ncol=1,dimnames=list(rownames(betas)[!is.na(betas)],colnames(betas)))
    total <- nrow(betas)
    type1 <- sum(probeTypes$Type=="I" & probeTypes$Name %in% rownames(betas))
    type2 <- sum(probeTypes$Type=="II" & probeTypes$Name %in% rownames(betas))

    plot(density(betas), main=main, xlab="Beta values", ylim=c(0,ymax), lwd=lwd, col=colors[1])
    lines(suppressWarnings(density(betas[rownames(betas) %in% probeTypes$Name[probeTypes$Type=="I"]],
                                   weights=rep(1/total,type1))), col=colors[2], lty=2, lwd=lwd)
    lines(suppressWarnings(density(betas[rownames(betas) %in% probeTypes$Name[probeTypes$Type=="II"]],
                                   weights=rep(1/total,type2))), col=colors[3], lty=2, lwd=lwd)
    legend(legendPos, c("All probes", "Infinium I", "Infinium II"), col=colors, lwd=lwd, bg="white",
           cex=cex.legend, lty=c(1,2,2))
}
