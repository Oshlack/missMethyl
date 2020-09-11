#' Generalised gene set testing for RNA-seq data
#'
#' Given a user defined list of gene sets, \code{gsaseq} will test whether 
#' significantly differentially expressed genes are enriched in these gene sets.
#' 
#' This function is a generalised version of \code{goana} and \code{kegga} from 
#' the \code{limma} package in that it can take a user-defined list of 
#' differentially expressed genes and perform gene set enrichment analysis, and 
#' is not limited to only testing GO and KEGG categories. It is not as flexible
#' as \code{goana} and \code{kegga}. Please note the vector of differentially
#' expressed genes and list of gene sets must be Entrez Gene IDs.
#' 
#' The \code{gsaseq} function will test for enrichment using a hypergeometric 
#' test if the \code{gene.length} parameter is NULL. If the 
#' \code{gene.length} parameter is supplied then the p-values are derived from 
#' Walllenius' noncentral hypergeometric distribution from the \code{BiasedUrn} 
#' package. Please note that the \code{gene.length} parameter must be in the 
#' same order and of the same length as the \code{universe} parameter.
#' 
#' 
#' @param sig.de Character vector of significant differentially expressed genes
#' to test for gene set enrichment. Must be Entrez Gene ID format.
#' @param universe Character vector of all genes analysed in the experiment. 
#' Must be Entrez Gene ID format.
#' @param collection A list of user specified gene sets to test. Can also be a
#' single character vector gene set. Gene identifiers must be Entrez Gene IDs.
#' @param plot.bias Logical, if true a plot showing gene length bias related to
#' differential expression will be displayed.
#' @param gene.length A vector containing the gene lengths for each gene in the 
#' same order as \code{universe}.
#' @param sort Logical, if TRUE then the output dataframe is sorted by p-value.
#'
#' @return A data frame with a row for each gene set and the following columns:
#' \item{N}{ number of genes in the gene set } \item{DE}{ number of genes that
#' are differentially expressed } \item{P.DE}{ p-value for over-representation
#' of the gene set } \item{FDR}{ False discovery rate, calculated using the
#' method of Benjamini and Hochberg (1995).  }
#' @author Belinda Phipson
#' @seealso \code{\link{goana},\link{kegga},\link{camera},\link{roast}}
#' @export gsaseq
#' 
#' 
#' @examples
#' \dontrun{ # to avoid timeout on Bioconductor build
#' library(org.Hs.eg.db)
#' # Use org.Hs.eg.db to extract GO terms
#' GOtoID <- suppressMessages(select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), 
#'                                   columns=c("ENTREZID","GO"), 
#'                                   keytype="ENTREZID"))
#' head(GOtoID)
#' 
#' # Define the universe as random sample of 20000 genes in the annotation
#' universe <- sample(unique(GOtoID$ENTREZID),20000)
#' 
#' # Randomly sample 500 genes as DE
#' de.genes <- sample(universe, 500)
#' 
#' # Generate random gene lengths for genes in universe
#' # This is based on the true distribution of log(gene length) of genes in the
#' # hg19 genome
#' logGL <- rnorm(length(universe),mean=7.9, sd=1.154)
#' genelength <- exp(logGL)
#' 
#' # Define a list of gene sets containing two GO terms
#' setname1 <- GOtoID$GO[1]
#' setname1
#' keep.set1 <- GOtoID$GO %in% setname1
#' set1 <- GOtoID$ENTREZID[keep.set1]
#' setname2 <- GOtoID$GO[2]
#' setname2
#' keep.set2 <- GOtoID$GO %in% setname2
#' set2 <- GOtoID$ENTREZID[keep.set2]
#' # Make the gene sets into a list
#' sets <- list(set1, set2)
#' names(sets) <- c(setname1,setname2)
#' 
#' # Test for enrichment of gene sets with no gene length bias
#' # The genes are randomly selected so we don't expect significant results
#' gsaseq(sig.de = de.genes, universe = universe, collection = sets)
#' 
#' # Test for enrichment of gene sets taking into account gene length bias
#' # Since the gene lengths are randomly generated this shouldn't make much 
#' # difference to the results
#' # Using log(gene length) or gene length doesn't make a difference to the 
#' # p-values because the probability weighting function is transformation
#' # invariant
#' gsaseq(sig.de = de.genes, univers = universe, collection = sets, 
#' gene.length = genelength)
#' }
#' 
gsaseq <- function(sig.de, universe, collection, plot.bias=FALSE, 
                   gene.length=NULL, sort = TRUE)
# Generalised version of goana with user-specified gene sets 
# Gene sets collections must be Entrez Gene ID
# Can take into account gene length bias
# Belinda Phipson
# 4 May 2020
{
    
    if(!is.vector(sig.de))
        stop("Input DE list is not a character vector")
    
    if(is.null(universe)){
        stop("Please supply the universe: all genes analysed in the experiment")
    }
    
    if(!is.null(gene.length)){
        if(length(gene.length)!=length(universe)){
            stop("Gene length and universe must be the same length")
        }
    }
    
    # Check collection is a list with character vectors
    if(!is.list(collection))
        collection <- list(collection=collection)
    collection <- lapply(collection, as.character)
    # Make sure gene set collections don't have any NAs
    collection <- lapply(collection, function(x) x[!is.na(x)])
    # Remove genes that are not in the universe from collections
    collection <- lapply(collection, function(x) x[x %in% universe])
    # Remove collections with no genes left after universe filter
    inUniv <- sapply(collection, function(x) length(x) > 0)
    collection <- collection[inUniv]
    
    test.de <- as.numeric(universe %in% sig.de)
    
    # Estimate prior probabilities
    if(!is.null(gene.length)){
        pwf <- .estimatePWF(D=test.de,bias=gene.length)
        if(plot.bias)
            .plotBiasSeq(D=test.de,bias=as.vector(gene.length))
        }
    
    results <- matrix(NA,ncol=4,nrow=length(collection))
    colnames(results) <- c("N","DE","P.DE","FDR")
    rownames(results) <- names(collection)
    results[,"N"] <- unlist(lapply(collection,length))
    results[,"DE"] <- unlist(lapply(collection, function(x) sum(sig.de %in% x)))
    Nuniverse <- length(universe)
    m <- length(sig.de)
    
    # Hypergeometric test with prior probabilities
    if(!is.null(gene.length)){
        for(i in 1:length(collection)){
            InSet <- universe %in% collection[[i]]
            pw.red <- sum(pwf[InSet])/results[i,"N"]
            pw.white <- sum(pwf[!InSet])/(Nuniverse-results[i,"N"])
            odds <- pw.red/pw.white
            results[i,"P.DE"] <- BiasedUrn::pWNCHypergeo(results[i,"DE"],
                                                         results[i,"N"],
                                                         Nuniverse-results[i,"N"],
                                                         m,
                                                         odds,
                                                         lower.tail=FALSE) + 
                BiasedUrn::dWNCHypergeo(results[i,"DE"],
                                        results[i,"N"],
                                        Nuniverse-results[i,"N"],
                                        m, odds)
        }
    }
    # Hypergeometric test without prior probabilities
    else{
        for(i in 1:length(collection)){
            results[i,"P.DE"] <- stats::phyper(q=results[i,"DE"]-0.5, m=m, 
                                               n=Nuniverse-m,
                                               k=results[i,"N"],
                                               lower.tail=FALSE)    
        }
    }
    results[,"FDR"] <- stats::p.adjust(results[,"P.DE"],method="BH")
    if(sort){
        o <- order(results[,"P.DE"])
        results[o,]
    }
    else
        data.frame(results)
}

.plotBiasSeq <- function(D,bias)
    # Plotting function to show gene level CpG density bias
    # Belinda Phipson
    # 4 March 2020
{
    o <- order(bias)
    splitf <- rep(1:100,each=200)[1:length(bias)]
    avgbias <- tapply(bias[o],factor(splitf),mean)
    sumDM <- tapply(D[o],factor(splitf),sum)
    propDM <- sumDM/table(splitf)
    graphics::par(mar=c(5,5,2,2))
    graphics::plot(avgbias,as.vector(propDM),
                   xlab="Gene length (binned)",
                   ylab="Proportion Differential Expression",cex.lab=1.5,
                   cex.axis=1.2)
    graphics::lines(stats::lowess(avgbias,propDM),col=4,lwd=2)
}
