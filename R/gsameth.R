gsameth <- function(sig.cpg, all.cpg=NULL, collection, array.type = c("450K","EPIC"), 
                    plot.bias=FALSE, prior.prob=TRUE, anno=NULL)
# Generalised version of gometh with user-specified gene sets 
# Gene sets collections must be Entrez Gene ID
# Can take into account probability of differential methylation 
# based on numbers of probes on array per gene.
# Belinda Phipson
# 10 February 2016
# Updated 18 September 2018
{
    
    if(!is.vector(sig.cpg))
        stop("Input CpG list is not a character vector")
    array.type <- match.arg(toupper(array.type),c("450K","EPIC"))
     
    # Get mapped entrez gene IDs from CpG probe names
    if(!is.null(anno)){
      out <- getMappedEntrezIDs(sig.cpg=sig.cpg,all.cpg=all.cpg,array.type=array.type,
                                anno)
    } else {
      out <- getMappedEntrezIDs(sig.cpg=sig.cpg,all.cpg=all.cpg,array.type=array.type)
    }
    sorted.eg.sig <- out$sig.eg
    eg.universe <- out$universe
    freq_genes <- out$freq
    test.de <- out$de
    
    # Check collection is a list with character vectors
    if(!is.list(collection))
            collection <- list(collection=collection)
    collection <- lapply(collection, as.character)
    # Make sure gene set collections don't have any NAs
    collection <- lapply(collection, function(x) x[!is.na(x)])
    # Make sure only collections with geneids present in universe are included
    inUniv <- lapply(collection, function(x) sum(eg.universe %in% x))
    collection <- collection[inUniv!=0]

    # Estimate prior probabilities
    if(prior.prob){
        pwf <- .estimatePWF(D=test.de,bias=as.vector(freq_genes))
        if(plot.bias)
            .plotBias(D=test.de,bias=as.vector(freq_genes))
    }
        
    results <- matrix(NA,ncol=4,nrow=length(collection))
    colnames(results) <- c("N","DE","P.DE","FDR")
    rownames(results) <- names(collection)
    results[,"N"] <- unlist(lapply(collection,length))
    results[,"DE"] <- unlist(lapply(collection, function(x) sum(sorted.eg.sig %in% x)))
    Nuniverse <- length(eg.universe)
    m <- length(sorted.eg.sig)
        
    # Hypergeometric test with prior probabilities
    if(prior.prob){
        for(i in 1:length(collection)){
            InSet <- eg.universe %in% collection[[i]]
            pw.red <- sum(pwf[InSet])/results[i,"N"]
            pw.white <- sum(pwf[!InSet])/(Nuniverse-results[i,"N"])
            odds <- pw.red/pw.white
            results[i,"P.DE"] <- BiasedUrn::pWNCHypergeo(results[i,"DE"],results[i,"N"],Nuniverse-results[i,"N"],m,odds,lower.tail=FALSE,precision=1E-32) 
        }
    }
    # Hypergeometric test without prior probabilities
    else{
        for(i in 1:length(collection)){
            results[i,"P.DE"] <- phyper(q=results[i,"DE"]-0.5,m=m,n=Nuniverse-m,k=results[i,"N"],lower.tail=FALSE)    
        }
    }
    results[,"FDR"] <- p.adjust(results[,"P.DE"],method="BH")
    results
}

topGSA <- function(gsa,number=20,sort=TRUE)
# Function to output the top 20 enriched gene sets from a gene set 
# enrichment analysis using gsameth.
# Belinda Phipson
# 12 February
{
    if(nrow(gsa) < number)
        number <- nrow(gsa)
    if(sort){
        o <- order(gsa[,"P.DE"])
        out <- gsa[o,]
        out[1:number,]
    }
    else
        gsa[1:number,]
}