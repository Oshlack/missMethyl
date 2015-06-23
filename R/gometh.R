gometh <- function(sig.cpg, all.cpg=NULL, plot.bias=FALSE, prior.prob=TRUE)
# Gene ontology testing for 450K methylation arrays based on goseq
# Takes into account probability of differential methylation based on
# numbers of probes on array per gene
# Belinda Phipson
# 28 January 2015
{
    # check input
    sig.cpg <- as.character(sig.cpg)
    sig.cpg <- sig.cpg[!is.na(sig.cpg)]
    
    # flatten 450k array annotation
    ann.450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    
    # get rid of the non-CpG sites
    strlen<-str_length(rownames(ann.450k))
    ann.keep<-ann.450k[strlen==10,]
    
    # get rid of CpGs that are not annotated
    missing<-ann.keep$UCSC_RefGene_Name==""
    ann.keep<-ann.keep[!missing,]
    
    # get individual gene names for each CpG
    geneslist<-strsplit(ann.keep$UCSC_RefGene_Name,split=";")
    names(geneslist)<-rownames(ann.keep)
    
    grouplist<-strsplit(ann.keep$UCSC_RefGene_Group,split=";")
    names(grouplist)<-rownames(ann.keep)
    
    flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))
    flat$symbol<-as.character(flat$symbol)
    flat$group <- as.character(flat$group)
    
    flat$cpg<- substr(rownames(flat),1,10)
        
    flat$alias <- alias2SymbolTable(flat$symbol)
    
    eg <- toTable(org.Hs.egSYMBOL2EG)
    m <- match(flat$alias,eg$symbol)
    flat$entrezid <- eg$gene_id[m]
    flat <- flat[!is.na(flat$entrezid),]
    
    # keep unique cpg by gene name annotation
    id<-paste(flat$cpg,flat$entrezid,sep=".")
    d <- duplicated(id)
    flat.u <- flat[!d,]
    
    if(is.null(all.cpg))
        all.cpg <- unique(flat.u$cpg)
    else{
        all.cpg <- as.character(all.cpg)
        all.cpg <- all.cpg[!is.na(all.cpg)]
        all.cpg <- unique(all.cpg)    
    }
    
    # map CpG sites to entrez gene id's
    sig.cpg <- unique(sig.cpg)
    m1 <- match(flat.u$cpg,sig.cpg)
    eg.sig <- flat.u$entrezid[!is.na(m1)]
    eg.sig <- unique(eg.sig)
    
    m2 <- match(flat.u$cpg,all.cpg)
    eg.all <- flat.u$entrezid[!is.na(m2)]
        
    freq_genes <- table(eg.all)
    eg.universe <- names(freq_genes)
        
    test.de <- as.integer(eg.universe %in% eg.sig)
    
    sorted.eg.sig <- eg.universe[test.de==1]
    
    # get gene-wise prior probabilities
    if(prior.prob){
        pwf <- .estimatePWF(D=test.de,bias=as.vector(freq_genes))
        if(plot.bias)
        	.plotBias(D=test.de,bias=as.vector(freq_genes))
    # go testing
        gst.go <- goana(sorted.eg.sig,universe=eg.universe,prior.prob=pwf)
    }
    else
        gst.go <- goana(sorted.eg.sig,universe=eg.universe)
    
    gst.go$FDR<-p.adjust(gst.go$P.DE,method="BH")
    gst.go
}

.plotBias <- function(D,bias)
# Plotting function to show gene level CpG density bias
# Belinda Phipson
# 5 March 2015
{
    o <- order(bias)
    splitf <- rep(1:100,each=200)[1:length(bias)]
    avgbias <- tapply(bias[o],factor(splitf),mean)
    sumDM <- tapply(D[o],factor(splitf),sum)
    propDM <- sumDM/table(splitf)
    par(mar=c(5,5,2,2))
    plot(avgbias,as.vector(propDM),xlab="Number of CpGs per gene", ylab="Proportion Differential Methylation",cex.lab=1.5,cex.axis=1.2)
    lines(lowess(avgbias,propDM),col=4,lwd=2)
}

.estimatePWF <- function(D,bias)
# An alternative to goseq function nullp, which is transformation invariant
# Belinda Phipson and Gordon Smyth
# 6 March 2015
{
    prior.prob <- bias
    o <- order(bias)
    prior.prob[o] <- tricubeMovingAverage(D[o],span=0.5,full.length=TRUE)
    prior.prob
}
