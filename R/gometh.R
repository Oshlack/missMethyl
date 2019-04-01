gometh <- function(sig.cpg, all.cpg=NULL, collection=c("GO","KEGG"), array.type = c("450K","EPIC"),
                   plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg = TRUE)
# Gene ontology testing or KEGG pathway analysis for Illumina methylation arrays based on goseq
# Takes into account probability of differential methylation based on
# numbers of probes on array per gene
# Belinda Phipson
# 28 January 2015. Last updated 29 March 2019.
# EPIC functionality contributed by Andrew Y.F. Li Yim
{
    if(!is.vector(sig.cpg))
        stop("Input CpG list is not a character vector")
    array.type <- match.arg(toupper(array.type),c("450K","EPIC"))
    collection <- match.arg(toupper(collection),c("GO","KEGG"))

    # Get mapped entrez gene IDs from CpG probe names
    if(!is.null(anno)){
      out <- getMappedEntrezIDs(sig.cpg=sig.cpg,all.cpg=all.cpg,array.type=array.type,
                                anno=anno)
    } else {
      out <- getMappedEntrezIDs(sig.cpg=sig.cpg,all.cpg=all.cpg,array.type=array.type)
    }
    sorted.eg.sig <- out$sig.eg
    eg.universe <- out$universe
    freq_genes <- out$freq
    test.de <- out$de
    equiv <- out$equiv

    # get gene-wise prior probabilities and perform testing
    if(prior.prob){
      if(equiv.cpg){
        pwf <- .estimatePWF(D=test.de,bias=as.vector(equiv))
        if(plot.bias)
            .plotBias(D=test.de,bias=as.vector(equiv))
        if(collection=="GO")
            gst <- limma:::goana(sorted.eg.sig,universe=eg.universe,prior.prob=pwf)
        if(collection=="KEGG")
            gst <- limma:::kegga(sorted.eg.sig,universe=eg.universe,prior.prob=pwf)
      }
      else{
        pwf <- .estimatePWF(D=test.de,bias=as.vector(freq_genes))
        if(plot.bias)
          .plotBias(D=test.de,bias=as.vector(freq_genes))
        if(collection=="GO")
          gst <- limma:::goana(sorted.eg.sig,universe=eg.universe,prior.prob=pwf)
        if(collection=="KEGG")
          gst <- limma:::kegga(sorted.eg.sig,universe=eg.universe,prior.prob=pwf)
      }
    }
    # Perform GO or KEGG testing without correcting for CpG density bias
    else{
        if(collection=="GO")
            gst <- limma:::goana(sorted.eg.sig,universe=eg.universe)
        if(collection=="KEGG")
            gst <- limma:::kegga(sorted.eg.sig,universe=eg.universe)
    }
    gst$FDR<-p.adjust(gst$P.DE,method="BH")
    gst
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
    plot(avgbias,as.vector(propDM),xlab="Number of CpGs per gene (binned)",
         ylab="Proportion Differential Methylation",cex.lab=1.5,cex.axis=1.2)
    lines(lowess(avgbias,propDM),col=4,lwd=2)
}

.estimatePWF <- function(D,bias)
# An alternative to goseq function nullp, which is transformation invariant
# Belinda Phipson and Gordon Smyth
# 6 March 2015
{
    prior.prob <- bias
    o <- order(bias)
    prior.prob[o] <- limma:::tricubeMovingAverage(D[o],span=0.5)
    prior.prob
}

.getFlatAnnotation <- function(array.type=c("450K","EPIC"),anno=NULL)
  # flatten 450k or EPIC array annotation
  # Jovana Maksimovic
  # 18 September 2018
  # Updated 18 September 2018
  # Modified version of Belida Phipson's .flattenAnn code
{
  if(is.null(anno)){
    if(array.type=="450K"){
      anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    } else {
      anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }
  }

  # get rid of the non-CpG sites
  ann.keep<-anno[grepl("^cg",anno$Name),]

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

  #flat$cpg <- rownames(flat)
  flat$alias <- suppressWarnings(limma:::alias2SymbolTable(flat$symbol))

  eg <- toTable(org.Hs.egSYMBOL2EG)
  m <- match(flat$alias,eg$symbol)
  flat$entrezid <- eg$gene_id[m]
  flat <- flat[!is.na(flat$entrezid),]

  # keep unique cpg by gene name annotation
  id<-paste(flat$cpg,flat$entrezid,sep=".")
  d <- duplicated(id)
  flat.u <- flat[!d,]
  flat.u
  # This randomly samples only 1 gene ID for multimapping CpGs
  #.reduceMultiMap(flat.u)
}

.reduceMultiMap <- function(flat){
  mm <- table(flat$cpg)
  mm <- names(mm[mm > 2])
  red <- tapply(flat$entrezid[flat$cpg %in% mm],
                 flat$cpg[flat$cpg %in% mm], sample, 1)
  key <- paste(flat$cpg,flat$entrezid,sep=".")
  rkey <- paste(names(red),red,sep = ".")
  w <- which(key %in% rkey)
  
  newmm <- flat[w,] 
  newflat <- flat[-which(flat$cpg %in% mm),]
  newflat <- rbind(newflat, newmm)
  newflat
}

getMappedEntrezIDs <- function(sig.cpg,all.cpg=NULL,array.type,anno=NULL)
  # From a list of CpG sites, obtain the Entrez Gene IDs that are used for testing pathway enrichment
  # Belinda Phipson
  # 10 February 2016
  # Updated 29 March 2019
{
  # check input
  sig.cpg <- as.character(sig.cpg)
  sig.cpg <- sig.cpg[!is.na(sig.cpg)]

  # Get annotaton in appropriate format
  if(is.null(anno)){
    flat.u <- .getFlatAnnotation(array.type)
  } else {
    flat.u <- .getFlatAnnotation(array.type,anno)
  }

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
  
  multimap <- data.frame(table(flat.u$cpg))
  multimap$Var1 <- as.character(multimap$Var1)
  m3 <- match(flat.u$cpg, multimap$Var1)
  flat.u$multimap <- multimap$Freq[m3]
  
  flat.u$inv.multimap <- 1/flat.u$multimap
  
  equivN <- tapply(flat.u$inv.multimap,flat.u$entrezid,sum)
  mm <- match(eg.universe,names(equivN))
  equivN <- equivN[mm]
  
  sig.flat <- flat.u[!is.na(m1),]
  
  fract <- data.frame(weight=pmin(tapply(1/sig.flat$multimap,sig.flat$entrezid,sum),1))
  
  m4 <- match(sorted.eg.sig,rownames(fract))
  fract.counts <- fract$weight[m4]
  
  out <- list(sig.eg = sorted.eg.sig, universe = eg.universe, 
              freq = freq_genes, equiv =  equivN, de = test.de, fract.counts = data.frame(sigid=sorted.eg.sig,frac=fract.counts))
  out
}

