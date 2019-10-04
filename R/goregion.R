goregion <- function(regions, all.cpg=NULL, collection=c("GO","KEGG"), array.type = c("450K","EPIC"),
                     plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg = TRUE)
  # Gene ontology testing or KEGG pathway analysis for differentially methylated regions
  # Takes into account probability of differential methylation based on
  # numbers of probes on array per gene
  # Jovana Maksimovic
  # 26 April 2019. Last updated 26 April 2019.
{
  collection <- match.arg(toupper(collection),c("GO","KEGG"))
  
  if(is.null(anno)){
    if(array.type=="450K"){
      anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    } else {
      anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }
  }
  
  if(!is.null(all.cpg)){
    anno <- anno[all.cpg,]
  }
  
  cpgs <- GRanges(seqnames = anno$chr, 
                  ranges = IRanges(start = anno$pos, 
                                   end = anno$pos),
                  strand = anno$strand,
                  name = anno$Name)
  
  overlaps <- findOverlaps(cpgs,regions)
  sig.cpg <- cpgs$name[from(overlaps)]
  
  result <- gometh(sig.cpg=sig.cpg, all.cpg=all.cpg, collection=collection, 
                      array.type=array.type, plot.bias=plot.bias, prior.prob=prior.prob, 
                      anno=anno, equiv.cpg=equiv.cpg)
  result
}  

gsaregion <- function(regions, all.cpg=NULL, collection, array.type = c("450K","EPIC"),
                   plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg = TRUE)
  # Generalised version of goregion with user-specified gene sets 
  # Gene sets collections must be Entrez Gene ID
  # Takes into account probability of differential methylation based on
  # numbers of probes on array per gene
  # Jovana Maksimovic
  # 26 April 2019. Last updated 26 April 2019.
{
  if(is.null(anno)){
    if(array.type=="450K"){
      anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    } else {
      anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }
  }
  
  if(!is.null(all.cpg)){
    anno <- anno[all.cpg,]
  }
  
  cpgs <- GRanges(seqnames = anno$chr, 
                  ranges = IRanges(start = anno$pos, 
                                   end = anno$pos),
                  strand = anno$strand,
                  name = anno$Name)
  
  overlaps <- findOverlaps(cpgs,regions)
  sig.cpg <- cpgs$name[from(overlaps)]
  
  result <- gsameth(sig.cpg=sig.cpg, all.cpg=all.cpg, collection=collection, 
                   array.type=array.type, plot.bias=plot.bias, prior.prob=prior.prob, 
                   anno=anno, equiv.cpg=equiv.cpg)
  result
}