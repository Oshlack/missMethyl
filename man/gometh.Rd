\name{gometh}
\alias{gometh}
\title{
Gene ontology testing for Ilumina methylation array data
}
\description{
Tests gene ontology enrichment for significant CpGs from Illumina's Infinium HumanMethylation450 or MethylationEPIC array, taking into account the differing number of probes per gene present on the array. 
}
\usage{
gometh(sig.cpg, all.cpg = NULL, collection = "GO", array.type = c("450K","EPIC"), plot.bias = FALSE, prior.prob = TRUE)
}
\arguments{
  \item{sig.cpg}{
      character vector of significant CpG sites to test for GO term enrichment
  }
  \item{all.cpg}{
      character vector of all CpG sites tested. Defaults to all CpG sites on the array.
  }
  \item{collection}{
      the collection of pathways to test. Options are "GO" and "KEGG". Defaults to "GO".
  }
  \item(array.type){
      the Illumina methylation array used. Options are "450K" or "EPIC". Defaults to "450K".
  }
  \item{plot.bias}{
      logical, if true a plot showing the bias due to the differing numbers of probes per gene will be displayed
  }
  \item{prior.prob}{
      logical, if true will take into account the probability of significant differentially methylation due to numbers of probes per gene. If false, a hypergeometric test is performed ignoring any bias in the data.
  }
}
\details{
This function takes a character vector of significant CpG sites, maps the CpG sites to Entrez Gene IDs, and tests for GO term or KEGG pathway enrichment using a hypergeometric test, taking into account the number of CpG sites per gene on the 450K/EPIC array. 

Geeleher et al. (2013) showed that a severe bias exists when performing gene set analysis for genome-wide methylation data that occurs due to the differing numbers of CpG sites profiled for each gene. \code{gometh} is based on the \code{goseq} method (Young et al., 2010) and calls the \code{goana} function for GO testing, or the \code{kegga} function for KEGG testing, both of which are from the \code{limma} package (Ritchie et al. 2015). If \code{prior.prob} is set to FALSE, then prior probabilities are not used and it is assumed that each gene is equally likely to have a significant CpG site associated with it.

Genes associated with each CpG site are obtained from the annotation package \code{IlluminaHumanMethylation450kanno.ilmn12.hg19} if the array type is "450K". For the EPIC array, the annotation package \code{IlluminaHumanMethylationEPICanno.ilm10b2.hg19} is used. In order to get a list which contains the mapped Entrez gene IDS, please use the \code{getMappedEntrezIDs} function.

\code{gometh} tests all GO or KEGG terms, and false discovery rates are calculated using the method of Benjamini and Hochberg (1995). 

The \code{limma} functions \code{topGO} and \code{topKEGG} can be used to display the top 20 most enriched pathways. 

For more generalised gene set testing where the user can specify the gene set/s of interest to be tested, please use the \code{gsameth} function.
}
\value{
A data frame with a row for each GO or KEGG term and the following columns:
  \item{Term}{
      GO term if testing GO pathways
  }
  \item{Ont}{
      ontology that the GO term belongs to if testing GO pathways. "BP" - biological process, "CC" - cellular component, "MF" - molecular function.
  }
  \item{Pathway}{
      the KEGG pathway being tested if testing KEGG terms.
  }
  \item{N}{
      number of genes in the GO or KEGG term
  }
  \item{DE}{
      number of genes that are differentially methylated
  }
  \item{P.DE}{
      p-value for over-representation of the GO or KEGG term term
  }
  \item{FDR}{
      False discovery rate
  }
}
\references{
Phipson, B., Maksimovic, J., and Oshlack, A. (2016). missMethyl: an R package for analysing methylation data from Illuminas HumanMethylation450 platform. \emph{Bioinformatics}, \bold{15};32(2), 286--8.

Geeleher, P., Hartnett, L., Egan, L. J., Golden, A., Ali, R. A. R., and Seoighe, C. (2013). Gene-set analysis is severely biased when applied to genome-wide methylation data. \emph{Bioinformatics}, \bold{29}(15), 1851--1857.

Young, M. D., Wakefield, M. J., Smyth, G. K., and Oshlack, A. (2010). Gene ontology analysis for RNA-seq: accounting for selection bias. \emph{Genome Biology}, 11, R14.

Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., and Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. \emph{Nucleic Acids Research}, gkv007.

Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. \emph{Journal of the Royal Statistical Society Series}, B, \bold{57}, 289-300.
}
\seealso{
\code{\link{goana},\link{kegga},\link{gsameth}}
}
\author{
Belinda Phipson
}
\examples{
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Randomly select 1000 CpGs to be significantly differentially methylated
sigcpgs <- sample(rownames(ann),1000,replace=FALSE)

# All CpG sites tested
allcpgs <- rownames(ann)

# GO testing with prior probabilities taken into account
# Plot of bias due to differing numbers of CpG sites per gene
gst <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "GO", plot.bias = TRUE, prior.prob = TRUE)

# Total number of GO categories significant at 5\% FDR
table(gst$FDR<0.05)

# Table of top GO results
topGO(gst)

# GO testing ignoring bias
gst.bias <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "GO", prior.prob=FALSE)

# Total number of GO categories significant at 5\% FDR ignoring bias
table(gst.bias$FDR<0.05)

# Table of top GO results ignoring bias
topGO(gst.bias)

# KEGG testing
kegg <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "KEGG", prior.prob=TRUE)

# Table of top KEGG results
topKEGG(kegg)
}
