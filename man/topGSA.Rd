\name{topGSA}
\alias{topGSA}
\title{
Get table of top 20 enriched pathways 
}
\description{
After using \code{gsameth}, calling topGSA will output the top 20 most significantly enriched pathways.
}
\usage{
topGSA(gsa, number = 20, sort = TRUE)
}
\arguments{
  \item{gsa}{
      matrix, from output of \code{gsameth}
  }
  \item{number}{
      scalar, number of pathway results to output. Default is 20
  }
  \item{sort}{
      logical, should the table be ordered by p-value. Default is TRUE.
  }      
}
\details{
    This function will output the top 20 most significant pathways from a pathway analysis using the \code{gsameth} function. The output is ordered by p-value.
}
\value{
A matrix ordered by P.DE, with a row for each gene set and the following columns:
  \item{N}{
  	number of genes in the gene set
  }
  \item{DE}{
  	number of genes that are differentially methylated
  }
  \item{P.DE}{
  	p-value for over-representation of the gene set
  }
  \item{FDR}{
  	False discovery rate, calculated using the method of Benjamini and Hochberg (1995).
  }
}
\seealso{
\code{\link{gsameth}}
}
\author{
Belinda Phipson
}
\examples{
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(org.Hs.eg.db)
library(limma)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Randomly select 1000 CpGs to be significantly differentially methylated
sigcpgs <- sample(rownames(ann),1000,replace=FALSE)

# All CpG sites tested
allcpgs <- rownames(ann)

# Use org.Hs.eg.db to extract a GO term
GOtoID <- toTable(org.Hs.egGO2EG)
setname1 <- GOtoID$go_id[1]
setname1
keep.set1 <- GOtoID$go_id %in% setname1
set1 <- GOtoID$gene_id[keep.set1]
setname2 <- GOtoID$go_id[2]
setname2
keep.set2 <- GOtoID$go_id %in% setname2
set2 <- GOtoID$gene_id[keep.set2]

# Make the gene sets into a list
sets <- list(set1, set2)
names(sets) <- c(setname1,setname2)

# Testing with prior probabilities taken into account
# Plot of bias due to differing numbers of CpG sites per gene
gst <- gsameth(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = sets, plot.bias = TRUE, prior.prob = TRUE)
topGSA(gst)

# Testing ignoring bias
gst.bias <- gsameth(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = sets, prior.prob = FALSE)
topGSA(gst.bias)
}
