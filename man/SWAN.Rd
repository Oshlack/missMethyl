\name{SWAN}
\alias{SWAN}
\alias{SWAN.default}
\alias{SWAN.MethyLumiSet}
\alias{SWAN.RGChannelSet}
\title{
Subset-quantile Within Array Normalisation for Illumina Infinium HumanMethylation450 BeadChips
}
\description{
Subset-quantile Within Array Normalisation (SWAN) is a within array normalisation method for the Illumina Infinium HumanMethylation450 platform. It allows Infinium I and II type probes on a single array to be normalized together.
}
\usage{
SWAN(data, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An object of class either \code{MethylSet}, \code{RGChannelSet} or \code{MethyLumiSet}.
}
  \item{verbose}{
Should the function be verbose?}
}
\details{
The SWAN method has two parts. First, an average quantile distribution is created using a subset of probes defined to be biologically similar based on the number of CpGs underlying the probe body. This is achieved by randomly selecting N Infinium I and II probes that have 1, 2 and 3 underlying CpGs, where N is the minimum number of probes in the 6 sets of Infinium I and II probes with 1, 2 or 3 probe body CpGs. If no probes have previously been filtered out e.g. sex chromosome probes, etc. N=11,303. This results in a pool of 3N Infinium I and 3N Infinium II probes. The subset for each probe type is then sorted by increasing intensity. The value of each of the 3N pairs of observations is subsequently assigned to be the mean intensity of the two probe types for that row or 'quantile'. This is the standard quantile procedure. The intensities of the remaining probes are then separately adjusted for each probe type using linear interpolation between the subset probes.
}
\value{
An object of class \code{MethylSet}
}
\references{
J Maksimovic, L Gordon and A Oshlack (2012).  \emph{SWAN: Subset quantile Within-Array Normalization for Illumina Infinium HumanMethylation450 BeadChips}.  Genome Biology 13, R44.
}
\author{
Jovana Maksimovic \email{jovana.maksimovic@mcri.edu.au}
}
\note{
SWAN uses a random subset of probes to perform the within-array normalization. In order to achive reproducible results, the seed needs to be set using \code{set.seed}.
}

\seealso{
\code{\linkS4class{RGChannelSet}} and
\code{\linkS4class{MethylSet}} as well as
\code{\linkS4class{MethyLumiSet}} and
\code{\linkS4class{IlluminaMethylationManifest}}.
}

\examples{
if (require(minfi) & require(minfiData)) {

  set.seed(100)
  datSwan1 <- SWAN(RGsetEx)
  
  dat <- preprocessRaw(RGsetEx)
  set.seed(100)
  datSwan2 <- SWAN(dat)
  
  head(getMeth(datSwan2)) == head(getMeth(datSwan1))
}

}

