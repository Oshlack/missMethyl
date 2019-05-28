\name{getINCs}
\alias{getINCs}
\title{
Extract intensity data for Illumina negative controls found on 450k or EPIC 
arrays.
}
\description{
Extracts the intensity data for the Illumina negative controls found on 450k 
or EPIC arrays and returns a matrix of M-values (log2 ratio of the green to red 
intensities).   
}
\usage{
getINCs(rgSet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rgSet}{
An object of class \code{RGChannelSet}.
}
}
\details{
The getINCs function extracts the intensity data for the INCs from an \linkS4class{RGChannelSet} object. The function retrieves both the green and red channel intensity values and returns the data as the log2 ratio of the green to red intensities. Essentially, the INCs are being treated like Type II probes for which the M-values are also given as the log2 ratio of the green to red intensities.
}
\value{
An matrix of M-values.
}
\author{
Jovana Maksimovic \email{jovana.maksimovic@mcri.edu.au}
}
\seealso{
\code{\linkS4class{RGChannelSet}}
}
\examples{
if (require(minfi) & require(minfiData)) {

  INCs <- getINCs(RGsetEx)
  head(INCs)
  dim(INCs)
}

}

