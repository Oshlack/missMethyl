\name{contrasts.varFit}
\alias{contrasts.varFit}
\title{
	Compute contrasts for a varFit object.
}
\description{
	Compute estimated coefficients, standard errors	and LogVarRatios for a given set of contrasts.
}
\usage{
contrasts.varFit(fit, contrasts=NULL)
}
\arguments{
  \item{fit}{
	list containing a linear model fit produced by \code{varFit}. The fit object should be of class \code{MArrayLM}.
}
  \item{contrasts}{
	numeric matrix with rows corresponding to coefficients in \code{fit} and columns containing contrasts. 
}
}
\details{
	This function calls the \code{contrasts.fit} function in \code{limma} to compute coefficients and standard errors for the specified contrasts corresponding to a linear model fit obtained from the \code{varFit} function. LogVarRatios are also computed in terms of the contrasts. A contrasts matrix can be computed using the \code{makeContrasts} function.
}
\value{
	A list object of the same class as \code{fit}.
}
\author{
	Belinda Phipson
}
\seealso{
	\code{varFit}, \code{contrasts.fit}, \code{makeContrasts}
}
\examples{
# Randomly generate data for a 3 group problem with 100 CpG sites and 4 arrays in each group. 

library(limma)

y<-matrix(rnorm(1200),ncol=12)

group<-factor(rep(c(1,2,3),each=4))
design<-model.matrix(~0+group)
colnames(design)<-c("grp1","grp2","grp3")

# Fit linear model for differential variability
vfit<-varFit(y,design)

# Specify contrasts
contr<-makeContrasts(grp2-grp1,grp3-grp1,grp3-grp2,levels=colnames(design))

# Compute contrasts from fit object
vfit.contr<-contrasts.varFit(vfit,contrasts=contr)

summary(decideTests(vfit.contr))

# Look at top table of results for first contrast

topVar(vfit.contr,coef=1)

}

