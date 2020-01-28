\name{missMethylnews}
\title{missMethyl News}
\encoding{UTF-8}

\section{Version 1.20.4}{\itemize{

\item Bug fix: Fixed bug that resulted in negative value for odds. Bug only occurred in rare cases when CpG universe was smaller than all CpGs on array AND if some large collections had more genes than the universe.    
}}

\section{Version 1.20.3}{\itemize{

\item Added \code{fract.counts} logical parameter to enable fractional counting of CpGs 
to be turned off in gometh, gsameth, goregion and gsaregion. This switches off 
the adjustment for CpGs that map to multiple genes. This is primarily intended 
for comparison with earlier versions and we strongly encourage it to be set to 
TRUE.
}}

\section{Version 1.20.0}{\itemize{

\item Updated gsameth to take into account that some CpGs map to multiple genes. For a small number of gene families, this previously caused their associated GO categories/gene sets to be erroneously overrepresented and thus highly significant.    
\item Updated gometh to call gsameth to perform statistical testing for GO categories and KEGG pathways; goana and kegga are no longer used.
\item Replaced use of the AnnotationDbi 'toTable' function with 'select' to speed up execution.
\item The 'topGSA' function is now used to look at the top ranked results from both gometh and gsameth. The vignette has been modified accordingly. 
}}

\section{Version 1.15.2}{\itemize{

\item Updated getMappedEntrezIDs, gometh and gsameth to to speed up execution by taking the array annotation in as an optional argument.
\item missMethyl now uses the latest IlluminaHumanMethylationEPICanno.ilm10b2.hg19 annotation by default for EPIC arrays.
}}


\section{Version 1.15.1}{\itemize{

\item Added getAdjusted function for extracting RUVm adjusted data for visualisation purposes
\item Updated vignette to demonstrate use of getAdjusted function
\item Vignette now includes an example of how to handle cases with RUVm where number of samples is greater than number of Illumina negative controls

}}

\section{Version 1.7.3}{\itemize{

\item modified gene set testing functions gometh and gsameth to accommodate EPIC arrays

}}

\section{Version 1.5.1}{\itemize{

\item New gene set testing function gsameth() added, kegg testing functionality added to gometh().

}}

\section{Version 1.1.10}{\itemize{

\item This version of the package contains functions to perform SWAN normalisation, RUV normalisation and differential methylation analysis, differential variability analysis and gene ontology testing for the 450K array.

}}

\section{Version 1.1.3}{\itemize{

\item
Added the gometh function to perform gene ontology analysis.

}}

\section{Version 1.1.2}{\itemize{

\item
Added functions to perform RUV normalization and differential methylation analysis.

}}

\section{Version 0.99.9}{\itemize{

\item
This version of the package has all calculations in matrix formulation and handles nuisance parameters appropriately.

}}

\section{Version 0.99.7}{\itemize{

\item
Added new function, contrasts.varFit. 

}}


\section{Version 0.99.6}{\itemize{

\item
Modified getLeveneResiduals 

}}

\section{Version 0.99.5}{\itemize{

\item
Modified the varFit function to accept DGEList objects for differential variability testing for RNA-Seq data.

}}

\section{Version 0.99.4}{\itemize{

\item
First version of package contains functions to perform SWAN normalisation and differential variability analysis for DNA methylation data from Illumina's Infinium HumanMethylation450 beadchip.

}}