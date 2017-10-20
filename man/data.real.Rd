\name{data.real}
\alias{data.real}
\docType{data}
\title{
Real data that include an OTU table, a taxonomy table, and a covariate table
}
\description{
The data are derived from a real gut microbiome study. The data contain an OTU table, a taxonomy table, and a covariate table. The data include 96 samples and 80 OTUs.   
}
\usage{data("data.real")}
\format{
data.real contains the following objects:
\itemize{
\item OTU.real: a matrix of OTU counts for 96 samples and 80 OTUs
\item Tax.real: a matrix of taxonomy table from Rank1 (kingdom level) to Rank6 (genus level)  
\item covariate.real: a matrix of three variables
}
}
\source{
Wu, Gary D., et al. "Linking long-term dietary patterns with gut microbial enterotypes." Science 334.6052 (2011): 105-108.
}
\examples{
data(data.real)
}
\keyword{datasets}
