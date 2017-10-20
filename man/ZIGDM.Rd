\name{ZIGDM}
\alias{ZIGDM}

\title{
Zero-Inflated Generalized Dirichlet Multinomial (ZIGDM) Tests
}
\description{
Different from the distribution-free QCAT and QCAT+GEE, the ZIGDM tests are based on a parametric model (ZIGDM) for multivariate taxon counts. 
The ZIGDM tests can not only detect differential mean but also differential dispersion level or presence-absence frequency in microbial compositions.
This function allows users to (a) perform ZIGDM tests on multivariate taxon counts; (b) perform ZIGDM tests on the taxonomic tree to localize the covariate-associated lineages; and (c) assess the overall association of the microbial community with the covariate of interest.
}
\usage{
ZIGDM(OTU, X4freq, X4mean, X4disp, test.type = "Mean", X.index, ZI.LB = 10, Tax = NULL, 
min.depth = 0, n.resample = NULL, fdr.alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OTU}{
a matrix contains counts with each row corresponds to a sample and each column corresponds to an OTU or a taxon. Column name is mandatory.
}
  \item{X4freq}{
a matrix contains covariates that link to presence-absence frequency in microbial compositions. Each column pertains to one variable. Set X4freq=NULL if only intercept term is needed. Samples in the OTU and X4freq matrices should be identical and in the same order.  
}
  \item{X4mean}{
a matrix contains covariates that link to mean abundance in microbial compositions. Each column pertains to one variable. Set X4mean=NULL if only intercept term is needed. Samples in the OTU and X4mean matrices should be identical and in the same order.  
}
  \item{X4disp}{
a matrix contains covariates that link to dispersion level in microbial compositions. Each column pertains to one variable. Set X4disp=NULL if only intercept term is needed. Samples in the OTU and X4disp matrices should be identical and in the same order.  
}
  \item{test.type}{
  If test.type = "Mean", the function will test for differential mean (Default). If test.type = "Disp", the function will test for differential dispersion. If test.type = "Freq", the function will test for differential presence-absence frequency. 
}
  \item{X.index}{
If test.type = "Mean", X.index is a vector indicates the columns in X4mean for the covariate(s) of interest.
If test.type = "Disp", X.index is a vector indicates the columns in X4disp for the covariate(s) of interest.
If test.type = "Freq", X.index is a vector indicates the columns in X4freq for the covariate(s) of interest.
}
  \item{ZI.LB}{
lower bound of zero counts for the taxon that needs zero-inflated model. 
The counts for a taxon is assumed to be zero-inflated if the number of zero observations for the taxon is greater than ZI.LB (Default is 10). If ZI.LB=NULL, GDM model (i.e. non-zero-inflated version of ZIGDM) will be applied.
}
  \item{Tax}{
a matrix define the taxonomic ranks with each row corresponds to an OTU or a taxon and each column corresponds to a rank (start from the higher taxonomic level, e.g., from kingdom to genus). Row name is mandatory and should be consistent with the column name of the OTU table,  Column name should be formatted as "Rank1", "Rank2", etc, representing the taxonomic levels from highest to lowest. 
}
  \item{min.depth}{
lower bound of sample read depth. Samples with read depth less than min.depth will be removed before the analysis.
}
  \item{n.resample}{
number of reamplings/permutations
}
  \item{fdr.alpha}{
desired false discovery rate for multiple tests on the lineages.
}
}
\value{
If Tax=NULL (Default), a test is performed using all the OTUs/taxa.

If Tax is provided, tests are performed for lineages derived from the taxonomic hierarchy. The output is a list that contains 3 components
\item{lineage.pval}{p-values for all lineages. By default (n.resample=NULL), only the asymptotic test will be performed. }
\item{sig.lineage}{a vector of significant lineages}
\item{global.pval}{p-values of the global tests}

}
\references{
Tang, Z.-Z. and Chen, G. (2017). Zero-Inflated Generalized Dirichlet Multinomial Regression Model for Microbiome Compositional Data Analysis. Submitted.

}
\author{
Zheng-Zheng Tang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(data.toy)
OTU.toy = data.toy$OTU.toy
Tax.toy = data.toy$Tax.toy
case = data.toy$covariate.toy
# the OTUs should be consistent between the OTU table and the taxonomy table
OTU.toy.reorder = OTU.toy[,match(rownames(Tax.toy), colnames(OTU.toy))]
# perform ZIGDM test for detecting differential dispersion
ZIGDM(OTU.toy.reorder, NULL, NULL, case, test.type = "Disp", 
1, ZI.LB = 10, Tax.toy, fdr.alpha = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ parametric association test }
