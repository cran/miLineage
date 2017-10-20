\name{miLineage-package}
\alias{miLineage-package}
\alias{miLineage}
\docType{package}
\title{
Association Tests for Microbial Lineages on a Taxonomic Tree
}
\description{
A variety of association tests for microbiome data analysis described in Tang Z.-Z. et al.(2017) <doi: 10.1093/bioinformatics/btw804> and Tang Z.-Z. & Chen G. (2017) <Submitted>.
miLineage allows users to (a) perform tests on multivariate taxon counts; (b) localize the covariate-associated lineages on the taxonomic tree; and (c) assess the overall association of the microbial community with the covariate of interest. 
}
\details{
\tabular{ll}{
Package: \tab miLineage\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2017-10-18\cr
License: \tab GPL (>=2)\cr
}
QCAT, QCAT_GEE, ZIGDM
}
\author{
Zheng-Zheng Tang

Maintainer: Zheng-Zheng Tang <tang@biostat.wisc.edu>

}
\references{
1. Tang, Z.-Z., Chen, G., Alekseyenko, A. V., and Li, H. (2017). A General Framework for Association Analysis of Microbial Communities on a Taxonomic Tree. Bioinformatics, 33(9), 1278-1285.\cr
2. Tang, Z.-Z. and Chen, G. (2017). Zero-Inflated Generalized Dirichlet Multinomial Regression Model for Microbiome Compositional Data Analysis. Submitted.
}

\keyword{ package }
\examples{
data(data.toy)
OTU.toy = data.toy$OTU.toy
Tax.toy = data.toy$Tax.toy
case = data.toy$covariate.toy
# the OTUs should be consistent between the OTU table and the taxonomy table
OTU.toy.reorder = OTU.toy[,match(rownames(Tax.toy), colnames(OTU.toy))]
# perform QCAT test for detecting differential mean
QCAT(OTU.toy.reorder, case, 1, Tax.toy, fdr.alpha=0.05)
# perform ZIGDM test for detecting differential dispersion
ZIGDM(OTU.toy.reorder, NULL, NULL, case, test.type = "Disp", 
1, ZI.LB = 10, Tax.toy, fdr.alpha = 0.05)

}
