\name{multibitTree-package} 
\alias{multibitTree-package}
\alias{multibitTree}
\docType{package}
\title{
Functions for Building a MulitbitTree and Searching
}
\description{
This package contains functions for loading sets of binary vectors
(fingerprints) into a search tree and perform fast similarity matching.

The measure of similarity is calculated by the Tanimoto-Coefficient.

The search tree's data structure and algorithms are based on a
combination of a 1DGrid, Multibit trees and XOR-folding
as discribed in the article referenced below.
}
\details{
\tabular{ll}{
Package: \tab multibitTree\cr
Type: \tab Package\cr
Version: \tab 1.7\cr
Date: \tab 2015-01-08\cr
License: \tab GPL-3\cr
}
}
\author{
Markus Koetter, iQbis consulting GmbH

Maintainer: Prof. Dr. Rainer Schnell <rainer.schnell@uni-due.de>
}
\references{
Kristensen et al. Algorithms for Molecular Biology 2010:
A tree-based method for the rapid screening of
chemical fingerprints, 
http://www.almob.org/content/5/1/9
}
\keyword{misc}
