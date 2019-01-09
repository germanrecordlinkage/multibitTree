\name{multibitTree.search}
\alias{multibitTree.search}
\title{
Search Fingerprints in MultibitTree
}
\description{
This function searches in a loaded MultibitTree. With a given query fingerprint and
Tanimoto coefficient the matching fingerprints will be returned.
}
\usage{
multibitTree.search(query, minTanimoto, size = 0, sort = FALSE)
}
\arguments{
  \item{query}{
  a character string consisting of the characters "0" and "1" representing a fingerprint to search for
}
  \item{minTanimoto}{
  a numeric value giving the lower bound of tanimoto coefficient to search for
}
  \item{size}{
  number of fingerprints that shall be returned in maximum
}
  \item{sort}{
  logical flag if the result shall be sorted, starting with the highest Tanimoto coefficient
}
}
\value{
The function returns a data.frame with two columns:
\item{fingerprint}{
  this column contains the matching fingerprint id
}
\item{tanimoto}{
  this column contains the corresponding Tanimoto coefficients
}
}
\seealso{
\code{\link{multibitTree.load}}, \code{\link{multibitTree.searchFile}}, \code{\link{multibitTree.statistics}}, \code{\link{multibitTree.unload}}
}
\examples{
## get name of example file with fingerprints in package directory

fileB <- file.path(path.package("multibitTree"), "extdata/B.csv")

## load fingerprints from file into memory

multibitTree.load(fileB)

## get name of second example file and open connection

fileA <- file.path(path.package("multibitTree"), "extdata/A.csv")

## search the first 10 prints from file A in loaded file B and print result

con <- file(fileA)
open(con)

for(i in 1:10) print(multibitTree.search(readLines(con, n=1), 0.8))

close(con)

## release memory

multibitTree.unload()
}
\keyword{misc}
