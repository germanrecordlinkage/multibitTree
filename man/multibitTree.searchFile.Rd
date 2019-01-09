\name{multibitTree.searchFile}
\alias{multibitTree.searchFile}
\title{
Search multiple Fingerprints from input file in MultibitTree
}
\description{
This function searches in a loaded MultibitTree. With a given input filename and
Tanimoto coefficient the matching fingerprints for all query-prints in the input
file will be returned.
}
\usage{
multibitTree.searchFile(filename, minTanimoto, resultFile = "", seperator = ",")
}
\arguments{
  \item{filename}{
  a character string containing the filename of the input file
}
  \item{minTanimoto}{
  a numeric value giving the lower bound of tanimoto coefficient to search for
}
  \item{resultFile}{
  an optional character string containing the filename of the result file
}
  \item{seperator}{
  an optional character string specifying the column seperator string for the result file
}
}
\value{
The function returns a data.frame with three columns. If a result file is specified the data.frame
will be empty and the results are written as csv-file instead.
\item{query}{
  this column contains the query fingerprint ids
}
\item{fingerprint}{
  this column contains the matching fingerprint ids
}
\item{tanimoto}{
  this column contains the corresponding Tanimoto coefficients
}
}
\seealso{
\code{\link{multibitTree.load}}, \code{\link{multibitTree.search}}, \code{\link{multibitTree.statistics}}, \code{\link{multibitTree.unload}}
}
\examples{
## get name of example file with fingerprints in package directory

fileB <- file.path(path.package("multibitTree"), "extdata/B.csv")

## load fingerprints from file into memory

multibitTree.load(fileB)

## get name of second example file and open connection

fileA <- file.path(path.package("multibitTree"), "extdata/A.csv")

## search all prints from file A in loaded file B and print result

print(multibitTree.searchFile(fileA, 0.8))

## search all prints from file A in loaded file B and store results in file C

multibitTree.searchFile(fileA, 0.8, "C.csv");

## release memory

multibitTree.unload()
}
\keyword{misc}
