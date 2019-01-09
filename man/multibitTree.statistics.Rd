\name{multibitTree.statistics}
\alias{multibitTree.statistics}
\title{
Print statistics of last search or searchFile operation
}
\description{
This function prints the statistic results of the last search or searchFile operation. 
}
\usage{
multibitTree.statistics()
}
\value{
The function returns a data.frame with three columns:
\item{Checkpoint}{
  this column contains the checkpoint name
}
\item{Count}{
  this column contains the corresponding meassured value
}
\item{Percentage}{
  this column contains the value Count in relation to the total number of searches
}
}
\seealso{
\code{\link{multibitTree.load}}, \code{\link{multibitTree.search}}, \code{\link{multibitTree.searchFile}}, \code{\link{multibitTree.unload}}
}
\examples{
## get name of example file with fingerprints in package directory

fileB <- file.path(path.package("multibitTree"), "extdata/B.csv")

## load fingerprints from file into memory

multibitTree.load(fileB)

## get name of second example file and open connection

fileA <- file.path(path.package("multibitTree"), "extdata/A.csv")

## search all prints from file A in loaded file B and print result

print(multibitTree.searchFile(fileA, 0.95))

## print statistics

multibitTree.statistics()

## release memory

multibitTree.unload()
}
\keyword{misc}
