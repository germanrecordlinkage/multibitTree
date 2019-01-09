\name{multibitTree.load}
\alias{multibitTree.load}
\title{
Load Fingerprints from File into MultibitTree
}
\description{
This function reads a set of fingerprints from a file and stores
them into the multibit search tree. If there is already a set
of prints loaded, it will be discarded.
}
\usage{
multibitTree.load(filename, threads = 1, size = 0, leafLimit = 8)
}
\arguments{
  \item{filename}{
  a character string giving the name of the file to load
}
  \item{threads}{
  the number of parallel threads that shall be used to construct the search tree
}
  \item{size}{
  the number of fingerprints that shall be loaded into the search tree (0 = all)
}
  \item{leafLimit}{
  the maximum number of fingerprints for which no further sub-tree shall be calculated
}
}
\value{
returns the number of fingerprints that could actually be loaded
}
\seealso{
\code{\link{multibitTree.search}}, \code{\link{multibitTree.unload}}
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
