\name{multibitTree.unload}
\alias{multibitTree.unload}
\title{
Discard MultibitTree from Memory
}
\description{
This function discards the MultibitTree data structure an frees the
used memory.
}
\usage{
multibitTree.unload()
}
\value{
This function always returns NULL.
}
\seealso{
\code{\link{multibitTree.load}}, \code{\link{multibitTree.search}}
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
