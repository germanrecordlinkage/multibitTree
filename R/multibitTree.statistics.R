multibitTree.statistics <-
function() {
  options("scipen"=16)
	result <- .Call(mbtStatisticsCall)
	return(data.frame(result))
}
