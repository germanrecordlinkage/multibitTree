multibitTree.load <-
function(filename, threads = 1, size = 0, leafLimit = 8) {
	result <- .Call(mbtLoadCall, filename, threads, size, leafLimit)
	return(result)
}
