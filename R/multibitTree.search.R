multibitTree.search <-
function(query, minTanimoto, size = 0, sort = FALSE) {
	result <- .Call(mbtSearchCall, query, minTanimoto, size, sort)
	return(data.frame(result))
}
