multibitTree.searchFile <-
function(filename, minTanimoto, resultFile = "", seperator = ",") {
	result <- .Call(mbtSearchFileCall, filename, minTanimoto, resultFile, seperator)
	return(data.frame(result))
}
