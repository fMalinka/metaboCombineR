sortDataFrameByRt <- function(data, mzprecision=2)
{
  mzprec <- 10^mzprecision
  mzs <- signif(as.numeric(regmatches(rownames(data), regexpr("\\d+(\\.\\d+)?",rownames(data), perl=TRUE)))*mzprec)/mzprec
  rts <- as.numeric(substring(regmatches(rownames(data), regexpr("T\\d+(\\.\\d+)?",rownames(data), perl=TRUE)),2))
  
  rtsOrder <- order(rts)
	newdata <- data[rtsOrder,]
  return(as.matrix(newdata))
}

getRTs <- function(data)
{
	return(as.numeric(substring(regmatches(rownames(data), regexpr("T\\d+(\\.\\d+)?",rownames(data), perl=TRUE)),2)))
}

getMZs <- function(data, mzprecision=2)
{
    mzprec <- 10^mzprecision
	return(as.numeric(signif(as.numeric(regmatches(rownames(data), regexpr("\\d+(\\.\\d+)?",rownames(data), perl=TRUE)))*mzprec)/mzprec))
}
	
sortMzByRt <- function(data, mzprecision=2)
{
  mzprec <- 10^mzprecision
  mzs <- signif(as.numeric(regmatches(rownames(data), regexpr("\\d+(\\.\\d+)?",rownames(data), perl=TRUE)))*mzprec)/mzprec
  rts <- as.numeric(substring(regmatches(rownames(data), regexpr("T\\d+(\\.\\d+)?",rownames(data), perl=TRUE)),2))
  
  rtsOrder <- order(rts)
  newdata <- data[rtsOrder,]
  newdata <- cbind(newdata, RT= rts[rtsOrder])
  newdata <- cbind(newdata, MZ= mzs[rtsOrder])
  return(as.matrix(newdata))
}

prepareExperiments <- function(listOfExperiment, listOfISindexes)
{
  allExps <- list()
  allIS <- list()
  i <- 1
  for(exp in listOfExperiment)
  {
    allExps <- append(allExps, list(sortMzByRt(exp)))
    allIS <- append(allIS, list(listOfISindexes[[i]]))
    i <- i + 1
  }
  return(list(experiments = allExps, ISindexes = allIS))
}

runMetaboCombiner <- function(listExperimens, mzprecision = 3, windowsize = 5)
{
	myclass <- new(metaboCombineR)
	finalmatrix <- myclass$run(listExperimens, mzprecision, windowsize)
	rm(myclass)
	return(finalmatrix)
}
