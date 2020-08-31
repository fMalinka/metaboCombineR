sortDataFrameByRt <- function(data, mzprecision=2)
{
    mzprec <- 10^mzprecision
    mzs <- trunc(as.numeric(regmatches(rownames(data),
    regexpr("\\d+(\\.\\d+)?",rownames(data), perl=TRUE)))*mzprec)/mzprec
    
    rts <- as.numeric(substring(regmatches(rownames(data),
    regexpr("T\\d+(\\.\\d+)?",rownames(data), perl=TRUE)),2))
    rtsOrder <- order(rts)
    newdata <- data[rtsOrder,]
    return(as.matrix(newdata))
}

sortDataFrameByMz <- function(data, mzprecision=2)
{    
    mzprec <- 10^mzprecision
    mzs <- trunc(as.numeric(regmatches(rownames(data),
    regexpr("\\d+(\\.\\d+)?",rownames(data), perl=TRUE)))*mzprec)/mzprec
    
    mzsOrder <- order(mzs)
    newdata <- data[mzsOrder,]
    return(as.matrix(newdata))
}

getRTs <- function(data)
{
    return(as.numeric(substring(regmatches(rownames(data),
    regexpr("T\\d+(\\.\\d+)?",rownames(data), perl=TRUE)),2)))
}

getMZs <- function(data, mzprecision=2)
{
    mzprec <- 10^mzprecision
    return(as.numeric(trunc(as.numeric(regmatches(rownames(data),
    regexpr("\\d+(\\.\\d+)?",rownames(data), perl=TRUE)))*mzprec)/mzprec))
}

sortMzByRt <- function(data, mzprecision=2)
{
    mzprec <- 10^mzprecision
    mzs <- trunc(as.numeric(regmatches(rownames(data),
    regexpr("\\d+(\\.\\d+)?",rownames(data), perl=TRUE)))*mzprec)/mzprec
    rts <- as.numeric(substring(regmatches(rownames(data),
    regexpr("T\\d+(\\.\\d+)?",rownames(data), perl=TRUE)),2))  
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

#check List
isValid <- function(dataframe)
{
    mzs <- rownames(dataframe)
    rownameFormat <- all(grepl(
        pattern = "^M(\\d+|\\d+\\.\\d+)T(\\d+|\\d+\\.\\d+)", mzs))
    valueFormat <- all(vapply(dataframe, is.numeric, logical(1)))
    rowDim <- nrow(dataframe) > 0
    colDim <- ncol(dataframe) > 1
    return(rownameFormat & valueFormat & rowDim & colDim)
}


runMetaboCombiner <- function(listExperimens, mzprecision = 3, windowsize = 5, algorithm="kmer")
{
    if(algorithm == "rtcor")
    {
        write(paste0("Algorithm: rtcorrectedAlignment"), stdout())
    }
    else if(algorithm == "kmer")
    {
        write(paste0("Algorithm: kmersAlignment"), stdout())
    }
    else
    {
        message(paste0("algorithm parameter '", algorithm, "' is unknown. Default is kmer."))
        algorithm="kmer"
        message("Algorithm: kmersAlignment")
    }
    myclass <- new(metaboCombineR)
    ilist  <- 1
    listValid <- rep(FALSE, length(listExperimens))
    for(elist in listExperimens)
    {
        if(isValid(elist))
        {
            write(paste0("experiment n.", ilist, " has ", ncol(elist),
                " samples and ", nrow(elist),
                " features. Format seems VALID."), stderr())
            listValid[ilist] <- TRUE
        }
        else
        {
            write(paste0("experiment n.", ilist, ". Format seems INVALID."),
                stderr())
        }
        ilist <- ilist + 1
    }
    if(all(listValid))
    {
        if(algorithm == "kmer")
        {
            finalmatrix <- myclass$run(listExperimens, mzprecision, windowsize)
        }
        else
        {
            rtwindowsize <- 50;
            finalmatrix <- myclass$runRT(listExperimens, mzprecision, rtwindowsize)
        }
    }
    else
    {
        write(paste0("ERROR: Check INPUT files!"), stderr())
        finalmatrix  <- NULL
    }
    rm(myclass)
    return(finalmatrix)
}

correctClassiqAling2 <- function(mzvect1, mzvect2, orderLimit = 50)
{
  swapIndex1 <- vector()
  swapIndex2 <- vector()
  i <- 1
  for(vect1 in mzvect1)
  {
    refind <- which(c(mzvect1[i], mzvect2[i]) == "-1")
    if(length(refind) == 1)
    {
      if(refind == 1)
      {
        if(mzvect2[i] %in% mzvect1[(i+1):min((i+orderLimit), length(mzvect1))])
        {
          thesameInd <- which(mzvect2[i] == mzvect1[(i+1):min((i+orderLimit), length(mzvect1))])
          mzvect1[thesameInd + i] <- "-1"
          swapIndex2 <- append(swapIndex2, i)
          swapIndex1 <- append(swapIndex1, thesameInd + i)
        }
      }
      if(refind == 2)
      {
        if(mzvect1[i] %in% mzvect2[(i+1):min((i+orderLimit), length(mzvect2))])
        {
          thesameInd <- which(mzvect1[i] == mzvect2[(i+1):min((i+orderLimit), length(mzvect2))])
          mzvect2[thesameInd + i] <- "-1"
          swapIndex2 <- append(swapIndex2, thesameInd + i)
          swapIndex1 <- append(swapIndex1, i)
        }
      }
    }
    i <- i + 1
  }
  return(list(swapIndex1, swapIndex2))
}
