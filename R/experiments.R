sortDataFrameByRt <- function(data, mzprecision=2)
{
    mzprec <- 10^mzprecision
    mzs <- signif(as.numeric(regmatches(rownames(data),
    regexpr("\\d+(\\.\\d+)?",rownames(data), perl=TRUE)))*mzprec)/mzprec
    
    rts <- as.numeric(substring(regmatches(rownames(data),
    regexpr("T\\d+(\\.\\d+)?",rownames(data), perl=TRUE)),2))
    rtsOrder <- order(rts)
    newdata <- data[rtsOrder,]
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
    return(as.numeric(signif(as.numeric(regmatches(rownames(data),
    regexpr("\\d+(\\.\\d+)?",rownames(data), perl=TRUE)))*mzprec)/mzprec))
}

sortMzByRt <- function(data, mzprecision=2)
{
    mzprec <- 10^mzprecision
    mzs <- signif(as.numeric(regmatches(rownames(data),
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
    valueFormat <- all(sapply(dataframe, is.numeric))
    rowDim <- nrow(dataframe) > 0
    colDim <- ncol(dataframe) > 0
    return(rownameFormat & valueFormat & rowDim & colDim)
}


runMetaboCombiner <- function(listExperimens, mzprecision = 3, windowsize = 5)
{
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
        finalmatrix <- myclass$run(listExperimens, mzprecision, windowsize)
    }
    else
    {
        write(paste0("ERROR: Check INPUT files!"), stderr())
        finalmatrix  <- NULL
    }
    rm(myclass)
    return(finalmatrix)
}
