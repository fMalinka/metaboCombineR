#' Sort experiment according to RT
#' @usage sortDataFrameByRt(experiment)
#' @param experiment an input data.frame representing experiments
#' @return Vector of MZs
#' @examples
#' data(metaboExp1)
#' sortDataFrameByRt(metaboExp1)
sortDataFrameByRt <- function(data)
{
    mzs <- as.numeric(regmatches(rownames(data),
           regexpr("\\d+(\\.\\d+)?",rownames(data), perl=TRUE)))
    
    rts <- as.numeric(substring(regmatches(rownames(data),
    regexpr("T\\d+(\\.\\d+)?",rownames(data), perl=TRUE)),2))
    rtsOrder <- order(rts)
    newdata <- data[rtsOrder,]
    return(as.matrix(newdata))
}

#' Sort experiment according to MZ
#' @usage sortDataFrameByMz(experiment)
#' @param experiment an input data.frame representing experiments
#' @return Data.frame sorted according to MZs
#' @examples
#' data(metaboExp1)
#' sortDataFrameByMz(metaboExp1)
sortDataFrameByMz <- function(data)
{    
    mzs <- as.numeric(regmatches(rownames(data),
    regexpr("\\d+(\\.\\d+)?",rownames(data), perl=TRUE)))
    
    mzsOrder <- order(mzs)
    newdata <- data[mzsOrder,]
    return(as.matrix(newdata))
}

#' Get vector of RT
#' @usage getRTs(experiment)
#' @param experiment an input data.frame representing experiments
#' @return Vector of RTs
#' @examples
#' data(metaboExp1)
#' getRTs(metaboExp1)
getRTs <- function(data)
{
    return(as.numeric(substring(regmatches(rownames(data),
    regexpr("T\\d+(\\.\\d+)?",rownames(data), perl=TRUE)),2)))
}

#' Get vector of MZ
#' @usage getMZs(experiment)
#' @param experiment an input data.frame representing experiments
#' @return Vector of MZs
#' @examples
#' data(metaboExp1)
#' getMZs(metaboExp1)
getMZs <- function(data)
{
    return(as.numeric(regmatches(rownames(data),
          regexpr("\\d+(\\.\\d+)?",rownames(data), perl=TRUE))))
}

#check for valid format of input
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

#' Run metaboCombineR
#' @usage runMetaboCombiner(listExperimens, windowsize, algorithm, trunc)
#' @param listExperimens a list of data.frames representing experiments
#' @param windowsize a number indicating window size parameter for 'rtcor', or k-mer size for 'k-mer' algorithm.
#' @param algorithm a string, 'rtcor' or 'kmer' algorithm
#' @param ppm a relative value of m/z difference for the peak/feature matching
#' @param abs an absolute value of m/z difference for the peak/feature matching
#' @param trunc a number of decimal places to be truncated for m/z peak/feature matching
#' @param s_match a score for matched features, in an alignment problem (default: 100)
#' @param s_del a penalization for deletion of features, in an alignment problem (default: 1)
#' @param s_ins a penalization for insertion of features, in an alignment problem (default: 1)
#' @param mzres reported mz values reported in pairwise alignment according to reference (ref), non-reference (nonref) experiments, or its average (avg). (default: avg)
#' @return Data.frame of final matrix
#' @details Only one parameter of 'ppm', 'abs', or 'trunc' can be specified.
#' @examples
#' data(metaboExp1)
#' data(metaboExp2)
#' data(metaboExp3)
#' data(metaboExp4)
#' runMetaboCombiner(list(metaboExp1, metaboExp2, metaboExp3, metaboExp4), algorithm = "kmer", trunc = 2)
runMetaboCombiner <- function(listExperimens, windowsize = 5, algorithm="kmer", ppm = -1, abs = -1, trunc = -1, s_match = 100, s_del = 1, s_ins = 1, mzres = "avg")
{
    matchMode <- NULL
    matchValue <- 0
    if(class(listExperimens) != "list")
    {
        stop("Argument \"listExperimens\" is not a List!")
    }
    if(length(listExperimens) == 0)
    {
        stop("List \"listExperimens\" does not contain any experiments!")
    }
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
    if(!is.numeric(ppm) | !is.numeric(abs) | !is.numeric(trunc))
    {
        stop("ppm, abs, and trunc parameters must be numeric!")
    }
    
    if(sum(c(ppm >= 0,  abs >= 0 , trunc >= 0)) >= 2)
    {
        stop("Only one parameter of ppm, abs, or trunc can be set up!")
    } 
    
    if(mzres != "ref" & mzres != "nonref" & mzres != "avg")
    {
        stop("Only values ref, nonref, or avg are supported for mzres!")
    }
    
    if(ppm != -1)
    {
        if(ppm >= 0)
        {
            if(is.null(matchMode))
            {
                matchMode <- "ppm"
                matchValue <- ppm
                write(paste0("ppm: ", ppm), stdout())
                } else {
                    stop("Only one parameter of ppm, abs, or trunc can be set up!")
                }
        } else
        {
            stop("ppm parameter must be > 0!")
        }
    } else if(abs != -1)
    {
        if(abs >= 0)
        {
            if(is.null(matchMode))
            {
                matchMode <- "abs"
                matchValue <- abs
                write(paste0("abs: ", abs), stdout())
                } else {
                    stop("Only one parameter of ppm, abs, or trunc can be set up!")
                }
        } else
        {
            stop("abs parameter must be > 0!")
        }
    } else if(trunc != -1)
    {
        if(trunc >= 0)
        {
            if(is.null(matchMode))
            {
                matchMode <- "trunc"
                matchValue <- trunc
                write(paste0("trunc: ", trunc), stdout())
                } else {
                    stop("Only one parameter of ppm, abs, or trunc can be set up!")
                }
        } else
        {
            stop("trunc parameter must be > 0!")
        }
    } else
    {
        stop("At least one of ppm, abs, or trunc parameters must be set up correctly!")
    }
    myclass <- new(metaboCombineR)
    ilist  <- 1
    isGroup <- rep(FALSE, length(listExperimens))
    listValid <- rep(FALSE, length(listExperimens))
    listGroup <- vector(mode = "list", length = length(listExperimens))
    for(elist in listExperimens)
    {
        if(nrow(elist) == 0)
        {
            stop(paste0("List n. ", ilist, " has 0 rows!"))
        }
      #group row present
        if(grepl("group", rownames(elist)[1], ignore.case=TRUE))
        {
          listGroup[ilist] <- list(elist[1,])
          elist <- elist[-1,]
          isGroup[ilist] <- TRUE
        }
        orn <- rownames(elist)
        elist <- apply(apply(elist,2, as.character), 2, as.numeric)
        rownames(elist) <- orn
        elist <- as.data.frame(elist)
        listExperimens[[ilist]] <- elist
        if(isValid(elist))
        {
            write(paste0("experiment n.", ilist, " has ", ncol(elist),
                " samples and ", nrow(elist),
                " features. Format seems VALID."), stderr())
            listValid[ilist] <- TRUE
        }
        else
        {
            stop(paste0("experiment n.", ilist, ". Format seems INVALID."))
        }
        ilist <- ilist + 1
    }
    if(!(all(isGroup)|(all(!isGroup))))
    {
        stop(paste0("There are inconsistencies in GROUP labels. The GROUP label is missing in experiments no. ", paste0(which(isGroup == FALSE), collapse=","), "."))
    }
    else if(all(listValid))
    {
        if(algorithm == "kmer")
        {
			if(windowsize <=3 | windowsize >= min(sapply(listExperimens, nrow)))
			{
			    stop("ERROR: windowsize parameter is out of range!")
			} else
			{
			    finalmatrix <- myclass$runKmersAlignment(listExperimens, matchValue, windowsize, matchMode, s_match, s_del, s_ins, mzres)
			}
        }
        else
        {
			if(windowsize < 0 | windowsize >= min(sapply(listExperimens, nrow)))
			{
			    stop("ERROR: windowsize parameter is out of range!")
			} else
			{
			    finalmatrix <- myclass$runRtcorrectedAlignment(listExperimens, matchValue, windowsize, matchMode, s_match, s_del, s_ins, mzres)
			}
        }
        
        if(all(isGroup))
        {
			gorder <- myclass$getAlignmentOrder()
			finalmatrix <- rbind(group = unlist(unlist(listGroup[gorder+1])), finalmatrix)
		}
    }
    else
    {
        stop("ERROR: Check INPUT files!")
    }
    rm(myclass)
    return(finalmatrix)
}