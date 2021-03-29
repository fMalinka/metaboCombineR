test_inputs_kmers <- function()
{
    data(metaboExp1)
    obs <- tryCatch(runMetaboCombiner(algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("argument \"listExperimens\" is missing, with no default", obs)
    obs <- tryCatch(runMetaboCombiner(list(), algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("List \"listExperimens\" does not contain any experiments!", obs)
    checkException(runMetaboCombiner(list(metaboExp1, obj), algorithm = "kmer", trunc = 2), "Input object not found!")
    obs <- tryCatch(runMetaboCombiner(metaboExp1, algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("Argument \"listExperimens\" is not a List!", obs)
    obs <- tryCatch(runMetaboCombiner(metaboExp1,metaboExp1, algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("Argument \"listExperimens\" is not a List!", obs)
    obs <- tryCatch(runMetaboCombiner(metaboExp1,metaboExp5, algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("Argument \"listExperimens\" is not a List!", obs)
}

test_malformedInputs_kmers <- function()
{
    data(metaboExp1)
    data(metaboExp2)
    dfempty <- data.frame()
    
    test <- metaboExp1
    rownames(test) <- NULL
    obs <- tryCatch(runMetaboCombiner(list(metaboExp2,test), algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("experiment n.2. Format seems INVALID.", obs)
    
    test <- metaboExp1
    rownames(test) <- seq(1, nrow(test))
    obs <- tryCatch(runMetaboCombiner(list(metaboExp2,test), algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("experiment n.2. Format seems INVALID.", obs)
    
    test <- metaboExp1
    rownames(test) <- paste0("M", seq(1, nrow(test)))
    obs <- tryCatch(runMetaboCombiner(list(metaboExp2,test), algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("experiment n.2. Format seems INVALID.", obs)
    
    test <- metaboExp1
    rownames(test) <- paste0("M", seq(1, nrow(test)), "T")
    obs <- tryCatch(runMetaboCombiner(list(metaboExp2,test), algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("experiment n.2. Format seems INVALID.", obs)
    
    test <- metaboExp1
    rownames(test) <- paste0("M", "T", seq(1, nrow(test)))
    obs <- tryCatch(runMetaboCombiner(list(metaboExp2,test), algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("experiment n.2. Format seems INVALID.", obs)
    
    
    test <- metaboExp1
    test <- rbind(group = rep("wt", ncol(test)), test)
    obs <- tryCatch(runMetaboCombiner(list(test,metaboExp1, metaboExp2), algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("There are inconsistencies in GROUP labels. The GROUP label is missing in experiments no. 2,3.", obs)
    obs <- tryCatch(runMetaboCombiner(list(metaboExp1, test, metaboExp2), algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("There are inconsistencies in GROUP labels. The GROUP label is missing in experiments no. 1,3.", obs)
    obs <- tryCatch(runMetaboCombiner(list(metaboExp1, metaboExp2, test), algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("There are inconsistencies in GROUP labels. The GROUP label is missing in experiments no. 1,2.", obs)
    
    obs <- tryCatch(runMetaboCombiner(list(dfempty), algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("List n. 1 has 0 rows!", obs)
    obs <- tryCatch(runMetaboCombiner(list(metaboExp1, dfempty, metaboExp2), algorithm = "kmer", trunc = 2), error=conditionMessage)
    checkIdentical("List n. 2 has 0 rows!", obs)
}

test_inputs_rtcor <- function()
{
    data(metaboExp1)
    obs <- tryCatch(runMetaboCombiner(algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("argument \"listExperimens\" is missing, with no default", obs)
    obs <- tryCatch(runMetaboCombiner(list(), algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("List \"listExperimens\" does not contain any experiments!", obs)
    checkException(runMetaboCombiner(list(metaboExp1, obj), algorithm = "rtcor", trunc = 2), "Input object not found!")
    obs <- tryCatch(runMetaboCombiner(metaboExp1, algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("Argument \"listExperimens\" is not a List!", obs)
    obs <- tryCatch(runMetaboCombiner(metaboExp1,metaboExp1, algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("Argument \"listExperimens\" is not a List!", obs)
    obs <- tryCatch(runMetaboCombiner(metaboExp1,metaboExp5, algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("Argument \"listExperimens\" is not a List!", obs)
}

test_malformedInputs_rtcor <- function()
{
    data(metaboExp1)
    data(metaboExp2)
    dfempty <- data.frame()
    
    test <- metaboExp1
    rownames(test) <- NULL
    obs <- tryCatch(runMetaboCombiner(list(metaboExp2,test), algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("experiment n.2. Format seems INVALID.", obs)
    
    test <- metaboExp1
    rownames(test) <- seq(1, nrow(test))
    obs <- tryCatch(runMetaboCombiner(list(metaboExp2,test), algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("experiment n.2. Format seems INVALID.", obs)
    
    test <- metaboExp1
    rownames(test) <- paste0("M", seq(1, nrow(test)))
    obs <- tryCatch(runMetaboCombiner(list(metaboExp2,test), algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("experiment n.2. Format seems INVALID.", obs)
    
    test <- metaboExp1
    rownames(test) <- paste0("M", seq(1, nrow(test)), "T")
    obs <- tryCatch(runMetaboCombiner(list(metaboExp2,test), algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("experiment n.2. Format seems INVALID.", obs)
    
    test <- metaboExp1
    rownames(test) <- paste0("M", "T", seq(1, nrow(test)))
    obs <- tryCatch(runMetaboCombiner(list(metaboExp2,test), algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("experiment n.2. Format seems INVALID.", obs)
    
    
    test <- metaboExp1
    test <- rbind(group = rep("wt", ncol(test)), test)
    obs <- tryCatch(runMetaboCombiner(list(test,metaboExp1, metaboExp2), algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("There are inconsistencies in GROUP labels. The GROUP label is missing in experiments no. 2,3.", obs)
    obs <- tryCatch(runMetaboCombiner(list(metaboExp1, test, metaboExp2), algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("There are inconsistencies in GROUP labels. The GROUP label is missing in experiments no. 1,3.", obs)
    obs <- tryCatch(runMetaboCombiner(list(metaboExp1, metaboExp2, test), algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("There are inconsistencies in GROUP labels. The GROUP label is missing in experiments no. 1,2.", obs)
    
    obs <- tryCatch(runMetaboCombiner(list(dfempty), algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("List n. 1 has 0 rows!", obs)
    obs <- tryCatch(runMetaboCombiner(list(metaboExp1, dfempty, metaboExp2), algorithm = "rtcor", trunc = 2), error=conditionMessage)
    checkIdentical("List n. 2 has 0 rows!", obs)
}