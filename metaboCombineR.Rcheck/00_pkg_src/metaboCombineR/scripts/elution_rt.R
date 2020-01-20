
sortMzByRt <- function(data)
{
  mzs <- trunc(as.numeric(regmatches(rownames(data), regexpr("\\d+(\\.\\d+)?",rownames(data), perl=TRUE)))*100)/100
  rts <- as.numeric(substring(regmatches(rownames(data), regexpr("T\\d+(\\.\\d+)?",rownames(data), perl=TRUE)),2))
  
  rtsOrder <- order(rts)
  newdata <- data[rtsOrder,]
  newdata <- cbind(newdata, RT= rts[rtsOrder])
  newdata <- cbind(newdata, MZ= mzs[rtsOrder])
  return(as.matrix(newdata))
}




mzb1 <- read.csv("/home/frantisek/BIOCEV/metaboDatasets/mzb1_parts/mzb1_wnorm_1.csv", row.names = 1, stringsAsFactors = FALSE)
mzb2  <- read.csv("/home/frantisek/BIOCEV/metaboDatasets/v2/mzb1v2_b.csv", row.names = 1)
mzb3 <- read.csv("/home/frantisek/BIOCEV/metaboDatasets/mzb1_parts/mzb1_norm_3.csv", row.names = 1)
#mzb4 <- read.csv("/home/frantisek/BIOCEV/metaboDatasets/mzb1_parts/mzb1_norm_4.csv", row.names = 1)

mzb1 <- mzb1[-1,]
mzb1_1 <- apply(mzb1, 2, as.numeric)
rownames(mzb1_1) <- rownames(mzb1)
mzb1 <- mzb1_1

mzb1sort <- sortMzByRt(mzb1)
rt1 <- mzb1sort[, "RT"]
subrt10 <- rt1[rt1 < 150 ]
subrt10
hist(subrt10, breaks = seq(96, 150))



a <- hist(rt1, breaks = seq(trunc(min(rt1)), trunc(max(rt1))+1, length.out = 4000))
sum(a$counts == 1)


oneIndex <- a$breaks[a$counts == 1]
intensityMean <- rep(0, length(oneIndex))
ii <- 1
for(i in oneIndex)
{
  intensityMean[ii] <- mean(mzb1sort[which(trunc(rt1) == i), -c(ncol(mzb1sort) -1, ncol(mzb1sort))], na.rm = TRUE)
  ii <- ii + 1
}

oneIndex2 <- a$breaks[a$counts > 10]
intensityMean2 <- rep(0, length(oneIndex2))
ii2 <- 1
for(i in oneIndex2)
{
  intensityMean2[ii2] <- mean(mzb1sort[which(trunc(rt1) == i), -c(ncol(mzb1sort) -1, ncol(mzb1sort))], na.rm = TRUE)
  ii2 <- ii2 + 1
}

summary(intensityMean)
summary(intensityMean2)

#trend prediction
trend <- lm(a$counts ~ a$breaks[1:(length(a$breaks)-1)])

var(rt1)
sd(rt1)
