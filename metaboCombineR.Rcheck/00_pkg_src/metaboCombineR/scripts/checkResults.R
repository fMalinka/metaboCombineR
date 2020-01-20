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

mzb1 <- read.csv("/home/frantisek/BIOCEV/metaboDatasets/mzb1_parts/mzb1_norm_1.csv", row.names = 1)
mzb2  <- read.csv("/home/frantisek/BIOCEV/metaboDatasets/mzb1_parts/mzb1_norm_2.csv", row.names = 1)
mzb3 <- read.csv("/home/frantisek/BIOCEV/metaboDatasets/mzb1_parts/mzb1_norm_3.csv", row.names = 1)
mzb4 <- read.csv("/home/frantisek/BIOCEV/metaboDatasets/mzb1_parts/mzb1_norm_4.csv", row.names = 1)
mzball  <- read.csv("../experiments/mzb1_all.csv", row.names = 1)
mzbres <- read.csv("../experiments/mymatrix_4.csv")

sort_mzb1 <- sortMzByRt(mzb1)
sort_mzb2 <- sortMzByRt(mzb2)
sort_mzb3 <- sortMzByRt(mzb3)
sort_mzb4 <- sortMzByRt(mzb4)
sort_mzball <- sortMzByRt(mzball)

###TODO
### test udelat obecny
### odchytit prave problemy s poradim viz 319


myres <- read.csv("../experiments/mymatrix_3.csv")
myres <- myres[,-1]
hres <- hist(myres$GAP)
hres

expCol <- ncol(sort_mzb1)-2
lastIn1 <- 0
lastIn2 <- 0
lastIn3 <- 0

testVector <- rep(FALSE, nrow(myres))
for(i in 1:nrow(myres))
{
  origin1 <- sort_mzb1[which(myres$MZ[i] == sort_mzb1[(lastIn1+1):nrow(sort_mzb1),"MZ"])[1], 1:expCol]
  origin2 <- sort_mzb2[which(myres$MZ[i] == sort_mzb2[(lastIn2+1):nrow(sort_mzb2),"MZ"])[1], 1:expCol]
  origin3 <- sort_mzb3[which(myres$MZ[i] == sort_mzb3[(lastIn3+1):nrow(sort_mzb3),"MZ"])[1], 1:expCol]
  
  if(!all(is.na(origin1)) && !all(is.na(myres_mzb1)))
  {
    lastIn1 <- which(myres$MZ[i] == sort_mzb1[(lastIn1+1):nrow(sort_mzb1),"MZ"])[[1]]
  }
  if(!all(is.na(origin2)) && !all(is.na(myres_mzb2)))
  {
    lastIn2 <- which(myres$MZ[i] == sort_mzb2[(lastIn2+1):nrow(sort_mzb2),"MZ"])[[1]]
  }
  if(!all(is.na(origin3)) && !all(is.na(myres_mzb3)))
  {
    lastIn3 <- which(myres$MZ[i] == sort_mzb3[(lastIn3+1):nrow(sort_mzb3),"MZ"])[[1]]
  }
  
  myres_mzb3 <- myres[i,1:6]
  myres_mzb1 <- myres[i,7:12]
  myres_mzb2 <- myres[i,13:18]
  
  l1 <- all(origin1 == myres_mzb1 || is.na(origin1) == is.na(myres_mzb1))
  l2 <- all(origin2 == myres_mzb2 || is.na(origin2) == is.na(myres_mzb2))
  l3 <- all(origin3 == myres_mzb3 || is.na(origin3) == is.na(myres_mzb3))
  if(all(l1,l2,l3, na.rm = TRUE))
  {
    testVector[i] <- TRUE
  }
}
#myres$MZ[1] %in% sort_mzb2[,"MZ"]
origin1 <- sort_mzb1[which(myres$MZ[1] == sort_mzb1[(lastIn1+1):nrow(sort_mzb1),"MZ"])[1], 1:expCol]
origin2 <- sort_mzb2[which(myres$MZ[1] == sort_mzb2[(lastIn2+1):nrow(sort_mzb2),"MZ"])[1], 1:expCol]
origin3 <- sort_mzb3[which(myres$MZ[1] == sort_mzb3[(lastIn3+1):nrow(sort_mzb3),"MZ"])[1], 1:expCol]

if(!all(is.na(origin1)))
{
  lastIn1 <- which(myres$MZ[1] == sort_mzb1[(lastIn1+1):nrow(sort_mzb1),"MZ"])[[1]]
}
if(!all(is.na(origin2)))
{
  lastIn2 <- which(myres$MZ[1] == sort_mzb2[(lastIn2+1):nrow(sort_mzb2),"MZ"])[[1]]
}
if(!all(is.na(origin3)))
{
  lastIn3 <- which(myres$MZ[1] == sort_mzb3[(lastIn3+1):nrow(sort_mzb3),"MZ"])[[1]]
}

myres_mzb3 <- myres[1,1:6]
myres_mzb1 <- myres[1,7:12]
myres_mzb2 <- myres[1,13:18]

l1 <- all(origin1 == myres_mzb1 || is.na(origin1) == is.na(myres_mzb1))
l2 <- all(origin2 == myres_mzb2 || is.na(origin2) == is.na(myres_mzb2))
l3 <- all(origin3 == myres_mzb3 || is.na(origin3) == is.na(myres_mzb3))
all(l1,l2,l3)
