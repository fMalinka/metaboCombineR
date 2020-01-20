#check elution order
mydata <- read.csv("../experiments/mymatrix_2_monointerpolator.csv")

rt <- mydata$RT
mz <- mydata$MZ
gap <- mydata$GAP


elutionOrder <- vector()
elutionNum <- vector()
elutioNumZero <- vector()
elutionIndex <- vector()
elutionRT <- vector()

for(i in 1:length(rt))
{
  if(!is.na(rt[i]))
  {
    rt[rt[i:length(rt)] < rt[i]+2]
    id <- i + which(rt[i:length(rt)] < rt[i]+2)
    elutionID <- i + which(mz[id] == mz[i])
    if(length(elutionID) > 0)
    {
      elutionOrder <- append(elutionOrder, mz[i])
      elutionRT <- append(elutionRT, rt[i])
      elutionNum <- append(elutionNum, length(elutionID))
      elutionIndex <- append(elutionIndex, i)
      elutioNumZero <- append(elutioNumZero, sum(gap[i:elutionID] == 0))
      rt[elutionID] <- NA
      mz[elutionID] <- NA
    }
    
  }
}

hist(elutionIndex)
hh <- hist(elutionIndex, breaks = seq(4,3626, length.out = 100))
hh

hh <- hist(elutionRT, breaks = seq(trunc(min(elutionRT)), trunc(max(elutionRT))+1, length.out = 4000))
hh

rt[rt[1:length(rt)] < rt[1]+10]
id <- which(rt[1:length(rt)] < rt[1]+10)

sum(mz[id] == 104.1)

which(mz[id] == 104.1)

