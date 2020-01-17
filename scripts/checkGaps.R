mydata <- read.csv("../experiments/mymatrix_4.csv")
print(paste0("number of peaks: ", nrow(mydata)))
ngaps <- hist(mydata$GAP, freq = TRUE)
print(paste0("number of peaks: ", nrow(mydata)))

for(i in 1:max(mydata$GAP))
{
  print(paste0(i, " gaps: #", sum(mydata$GAP == i)))
}


nnullgapsIndex <- mydata$GAP == 0
x <- mydata[nnullgapsIndex, 1:(ncol(mydata)-3)]
plot(1:ncol(x),x[1,], type = "l", xlab = c(1, ncol(x)))


exp1 <- 1
exp2 <- 7
exp3 <- 13
exp4 <- 19

print(paste0("gaps in exp1: ", sum(is.na(mydata[,exp1])), " of ", nrow(mydata)))
print(paste0("gaps in exp2: ", sum(is.na(mydata[,exp2])), " of ", nrow(mydata)))
print(paste0("gaps in exp3: ", sum(is.na(mydata[,exp3])), " of ", nrow(mydata)))
print(paste0("gaps in exp4: ", sum(is.na(mydata[,exp4])), " of ", nrow(mydata)))

rt1 <- read.csv("rt1.csv")
rt1 <- rt1[,-1]

rt2 <- read.csv("rt2.csv")
rt2 <- rt2[,-1]
mysplineModel <- smooth.spline(x = rt2, y = rt1, spar = 0.0)
predict(mysplineModel, c(97.27))
predict(mysplineModel, c(rt2[1:5]))

plot(loess(rt1 ~ rt2, span=0.2), type="l")
plot(predict(loess(rt1 ~ rt2, span=0.2), ), type="l")

m <-loess(rt1 ~ rt2, span=0.1)

plot(m, type= "l")

predict(m, rt2[1:5])
rt2[1:5]

 mm  <- lm(rt1 ~ rt2)
mm

predict(mm, 97.27)
97.27*0.9995 + 0.7535

plot(loess(mrt1 ~ mrt2, span=0.2), type="l")
mm <- loess(mrt1 ~ mrt2, span=0.2)
predict(mm, 97.27)

 mm <- smooth.spline(mrt1 ~ mrt2)
predict(mm, 97.27)
predict(mm, mrt1)

f <- approx(mrt2, mrt1)
plot(f)

f <- approxfun(mrt2, mrt1)
plot(f)


f <- spline(mrt2, mrt1)
plot(spline(mrt2, mrt1))

f <- splinefun(mrt2, mrt1)
f(97.27)

f<- loess(mrt1 ~ mrt2)
predict(f,97.27)


f <- approxfun(rt2, rt1)
f(mrt2)
mrt1


g <- splinefun(rt2, rt1, method = "fmm", ties = mean)
g(mrt2)
mrt1

library(splines)

x <- 0:100
a <- ns(x, df=5)
plot(a)


require(graphics); require(stats)
ispl <- interpSpline( women$height, women$weight )

a <- interpSpline(mrt2[1:4], mrt1[1:4])
plot(1:max(mrt2[1:4]), predict(a, 1:max(mrt2[1:4]))$y, type = "l")


g <- splinefun(mrt2[1:4], mrt1[1:4], method = "fmm")
plot(1:max(mrt2[1:4]), g(1:max(mrt2[1:4])), type = "l")

g <- approxfun(mrt2[1:4], mrt1[1:4])
plot(1:max(mrt2[1:4]), g(1:max(mrt2[1:4])), type = "l")

g(mrt2[1:4])
mrt1


a <- smooth.spline(mrt2[1:4], mrt1[1:4])


length(rt1)
length(unique(rt1))

rt1[1:10]
rt2[1:10]


nrt2 <- rt2[!duplicated(rt2)]
nrt1 <- rt1[!duplicated(rt2)]

mymodel <- interpSpline(nrt2, nrt1, bSpline = FALSE)
plot(1:max(nrt2), predict(mymodel, 1:max(nrt2))$y, type = "l")
points(nrt2, nrt1)

g <- splinefun(nrt2, nrt1, method = "fmm")
plot(1:max(nrt2), g(1:max(nrt2)), type = "l")

g <- approxfun(nrt2, nrt1, method = "linear", rule = 2:1)
plot(1:max(nrt2), g(1:max(nrt2)), type = "l")
points(nrt2, nrt1)

#final model

nrt2 <- rt2[!duplicated(rt2)]
nrt1 <- rt1[!duplicated(rt2)]

mymodel <- smooth.spline(nrt1 ~ nrt2, all.knots = TRUE, df=2)
plot(1:max(nrt2), predict(mymodel, 1:max(nrt2))$y, type = "l")
points(nrt2, nrt1)

mymodel <- smooth.spline(nrt1 ~ nrt2, all.knots = TRUE, df=2)
plot(1:max(nrt2[1:4]), predict(mymodel, 1:max(nrt2[1:4]))$y, type = "l")
points(nrt2, nrt1)

g <- approxfun(nrt2, nrt1, method = "linear", rule = 2)
plot(1:max(nrt2), g(1:max(nrt2)), type = "l")
points(nrt2, nrt1)

