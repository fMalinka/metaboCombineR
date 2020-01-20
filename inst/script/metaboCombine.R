#load data
#mzb1 <- read.csv("/home/frantisek/BIOCEV/metaboDatasets/mzb1_parts/mzb1_norm_1.csv", row.names = 1)
#mzb2  <- read.csv("/home/frantisek/BIOCEV/metaboDatasets/mzb1_parts/mzb1_norm_2.csv", row.names = 1)
#mzb3 <- read.csv("/home/frantisek/BIOCEV/metaboDatasets/mzb1_parts/mzb1_norm_3.csv", row.names = 1)
#mzb4 <- read.csv("/home/frantisek/BIOCEV/metaboDatasets/mzb1_parts/mzb1_norm_4.csv", row.names = 1)

#/media/frantisek/data/BIOCEV/metaboDatasets/mzb1_parts
mzb1 <- read.csv("data/mzb1_norm_1.csv", row.names = 1)
mzb2  <- read.csv("data/mzb1_norm_2.csv", row.names = 1)
mzb3 <- read.csv("data/mzb1_norm_3.csv", row.names = 1)
mzb4 <- read.csv("data/mzb1_norm_4.csv", row.names = 1)


getMZ <- function(datatable)
{
	return(sub("T.*", "", rownames(datatable)))
}


library(metaboCombineR)
mz <- runMetaboCombiner(list(mzb3, mzb4), mzprecision = 2)
vectormz <- getMZ(mz)





###### testing part ######


test <- function(orig, result)
{
	mzo <- getMZs(orig)
	toret <- vector()
	for(i in 1:length(mzo))
	{
		if(!any(grep(mzo[i], rownames(result))))
		{
#			print(paste0(mzo[i], " neni"))
			toret <- append(toret, mzo[i])
		}
	}
return(toret)
}





#mztable <- table(as.numeric(sub("T.*", "",rownames(mz))))
#mzsampe <- table(union(getMZs(mzb1), getMZs(mzb2)))

#library(utils)
library(metaboCombineR)
mz <- runMetaboCombiner(list(mzb3, mzb4), mzprecision = 2)






mz <- runMetaboCombiner(list(mzb1, mzb1), mzprecision = 2)
write.csv(mz, "mz.csv")
mz <- runMetaboCombiner(list(mzb2, mzb2), mzprecision = 2)
write.csv(mz, "mz.csv")
mz <- runMetaboCombiner(list(mzb3, mzb3), mzprecision = 2)
write.csv(mz, "mz.csv")
mz <- runMetaboCombiner(list(mzb4, mzb4), mzprecision = 2)
write.csv(mz, "mz.csv")

mz <- runMetaboCombiner(list(mzb1, mzb2), mzprecision = 2)
test(mzb1, mz)
test(mzb2, mz)

mz <- runMetaboCombiner(list(mzb1, mzb3), mzprecision = 2)
test(mzb1, mz)
test(mzb3, mz)

mz <- runMetaboCombiner(list(mzb1, mzb4), mzprecision = 2)
test(mzb1, mz)
test(mzb4, mz)

mz <- runMetaboCombiner(list(mzb2, mzb4), mzprecision = 2)
test(mzb2, mz)
test(mzb4, mz)

mz <- runMetaboCombiner(list(mzb3, mzb4), mzprecision = 2)
test(mzb3, mz)
test(mzb4, mz)

mz <- runMetaboCombiner(list(mzb1, mzb2, mzb3, mzb4), mzprecision = 2)
write.csv(mz, "mz.csv")


