\name{getRTs}
\alias{getRTs}
\title{Extract RT values from experiment in mataboAnalyst format.}
\description{
This function extract all RT values from given dataframe that represents
an experiment.
}
\usage{
getRTs(data)
}
\arguments{
\item{data}{
a data.frame where columns represents samples and rows represents features.
}
}
\details{
The experiment data.frame must be in metaboAnalyst format, i.e. columns
represent samples and rows represent features. Description of features
(MZ[mz]T[RT]) is gived as 'rownames'.
}
\value{
a numeric vector of rt values
}
\references{
\url{https://github.com/fMalinka/metaboCombiner}
}
\author{
Frantisek Malinka (frantisek.malinka@gmail.com)
}
\seealso{
Handling elution order and retention time swaps in large-scale multi-batch
untargeted LC/MS studies
}
\examples{
##load experiment
data(metaboExp1)
getRTs(metaboExp1)
}
\keyword{mz}
