\name{getMZs}
\alias{getMZs}
\title{Extract MZ values from experiment in mataboAnalyst format.}
\description{
This function extract all MZ values from given dataframe that represents
an experiment.
}
\usage{
getMZs(data, mzprecision = 2)
}
\arguments{
\item{data}{
a data.frame where columns represents samples and rows represents features.
}
\item{mzprecision}{
an integer truncates mz value on given number of digits
}
}
\details{
The experiment data.frame must be in metaboAnalyst format, i.e. columns
represent samples and rows represent features. Description of features
(MZ[mz]T[RT]) is gived as 'rownames'.
}
\value{
a numeric vector of mz values
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
getMZs(metaboExp1)
}
\keyword{mz}
