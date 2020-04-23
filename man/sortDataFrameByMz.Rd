\name{sortDataFrameByMz}
\alias{sortDataFrameByMz}
\title{Sort experimet data.frame by MZ}
\description{
This function returns a data.frame where rows are sorted
by the corresponding mz values.
}
\usage{
sortDataFrameByMz(data, mzprecision = 2)
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
Sorted data.frame.
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
\keyword{data.frame}
\keyword{experiment}
