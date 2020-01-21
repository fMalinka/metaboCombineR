\name{metaboExp2}
\alias{metaboExp2}
\docType{data}
\title{
MetaboCombineR experiment example no. 1
}
\description{
This is the second example of LC/MS dataset in metaboAnalyst format.
}
\usage{data("metaboExp2")}
\format{
A data frame with 2500 observations on the following 6 variables.
\describe{
\item{\code{X064.EPL03_m_Mzb1_ESI.mzML}}{a numeric vector    sample no. 7}
\item{\code{X064.EPL05_m_Mzb1_ESI.mzML}}{a numeric vector    sample no. 8}
\item{\code{X064.EPL07_f_Mzb1_ESI.mzML}}{a numeric vector    sample no. 9}
\item{\code{X64Z.61018_m_WT_ESI.mzML}}{a numeric vector    sample no. 10}
\item{\code{X64Z.61019_m_WT_ESI.mzML}}{a numeric vector    sample no. 11}
\item{\code{X64Z.61020_m_WT_ESI.mzML}}{a numeric vector    sample no. 12}
}
}
\details{
Columns represent samples and rows represent features.
Rownames are in the following format: M[mz]R[rt] where
[mz] is mz value and [rt] is rt value of the given feature,
both are real numbers.
}
\source{
\url{https://github.com/fMalinka/metaboCombiner}
}
\references{
Handling elution order and retention time swaps in large-scale multi-batch
untargeted LC/MS studies
}
\examples{
data(metaboExp2)
}
\keyword{datasets}
