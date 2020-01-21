\name{metaboExp1}
\alias{metaboExp1}
\docType{data}
\title{
MetaboCombineR experiment example no. 1
}
\description{
This is the first example of LC/MS dataset in metaboAnalyst format.
}
\usage{data("metaboExp1")}
\format{
A data frame with 2288 observations on the following 6 variables.
\describe{
\item{\code{X064.EPK83_m_Mzb1_ESI.mzML}}{a numeric vector    sample no. 1}
\item{\code{X064.EPK88_m_Mzb1_ESI.mzML}}{a numeric vector    sample no. 2}
\item{\code{X064.EPK94_m_Mzb1_ESI.mzML}}{a numeric vector    sample no. 3}
\item{\code{X64Z.61012_m_WT_ESI.mzML}}{a numeric vector    sample no. 4}
\item{\code{X64Z.61014_m_WT_ESI.mzML}}{a numeric vector    sample no. 5}
\item{\code{X64Z.61015_m_WT_ESI.mzML}}{a numeric vector    sample no. 6}
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
data(metaboExp1)
}
\keyword{datasets}
