\name{metaboExp3}
\alias{metaboExp3}
\docType{data}
\title{
MetaboCombineR experiment example no. 3
}
\description{
This is the third example of LC/MS dataset in metaboAnalyst format.
}
\usage{data("metaboExp3")}
\format{
A data frame with 2578 observations on the following 6 variables.
\describe{
\item{\code{X064.EPL09_f_Mzb1_ESI.mzML}}{a numeric vector    sample no. 13}
\item{\code{X064.EPL10_m_Mzb1_ESI.mzML}}{a numeric vector    sample no. 14}
\item{\code{X064.EPL12_m_Mzb1_ESI.mzML}}{a numeric vector    sample no. 15}
\item{\code{X64Z.61023_f_WT_ESI.mzML}}{a numeric vector    sample no. 16}
\item{\code{X64Z.61024_f_WT_ESI.mzML}}{a numeric vector    sample no. 17}
\item{\code{X64Z.61025_f_WT_ESI.mzML}}{a numeric vector    sample no. 18}
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
data(metaboExp3)
}
\keyword{datasets}
