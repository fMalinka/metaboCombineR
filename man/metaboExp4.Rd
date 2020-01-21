\name{metaboExp4}
\alias{metaboExp4}
\docType{data}
\title{
MetaboCombineR experiment example no. 4
}
\description{
This is the fourth example of LC/MS dataset in metaboAnalyst format.
}
\usage{data("metaboExp4")}
\format{
A data frame with 2442 observations on the following 6 variables.
\describe{
\item{\code{X064.EPL15_f_Mzb1_ESI.mzML}}{a numeric vector    sample no. 19}
\item{\code{X064.EPL16_f_Mzb1_ESI.mzML}}{a numeric vector    sample no. 20}
\item{\code{X064.EPL34_f_Mzb1_ESI.mzML}}{a numeric vector    sample no. 21}
\item{\code{X64Z.61026_f_WT_ESI.mzML}}{a numeric vector    sample no. 22}
\item{\code{X64Z.61027_f_WT_ESI.mzML}}{a numeric vector    sample no. 23}
\item{\code{X64Z.61028_f_WT_ESI.mzML}}{a numeric vector    sample no. 24}
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
data(metaboExp4)
}
\keyword{datasets}
