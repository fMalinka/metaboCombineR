\name{runMetaboCombiner}
\alias{runMetaboCombiner}
\title{Run metaboCombineR algorithm}
\description{
Description of metaboCombiner algorithm.
}
\usage{
runMetaboCombiner(listExperimens, mzprecision, windowsize)
}
\arguments{
\item{listExperimens}{list of experiments}
\item{mzprecision}{number of digits of mz values}
\item{windowsize}{size of windows}
}
\details{
The method simply searches for a given fragment ion in an xcmsFragment
object type given a certain ppm error window
}
\value{
A list frame with the following elements:

\item{colOntologyPath}{Path to DOT column ontology in OBO format}
\item{rowOntologyPath}{Path to Gene ontology row ontology in OBO
format (simplified version)}
\item{datamatrix}{Input binary DOT matrix (subset)}
\item{colOntologyDesc}{List of column anotation}
\item{rowOntologyDesc}{List of row anotation}
}
\seealso{
Look at our readme file at
\url{https://github.com/fMalinka/metaboCombiner}
}
\examples{
myExample <- runMetaboCombiner()
}
\author{Frantisek Malinka}
