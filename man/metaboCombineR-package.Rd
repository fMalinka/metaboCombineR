\name{metaboCombineR-package}
\alias{metaboCombineR-package}
\alias{metaboCombineR}
\docType{package}
\title{
metaboCombineR: A package for aligning LC/MS experiments
}
\description{
A more detailed description of what the package does. A length
of about one to five lines is recommended.
}
\details{
This section should provide a more detailed overview of how to use the
package, including the most important functions.
}
\author{
Your Frantisek Malinka, frantisek.malinka@gmail.com.
Maintainer: Frantisek Malinka <frantisek.malinka@gmail.com>
}
\references{
This optional section can contain literature or other references for
background information.
}
\keyword{ package }
\seealso{
Optional links to other man pages
}
\examples{
data(metaboExp1)
data(metaboExp2)
data(metaboExp3)
data(metaboExp4)
runMetaboCombiner(list(metaboExp1, metaboExp2, metaboExp3, metaboExp4),
mzprecision=3)
}
