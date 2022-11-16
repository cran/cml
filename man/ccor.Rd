\name{ccor}
\alias{ccor}
\title{
Canonical Correlations
}
\description{
Computes canonical correlations for two sets of multivariate data \code{x} and \code{y}.
}
\usage{
ccor(x, y)
}
\arguments{
  \item{x}{
the first multivariate dataset.
}
  \item{y}{
the second multivariate dataset.
}
}
\value{
a list of the following components:
\item{cancor }{a vector of canonical correlations.}
\item{xcoef }{a matrix, each column of which is the vector of coefficients of x to produce the corresponding canonical covariate.}
\item{ycoef }{a matrix, each column of which is the vector of coefficients of y to produce the corresponding canonical covariate.}
}
\author{
Anh Tuan Bui
}
\examples{
ccor(iris[,1:2], iris[,3:4])
}
