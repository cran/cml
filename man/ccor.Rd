\name{ccor}
\alias{ccor}
\title{
Canonical Correlation
}
\description{
Find the canonical correlations of two multivariate datasets.
}
\usage{
ccor(x, y)
}
\arguments{
  \item{x}{ the first dataset.}
  \item{y}{ the second dataset.}
}
\value{
\item{cancor }{ canonical correlations.}
\item{xcoef }{ the columns of this matrix contain the coefficients for the first dataset.}
\item{ycoef }{the columns of this matrix contain the coefficients for the second dataset.}}
\author{
Anh Tuan Bui
}
\examples{
ccor(iris[,1:2], iris[,3:4])
}
