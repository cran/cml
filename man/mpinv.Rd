\name{mpinv}
\alias{mpinv}
\title{
Moore-Penrose Inverse
}
\description{
Computes the Moore-Penrose inverse (a.k.a., generalized inverse or pseudoinverse) of a matrix based on singular-value decomposition (SVD).
}
\usage{
mpinv(A, eps = sqrt(.Machine$double.eps))
}
\arguments{
  \item{A}{
a matrix of real numbers.
}
  \item{eps}{
a threshold (to be multiplied with the largest singular value) for dropping SVD parts that correspond to small singular values.
}
}
\value{
the Moore-Penrose inverse.
}
\author{
Anh Tuan Bui
}
\examples{
mpinv(2*diag(4))
}
\keyword{generalized inverse}
\keyword{pseudoinverse}
