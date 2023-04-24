\name{condMDSeigen}
\alias{condMDSeigen}
\title{
Conditional Multidimensional Scaling With Closed-Form Solution
}
\description{
Provides a closed-form solution for conditional multidimensional scaling, based on multiple linear regression and eigendecomposition.}
\usage{
condMDSeigen(d, V, u.dim, method = c('matrix', 'vector'))
}
\arguments{
  \item{d}{
a \code{dist} object of \emph{N} entities.
}
  \item{V}{
an \emph{N}x\emph{q} matrix of \emph{q} manifold auxiliary parameter values of the \emph{N} entities.
}
  \item{u.dim}{
the embedding dimension.
}
  \item{method}{
if \code{matrix}, there are no restrictions for the B matrix . If \code{vector}, the B matrix is restricted to be diagonal.
}
}
\value{
\item{U }{the embedding result.}
\item{B }{the estimated \code{B} matrix.}
\item{eig }{the computed eigenvalues.}
\item{stress }{the corresponding normalized conditional stress value of the solution.}
}
\references{
Bui, A. T. (2022). A Closed-Form Solution for Conditional Multidimensional Scaling. Pattern Recognition Letters 164, 148-152. https://doi.org/10.1016/j.patrec.2022.11.007
}
\author{
Anh Tuan Bui
}
\seealso{
\code{\link{condMDS}}, \code{\link{condIsomap}}
}
\examples{
# see help(cml)
}
