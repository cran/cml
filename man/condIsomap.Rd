\name{condIsomap}
\alias{condIsomap}

\title{
Conditional ISOMAP
}

\description{
Finds a low-dimensional manifold embedding of a given distance/dissimilarity matrix, conditioning on available manifold information. The method applies conditional MDS (see \link{condMDS}) to a graph distance matrix computed for the given distances/dissimilarities, using the \code{isomap{vegan}} function.
}

\usage{
condIsomap(d, V, u.dim, epsilon = NULL, k, W,
           method = c('matrix', 'vector'), exact = TRUE,
           it.max = 1000, gamma = 1e-05,
           init = c('none', 'eigen', 'user'),
           U.start, B.start, ...) }

\arguments{
  \item{d}{
a distance/dissimilarity matrix of \emph{N} entities (or a \code{dist} object).
}
  \item{V}{
an \emph{N}x\emph{q} matrix of \emph{q} manifold auxiliary parameter values of the \emph{N} entities.
}
  \item{u.dim}{
the embedding dimension.
}
  \item{epsilon}{
shortest dissimilarity retained.
}
  \item{k}{
Number of shortest dissimilarities retained for a point. If both \code{epsilon} and \code{k} are given, \code{epsilon} will be used.
}
\item{W}{
an \emph{N}x\emph{N} symmetric weight matrix. If not given, a matrix of ones will be used.
}
  \item{method}{
if \code{matrix}, there are no restrictions for the B matrix . If \code{vector}, the B matrix is restricted to be diagonal. The latter is more efficient for large \emph{q}.
}
  \item{exact}{
only relevant if \code{W} is not given. In this case, if \code{exact == FALSE}, \code{U} is updated by the large-\emph{N} approximation formula.
}
  \item{it.max}{
the max number of conditional SMACOF iterations.
}
  \item{gamma}{
conditional SMACOF stops early if the reduction of normalized conditional stress is less than \code{gamma}
}
  \item{init}{
initialization method.
}
  \item{U.start}{
user-defined starting values for the embedding (when \code{init = 'user'})
}
  \item{B.start}{
starting \code{B} matrix.
}
  \item{...}{
other arguments for the \code{isomap{vegan}} function.
}
}

\value{
\item{U }{the embedding result.}
\item{B }{the estimated \code{B} matrix.}
\item{stress }{Normalized conditional stress value.}
\item{sigma }{the conditional stress value at each iteration.}
\item{init }{the value of the \code{init} argument.}
\item{U.start }{the starting values for the embedding.}
\item{B.start }{starting values for the \code{B} matrix.}
\item{method }{the value of the \code{method} argument.}
\item{exact }{the value of the \code{exact} argument.}
}
\references{
Bui, A.T. (2021). Dimension Reduction with Prior Information for Knowledge Discovery. arXiv:2111.13646. https://arxiv.org/abs/2111.13646.

Bui, A. T. (2022). A Closed-Form Solution for Conditional Multidimensional Scaling. Pattern Recognition Letters 164, 148-152. https://doi.org/10.1016/j.patrec.2022.11.007
}
\author{
Anh Tuan Bui
}

\seealso{
\code{\link{condMDS}}, \code{\link{condMDSeigen}}
}
\examples{
# see help(cml)
}
