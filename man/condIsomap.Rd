\name{condIsomap}
\alias{condIsomap}

\title{
Conditional ISOMAP
}

\description{
Finds a low-dimensional manifold embedding of a given distance/dissimilarity matrix, conditioning on auxiliary manifold parameters. The method applies conditional MDS (see \link{condMDS}) to a graph distance matrix computed for the given distances/dissimilarities, using the \code{isomap{vegan}} function.
}

\usage{
condIsomap(d, V, u.dim, epsilon = NULL, k, W = NULL, method = c('matrix', 'vector'),
           it.max = 1000, gamma = 1e-03, U.start = NULL, B.start = NULL,...) }

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
an \emph{N}x\emph{N} weight matrix. If not given, a matrix of ones will be used.
}
  \item{method}{
if \code{matrix}, conditional SMACOF with a full B matrix will be used. If \code{vector}, conditional SMACOF with a diagonal B matrix will be used. The latter is more efficient for large \emph{q}.
}
  \item{it.max}{
the max number of conditional SMACOF iterations.
}
  \item{gamma}{
conditional SMACOF stops early if the conditional stress reduction is less than \code{gamma}
}
  \item{U.start}{
starting values for the embedding.
}
  \item{B.start}{
starting values for the \code{B} matrix.
}
  \item{...}{
other arguments for the \code{isomap{vegan}} function.
}
}

\value{
\item{U }{the embedding result.}
\item{B }{the estimated \code{B} matrix.}
\item{sigma }{the conditional stress value at each iteration.}
\item{U.start }{the starting values for the embedding.}
\item{B.start }{starting values for the \code{B} matrix, if \code{method="matrix"}.}
\item{b.start }{the starting values of the \code{B} matrix (which is restricted to be diagonal), if \code{method="vector"}.}
}
\references{
Bui, A.T. (2021). Conditional Manifold Learning. arXiv:2111.13646 [stat.ML]. https://arxiv.org/abs/2111.13646
}
\author{
Anh Tuan Bui
}

\seealso{
\code{\link{condMDS}}
}
\examples{
# see help(cml)
}