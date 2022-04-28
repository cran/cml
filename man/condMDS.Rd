\name{condMDS}
\alias{condMDS}

\title{
Conditional Multidimensional Scaling
}

\description{
Finds a low-dimensional embedding of a given distance/dissimilarity matrix, conditioning on auxiliary manifold parameters. The method uses conditional SMACOF to optimize its conditional stress function.
}

\usage{
condMDS(d, V, u.dim, W = NULL, method = c('matrix', 'vector'),
        it.max = 1000, gamma = 1e-05,
        init = c('none', 'user'), U.start = NULL,
        B.start = NULL, b.start = NULL) }

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
\item{W}{
an \emph{N}x\emph{N} symmetric weight matrix. If not given, a matrix of ones will be used.
}
  \item{method}{
if \code{matrix}, conditional SMACOF with a full B matrix will be used. If \code{vector}, conditional SMACOF with a diagonal B matrix will be used. The latter is more efficient for large \emph{q}.
}
  \item{it.max}{
the max number of conditional SMACOF iterations.
}
  \item{gamma}{
conditional SMACOF stops early if the reduction of normalized conditional stress is less than \code{gamma}
}
  \item{init}{
method to initialize the starting values for the embedding and B (or b).
}
  \item{U.start}{
user-defined starting values for the embedding (when \code{U.init = 'user'})
}
  \item{B.start}{
starting \code{B} matrix.
}
  \item{b.start}{starting diagonal values of the \code{B} matrix (which is restricted to be diagonal).
}
}

\value{
\item{U }{the embedding result.}
\item{B }{the estimated \code{B} matrix.}
\item{stress }{Normalized conditional stress value.}
\item{sigma }{the conditional stress value at each iteration.}
\item{init }{method to initialize the starting values for the embedding and B (or b).}
\item{U.start }{the starting values for the embedding.}
\item{B.start }{starting values for the \code{B} matrix, if \code{method="matrix"}.}
\item{b.start }{the starting values of the \code{B} matrix (which is restricted to be diagonal), if \code{method="vector"}.}
}
\references{
Bui, A.T. (2021). Dimension Reduction with Prior Information for Knowledge Discovery. arXiv:2111.13646. URL: https://arxiv.org/abs/2111.13646.
}
\author{
Anh Tuan Bui
}

\seealso{
\code{\link{condIsomap}}
}
\examples{
# see help(cml)
}
