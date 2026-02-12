\name{condSmacof_mer}
\alias{condSmacof_mer}
\title{
Conditional SMACOF with incomplete conditioning data (missing entire row)
}
\description{
Conditional SMACOF with missing known feature values. If a row has a missing value in a known feature, the algorithm treats all other known feature values of that row as missing. The algorithm also imputes the actual missing value(s). Intended for internal usage.
}
\usage{
condSmacof_mer(d, V, u.dim, W,
           method = c('matrix', 'vector'), exact = TRUE,
           it.max = 1000, gamma = 1e-05,
           init = c('none', 'eigen', 'user'),
           U.start, B.start, V.tilde.start)
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
\item{W}{
an \emph{N}x\emph{N} symmetric weight matrix. If not given, a matrix of ones will be used.
}
  \item{method}{
if \code{matrix}, there are no restrictions for the B matrix . If \code{vector}, the B matrix is restricted to be diagonal. The latter is more efficient for large \emph{q}. \code{method = 'vector'} is not available in this version.
}
  \item{exact}{
If \code{FALSE}, use the large-\emph{N} approximation formula to update the embedding.
}
  \item{it.max}{
the max number of iterations.
}
  \item{gamma}{
stops early if the reduction of normalized conditional stress is less than \code{gamma}
}
  \item{init}{
initialization method. Except for \code{init = 'user'}, in which case the user directly provides starting values, the other cases apply condMDS with \code{init = 'none'}) (if \code{init = 'none'}) or condMDSeigen (if \code{init = 'eigen'}) to the complete data to initialize \code{B} and part of the embedding.
}
  \item{U.start}{
user-defined starting values for the embedding (when \code{init = 'user'}).
}
  \item{B.start}{
user-defined starting values for the \code{B} matrix (when \code{init = 'user'}).
}
  \item{V.tilde.start}{
user-defined starting values for the \code{V.tilde} matrix (when \code{init = 'user'}).
}
}
\value{
\item{U }{the embedding result.}
\item{B }{the estimated \code{B} matrix.}
\item{stress }{normalized conditional stress value.}
\item{sigma }{the conditional stress value at each iteration.}
\item{init }{the value of the \code{init} argument.}
\item{U.start }{the starting values for the embedding.}
\item{B.start }{the starting values for the \code{B} matrix.}
\item{V.tilde.start }{the starting values for \code{V.tilde.start}.}
\item{V.hat }{imputed \code{V} matrix, if \code{V} contains rows with missing values.}
\item{method }{the value of the \code{method} argument.}
}
\references{
Bui, A.T. (2024). Dimension Reduction with Prior Information for Knowledge Discovery. IEEE Transactions on Pattern Analysis and Machine Intelligence, 46, 3625-3636. https://doi.org/10.1109/TPAMI.2023.3346212

Bui, A. T. (2022). A Closed-Form Solution for Conditional Multidimensional Scaling. Pattern Recognition Letters 164, 148-152. https://doi.org/10.1016/j.patrec.2022.11.007

Bui, A.T. (2026). Conditional Multidimensional Scaling for Incomplete Conditioning Data. Journal of Multivariate Analysis. Accepted.
}
\author{
Anh Tuan Bui
}
