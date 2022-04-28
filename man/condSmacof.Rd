\name{condSmacof}
\alias{condSmacof_matrix}
\alias{condSmacof_vector}
\title{
Conditional SMACOF
}
\description{
Conditional SMACOF algorithms. Intended for internal usage.
}
\usage{
condSmacof_matrix(d, V, u.dim, W = NULL, it.max = 1000, gamma = 1e-05,
                  init = c('none', 'user'), U.start, B.start)
condSmacof_vector(d, V, u.dim, W = NULL, it.max = 1000, gamma = 1e-05,
                  init = c('none', 'user'), U.start, b.start)
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
  \item{it.max}{
the max number of conditional SMACOF iterations.
}
  \item{gamma}{
conditional SMACOF stops early if the  reduction of normalized conditional stress is less than \code{gamma}
}
  \item{init}{
method to initialize the starting values for the embedding and B (or b).
}
  \item{U.start}{
user-defined starting values for the embedding (when \code{U.init = 'user'})
}
  \item{B.start}{\code{condSmacof_matrix()} only; starting \code{B} matrix.
}
  \item{b.start}{for \code{condSmacof_vector()} only; starting diagonal values of the \code{B} matrix (which is restricted to be diagonal).
}
}
\details{
\code{condSmacof_matrix()} uses a full \code{B} matrix, whereas \code{condSmacof_vector()} restricts \code{B} to be a diagonal matrix and uses vector operations on the diagonal elements of \code{B}. The latter is therefore faster when \emph{q} is large.
}
\value{
\item{U }{the embedding result.}
\item{B }{the estimated \code{B} matrix.}
\item{stress }{Normalized conditional stress value.}
\item{sigma }{the conditional stress value at each iteration.}
\item{init }{method to initialize the starting values for the embedding and B (or b).}
\item{U.start }{the starting values for the embedding.}
\item{B.start }{the starting \code{B} matrix, if using \code{condSmacof_matrix()}}
\item{b.start }{the starting diagonal values of the \code{B} matrix (which is restricted to be diagonal), if using \code{condSmacof_vector()}.}

}
\references{
Bui, A.T. (2021). Dimension Reduction with Prior Information for Knowledge Discovery. arXiv:2111.13646. URL: https://arxiv.org/abs/2111.13646.
}
\author{
Anh Tuan Bui
}
