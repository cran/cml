\name{condDist}
\alias{condDist}
\alias{condDist2}
\title{
Conditional Euclidean distance
}
\description{
Internal functions.
}
\usage{
condDist(U, V.tilde, one_n_t=t(rep(1,nrow(U))))
condDist2(U, V.tilde2, one_n_t=t(rep(1,nrow(U))))
}
\arguments{
  \item{U}{ the embedding \code{U}}
  \item{V.tilde}{ \code{= V \%*\% B}}
  \item{V.tilde2}{ \code{= V \%*\% b^2*t(V)}}
  \item{one_n_t}{ \code{= t(rep(1,nrow(U)))}}
}
\value{
a \code{dist} object.
}
\references{
Bui, A.T. (2024). Dimension Reduction with Prior Information for Knowledge Discovery. IEEE Transactions on Pattern Analysis and Machine Intelligence, 46, 3625-3636. https://doi.org/10.1109/TPAMI.2023.3346212
}
\author{
Anh Tuan Bui
}
