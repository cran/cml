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
condDist(U, V.tilda, one_n_t=t(rep(1,nrow(U))))
condDist2(U, V.tilda2, one_n_t=t(rep(1,nrow(U))))
}
\arguments{
  \item{U}{ the embedding \code{U}}
  \item{V.tilda}{ \code{= V \%*\% B}}
  \item{V.tilda2}{ \code{= V \%*\% b^2*t(V)}}
  \item{one_n_t}{ \code{= t(rep(1,nrow(U)))}}
}
\value{
a \code{dist} object.
}
\references{
Bui, A.T. (2021). Dimension Reduction with Prior Information for Knowledge Discovery. arXiv:2111.13646. https://arxiv.org/abs/2111.13646.
}
\author{
Anh Tuan Bui
}
