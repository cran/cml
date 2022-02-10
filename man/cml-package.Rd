\name{cml-package}
\alias{cml-package}
\alias{cml}
\title{
Conditional Manifold Learning
}
\description{
Find a low-dimensional embedding of high-dimensional data, conditioning on auxiliary manifold information. The current version supports conditional MDS and conditional ISOMAP.

Please cite this package as follows:
Bui, A.T. (2021). Conditional Manifold Learning. arXiv.
}
\details{
Brief descriptions of the main functions of the package are provided below:

\code{condMDS()}: is the conditional MDS method, which uses conditional SMACOF to optimize its conditional stress function.

\code{condIsomap()}: is the conditional ISOMAP method, which is basically conditional MDS applying to graph distances (i.e., estimated geodesic distances) of the given distances/dissimilarities.
}

\references{
Bui, A.T. (2021). Conditional Manifold Learning. arXiv:2111.13646 [stat.ML]. https://arxiv.org/abs/2111.13646
}
\author{
Anh Tuan Bui

Maintainer: Anh Tuan Bui <atbui@u.northwestern.edu>
}
\examples{
# generate toy data
N = 100
u.dim = 1
v.dim = 2
u.true = matrix(runif(N*u.dim), N, u.dim)
v.true = matrix(runif(N*v.dim), N, v.dim)
noise = rnorm(N,0,.1)
d = dist(cbind(u.true, v.true, noise))

# Conditional MDS
u.cmds = condMDS(d, v.true, 1)$U

# Conditional ISOMAP
u.cisomap = condIsomap(d, v.true, 1, k = 20)$U
}
