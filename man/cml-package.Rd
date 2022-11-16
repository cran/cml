\name{cml-package}
\alias{cml-package}
\alias{cml}
\title{
Conditional Manifold Learning
}
\description{
Finds a low-dimensional embedding of high-dimensional data, conditioning on available manifold information. The current version supports conditional MDS (based on either conditional SMACOF or closed-form solution) and conditional ISOMAP.

Please cite this package as follows:

Bui, A.T. (2021). Dimension Reduction with Prior Information for Knowledge Discovery. arXiv:2111.13646. https://arxiv.org/abs/2111.13646

Bui, A. T. (2022). A Closed-Form Solution for Conditional Multidimensional Scaling. Pattern Recognition Letters. https://doi.org/10.1016/j.patrec.2022.11.007
}
\details{
Brief descriptions of the main functions of the package are provided below:

\code{condMDS()}: is the conditional MDS method, which uses conditional SMACOF to optimize its conditional stress objective function.

\code{condMDSeigen()}: is the conditional MDS method, which uses a closed-form solution based on multiple linear regression and eigendecomposition.

\code{condIsomap()}: is the conditional ISOMAP method, which is basically conditional MDS applying to graph distances (i.e., estimated geodesic distances) of the given distances/dissimilarities.
}

\references{
Bui, A.T. (2021). Dimension Reduction with Prior Information for Knowledge Discovery. arXiv:2111.13646. https://arxiv.org/abs/2111.13646.

Bui, A. T. (2022). A Closed-Form Solution for Conditional Multidimensional Scaling. Pattern Recognition Letters. https://doi.org/10.1016/j.patrec.2022.11.007

}
\author{
Anh Tuan Bui

Maintainer: Anh Tuan Bui <atbui@u.northwestern.edu>
}
\examples{
## Generate car-brand perception data
factor.weights <- c(90, 88, 83, 82, 81, 70, 68)/562
N <- 100
set.seed(1)
data <- matrix(runif(N*7), N, 7)
colnames(data) <- c('Quality', 'Safety', 'Value',	'Performance', 'Eco', 'Design', 'Tech')
rownames(data) <- paste('Brand', 1:N)
data.hat <- data + matrix(rnorm(N*7), N, 7)*data*.05
data.weighted <- t(apply(data, 1, function(x) x*factor.weights))
d <- dist(data.weighted)
d.hat <- d + rnorm(length(d))*d*.05

## The following examples use the first 4 factors as known features
# Conditional MDS based on conditional SMACOF
u.cmds = condMDS(d.hat, data.hat[,1:4], 3, init='none')
u.cmds$B # compare with diag(factor.weights[1:4])
ccor(data.hat[,5:7], u.cmds$U)$cancor # canonical correlations
vegan::procrustes(data.hat[,5:7], u.cmds$U, symmetric = TRUE)$ss # Procrustes statistic

# Conditional MDS based on the closed-form solution
u.cmds = condMDSeigen(d.hat, data.hat[,1:4], 3)
u.cmds$B # compare with diag(factor.weights[1:4])
ccor(data.hat[,5:7], u.cmds$U)$cancor # canonical correlations
vegan::procrustes(data.hat[,5:7], u.cmds$U, symmetric = TRUE)$ss # Procrustes statistic

# Conditional MDS based on conditional SMACOF,
# initialized by the closed-form solution
u.cmds = condMDS(d.hat, data.hat[,1:4], 3, init='eigen')
u.cmds$B # compare with diag(factor.weights[1:4])
ccor(data.hat[,5:7], u.cmds$U)$cancor # canonical correlations
vegan::procrustes(data.hat[,5:7], u.cmds$U, symmetric = TRUE)$ss # Procrustes statistic

# Conditional ISOMAP
u.cisomap = condIsomap(d.hat, data.hat[,1:4], 3, k = 20, init='eigen')
u.cisomap$B # compare with diag(factor.weights[1:4])
ccor(data.hat[,5:7], u.cisomap$U)$cancor
vegan::procrustes(data.hat[,5:7], u.cisomap$U, symmetric = TRUE)$ss
}
