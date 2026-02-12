\name{cml-package}
\alias{cml-package}
\alias{cml}
\title{
Conditional Manifold Learning
}
\description{
Find a low-dimensional embedding of high-dimensional data, conditioning on auxiliary manifold information. The current version supports conditional MDS and conditional ISOMAP.

Please cite this package as follows:
Bui, A.T. (2024). Dimension Reduction with Prior Information for Knowledge Discovery. IEEE Transactions on Pattern Analysis and Machine Intelligence, 46, 3625-3636. https://doi.org/10.1109/TPAMI.2023.3346212

}
\details{
Brief descriptions of the main functions of the package are provided below:

\code{condMDS()}: is the conditional MDS method, which uses conditional SMACOF to optimize its conditional stress function.

\code{condIsomap()}: is the conditional ISOMAP method, which is basically conditional MDS applying to graph distances (i.e., estimated geodesic distances) of the given distances/dissimilarities.
}

\references{
Bui, A.T. (2024). Dimension Reduction with Prior Information for Knowledge Discovery. IEEE Transactions on Pattern Analysis and Machine Intelligence, 46, 3625-3636. https://doi.org/10.1109/TPAMI.2023.3346212

Bui, A. T. (2022). A Closed-Form Solution for Conditional Multidimensional Scaling. Pattern Recognition Letters 164, 148-152. https://doi.org/10.1016/j.patrec.2022.11.007

Bui, A.T. (2026). Conditional Multidimensional Scaling for Incomplete Conditioning Data. Journal of Multivariate Analysis. Accepted.
}
\author{
Anh Tuan Bui

Maintainer: Anh Tuan Bui <atbui@u.northwestern.edu>
}
\examples{
# Generate car-brand perception data
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

# Conditional MDS, using the first 4 factors as known features
u.cmds = condMDS(d.hat, data.hat[,1:4], 3, init='eigen')
u.cmds$B # compare with diag(factor.weights[1:4])
ccor(data.hat[,5:7], u.cmds$U)$cancor # canonical correlations
vegan::procrustes(data.hat[,5:7], u.cmds$U, symmetric = TRUE)$ss # Procrustes statistic

# Conditional ISOMAP, using the first 4 factors as known features
u.cisomap = condIsomap(d.hat, data.hat[,1:4], 3, k = 20, init='eigen')
u.cisomap$B # compare with diag(factor.weights[1:4])
ccor(data.hat[,5:7], u.cisomap$U)$cancor
vegan::procrustes(data.hat[,5:7], u.cisomap$U, symmetric = TRUE)$ss

# with missing values in V
V = data.hat[,1:4]
V[1, 1] = NA
u.cmds2 = condMDS(d.hat, V, 3, init='eigen')
u.cmds2$B # compare with diag(factor.weights[1:4])
ccor(data.hat[,5:7], u.cmds2$U)$cancor # canonical correlations
vegan::procrustes(data.hat[,5:7], u.cmds$U, symmetric = TRUE)$ss # Procrustes statistic
u.cmds2$V.hat[1,1] # imputed value for V[1, 1]; the ground truth is data[1,1]
}
