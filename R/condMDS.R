condMDS <- function (d, V, u.dim, W = NULL, method = c('matrix', 'vector'),
                     it.max = 1000, gamma = 1e-03, U.start = NULL, B.start = NULL)
{
  d <- as.dist(d)
  V <- as.matrix(V)
  method <- match.arg(method, c('matrix', 'vector'))
  if (method == 'matrix') {
    return(condSmacof_matrix(d, V, u.dim, W = W, it.max = it.max, gamma = gamma,
                             U.start = U.start, B.start = B.start))
  } else if (method == 'vector') {
    return(condSmacof_vector(d, V, u.dim, W = W, it.max = it.max, gamma = gamma,
                             U.start = U.start, b.start = diag(B.start)))
  }
}
