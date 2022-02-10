condIsomap <- function (d, V, u.dim, epsilon = NULL, k, W = NULL,
                        method = c('matrix', 'vector'),
                        it.max = 1000, gamma = 1e-05,
                        init = c('none', 'user'),
                        U.start = NULL, B.start = NULL, b.start = NULL, ...)
{
  if (is.null(epsilon)) {
    d <- vegan::isomapdist(d, k=k,...)
  } else {
    d <- vegan::isomapdist(d, epsilon=epsilon,...)
  }

  V <- as.matrix(V)
  q <- dim(V)[2]
  if (q > 1) {
    method <- match.arg(method, c('matrix', 'vector'))
  } else if (q == 1) {
    method <- 'vector'
  } else {
    stop('V should be a matrix!')
  }

  if (method == 'matrix') {
    return(condSmacof_matrix(d, V, u.dim, W = W, it.max=it.max, gamma=gamma,
                             init = init, U.start = U.start, B.start = B.start))
  } else if (method == 'vector') {
    return(condSmacof_vector(d, V, u.dim, W = W, it.max=it.max, gamma=gamma,
                             init = init, U.start=U.start, b.start = b.start))
  }
}
