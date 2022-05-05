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

  condMDS(d = d, V = V, u.dim = u.dim, W = W,
          method = method, it.max = it.max, gamma = gamma,
          init = init, U.start = U.start,
          B.start = B.start, b.start = b.start)
}
