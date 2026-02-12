condIsomap <- function (d, V, u.dim, epsilon = NULL, k, W,
                        method = c('matrix', 'vector'), exact = TRUE,
                        it.max = 1000, gamma = 1e-05,
                        init = c('none', 'condSmacof', 'eigen', 'user'),
                        U.start, B.start,
                        V.tilde.start = NULL, ...)
{
  if (is.null(epsilon)) {
    d <- vegan::isomapdist(d, k=k,...)
  } else {
    d <- vegan::isomapdist(d, epsilon=epsilon,...)
  }

  condMDS(d = d, V = V, u.dim = u.dim, W = W,
          method = method, exact = exact, it.max = it.max, gamma = gamma,
          init = init, U.start = U.start, B.start = B.start,
          V.tilde.start = V.tilde.start)
}
