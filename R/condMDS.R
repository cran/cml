condMDS <- function (d, V, u.dim, W,
                     method = c('matrix', 'vector'), exact = TRUE,
                     it.max = 1000, gamma = 1e-05,
                     init = c('none', 'eigen', 'user'),
                     U.start, B.start, V.tilde.start = NULL)
{
  if (any(is.na(V))) {
    cat("V contains missing values and has been imputed (see <output_name>$V.hat).\n\n")
    return(condSmacof_mer(d, V, u.dim, W = W, method = method, exact = exact,
                          it.max = it.max, gamma = gamma, init = init,
                          U.start = U.start, B.start = B.start,
                          V.tilde.start = V.tilde.start))
  } else {
    return(condSmacof(d, V, u.dim, W = W, method = method, exact = exact,
                      it.max = it.max, gamma = gamma, init = init,
                      U.start = U.start, B.start = B.start))
  }
}
