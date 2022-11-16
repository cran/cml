mpinv <- function (A, eps = sqrt(.Machine$double.eps))
{
  A <- as.matrix(A)
  A.dim <- dim(A)
  A.svd <- svd(A)
  tmp <- sum(A.svd$d > eps*abs(A.svd$d[1]))
  if (tmp <= 0) {
    matrix(0, A.dim[2], A.dim[1])
  } else {
    A.svd$v[, 1:tmp] %*% (t(A.svd$u[,1:tmp])/A.svd$d[1:tmp])
  }
}
