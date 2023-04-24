ccor <- function(x, y) {

  x = as.matrix(x)
  y = as.matrix(y)
  d <- min(ncol(x),ncol(y))

  cov.x.inv = solve(cov(x))
  cov.y.inv = solve(cov(y))
  cov.xy = cov(x,y)
  cov.yx = t(cov.xy)
  x.eig = eigen(cov.x.inv %*% cov.xy %*% cov.y.inv %*% cov.yx)
  y.eig = eigen(cov.y.inv %*% cov.yx %*% cov.x.inv %*% cov.xy)

  list(cancor = sqrt(x.eig$values[1:d]),
       xcoef = x.eig$vectors[,1:d],
       ycoef = y.eig$vectors[,1:d])
}
