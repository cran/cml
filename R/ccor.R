ccor <- function(x, y) {

  x = as.matrix(x)
  y = as.matrix(y)
  cov.x.inv = solve(cov(x))
  cov.y.inv = solve(cov(y))
  cov.xy = cov(x,y)
  cov.yx = t(cov.xy)
  x.eig = eigen(cov.x.inv %*% cov.xy %*% cov.y.inv %*% cov.yx)
  y.eig = eigen(cov.y.inv %*% cov.yx %*% cov.x.inv %*% cov.xy)

  list(cancor = x.eig$values,
       xcoef = x.eig$vectors,
       ycoef = y.eig$vectors)
}
