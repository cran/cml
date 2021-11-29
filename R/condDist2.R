condDist2 <- function(U, V.tilda2, one_n_t=t(rep(1,nrow(U)))) {
  du2 <- dist(U)^2
  tmp <- diag(V.tilda2) %*% t(rep(1,nrow(U)))
  dv2 <- as.dist(tmp + t(tmp) - 2*V.tilda2)
  sqrt(du2 + dv2)
}