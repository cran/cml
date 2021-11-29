condDist <- function(U, V.tilda, one_n_t=t(rep(1,nrow(U)))) {
  du2 <- dist(U)^2
  VV <- V.tilda %*% t(V.tilda)
  tmp <- diag(VV) %*% one_n_t
  dv2 <- as.dist(tmp + t(tmp) - 2*VV)
  sqrt(du2 + dv2)
}