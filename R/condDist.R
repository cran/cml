condDist <- function(U, V.tilde, one_n_t=t(rep(1,nrow(U)))) {
  du2 <- dist(U)^2
  V.tilde2 <- V.tilde %*% t(V.tilde)
  tmp <- diag(V.tilde2) %*% one_n_t
  dv2 <- as.dist(tmp + t(tmp) - 2*V.tilde2)
  sqrt(du2 + dv2)
}
