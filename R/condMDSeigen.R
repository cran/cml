condMDSeigen <- function(d, V, u.dim, method = c('matrix', 'vector'))
{
  d <- as.dist(d)
  V <- as.matrix(V)
  q <- ncol(V)
  N <- nrow(V)
  if (u.dim > (N - 1 - q))
    stop("Max u.dim is N-1-ncol(V)!")

  # estimate B
  method <- match.arg(method, c('matrix', 'vector'))
  if (q > 1 & method == 'matrix') {
    v.eig <- eigen(cov(V))
    R <- v.eig$vectors %*% diag((v.eig$values)^(-.5)) %*% t(v.eig$vectors)
    V.whitened = V %*% R
  } else {
    R = diag(q)
    V.whitened = V
  }
  d2.V <- sapply(1:q, function(k) c(dist(V.whitened[,k])^2))
  d2 <- d^2
  tmp <- lm(c(d2) ~ ., data=as.data.frame(d2.V))

  betas <- tmp[[1]][-1]
  if (any(betas < 0)) {
    cat('Negative beta coefficient(s) found and set to 0!')
    betas[which(betas < 0)] <- 0
  }

  # estimate U
  if (q > 1) {
    B <- R %*% diag(sqrt(betas))
    V.tilda <- V %*% B
  } else {
    B <- sqrt(betas)
    V.tilda <- V * B
  }
  M <- diag(N) - matrix(1/N,N,N)
  E <- eigen(M %*% (as.matrix(-d2/2) - V.tilda %*% t(V.tilda)) %*% M, symmetric = TRUE)
  ev <- E$values[1:u.dim]
  evec <- E$vectors[, 1:u.dim, drop = FALSE]
  k1 <- sum(ev > 0)
  if (k1 < u.dim) {
    warning(gettextf("only %d of the first %d eigenvalues are > 0",
                     k1, u.dim), domain = NA)
    evec <- evec[, 1:k1, drop = FALSE]
    ev <- ev[1:k1]
  }
  U <- evec * rep(sqrt(ev), each = N)

  list(U = U, B = B, eig = ev,
       stress = sum((d - condDist(U, V.tilda))^2)/sum(d2))
}
