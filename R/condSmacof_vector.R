condSmacof_vector <- function(d, V, u.dim, W = NULL, it.max = 1000, gamma = 1e-03,
                              U.start = NULL, b.start = NULL)
{
  p <- u.dim
  q <- dim(V)[2]
  N <- nrow(as.matrix(d))
  if (p > (N - 1 - q))
    stop("Max u.dim is N-1-vdim!")

  tV <- t(V)
  one_n_t <- t(rep(1,N))

  if (is.null(W)) {
    W <- matrix(1, N, N)
    W[is.na(as.matrix(d))] <- 0
  }
  w <- as.dist(W)
  H <- diag(rowSums(W)) - W
  g <- diag(x = tV %*% H %*% V)

  if (is.null(U.start)) {
    U.start <- matrix(runif(N * p, min = -1), ncol = p)
  }
  U <- U.start

  if (is.null(b.start))
    b <- rep(1,q)

  b2 <- b^2
  V.tilda2 <- V %*% (b2*tV)

  eta.d <- sum(w*d^2)
  dz <- condDist2(U, V.tilda2, one_n_t)
  sigma <- rep(0, it.max)
  sigma[1] <- eta.d +
    sum(diag(t(U)%*%H%*%U)) + sum(b2*g) -
    2*sum(w*d*dz)

  if (is.null(W)) {
    for (iter in 2:it.max) {
      cz <- -w*d/dz
      cz[which(cz == Inf)] <- 0
      cz <- as.matrix(cz)
      cz_diag <- -rowSums(cz)
      cz <- cz + diag(cz_diag)
      U <- cz %*% U / N
      t = diag(tV %*% cz %*% V)
      b <- t*b/g
      b2 <- b^2
      V.tilda2 <- V %*% (b2*tV)
      dz <- condDist2(U, V.tilda2, one_n_t)
      sigma[iter] <- eta.d +
        sum(diag(t(U)%*%H%*%U)) + sum(b2*g) -
        2*sum(w*d*dz)
      if (sigma[iter - 1] - sigma[iter] < gamma)
        break()
    }
  } else {
    Hp <- solve(H + 1) - N^(-2)
    for (iter in 2:it.max) {
      cz <- -w*d/dz
      cz[which(cz == Inf)] <- 0
      cz <- as.matrix(cz)
      cz_diag <- -rowSums(cz)
      cz <- cz + diag(cz_diag)
      U <- Hp %*% cz %*% U
      t = diag(tV %*% cz %*% V)
      b <- t*b/g
      b2 <- b^2
      V.tilda2 <- V %*% (b2*tV)
      dz <- condDist2(U, V.tilda2, one_n_t)
      sigma[iter] <- eta.d +
        sum(diag(t(U)%*%H%*%U)) + sum(b2*g) -
        2*sum(w*d*dz)
      if (sigma[iter - 1] - sigma[iter] < gamma)
        break()
    }
  }

  if (iter == it.max) {
    warning("Max iteration reached!")
  } else {
    sigma <- sigma[1:iter]
  }

  list(U = U, b = b, sigma = sigma, U.start = U.start, b.start = b.start)
}


