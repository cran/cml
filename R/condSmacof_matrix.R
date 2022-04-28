condSmacof_matrix <- function (d, V, u.dim, W = NULL, it.max = 1000, gamma = 1e-05,
                               init = c('none', 'user'),
                               U.start, B.start)
{
  q <- ncol(V)
  N <- nrow(V)
  if (u.dim > (N - 1 - q))
    stop("Max u.dim is N-1-ncol(V)!")

  one_n_t <- t(rep(1,N))

  init <- match.arg(init, c('none', 'user'))
  if (init == 'none') {
    B.start <- diag(rep(1,q))
    U.start <- matrix(runif(N * u.dim, min = -1), ncol = u.dim)
  } else {
    if (is.null(B.start) | is.null(U.start)) {
      stop('U.start and B.start need to be provided when init = "user"!')
    } else {
      B.start <- as.matrix(B.start)
      U.start <- as.matrix(U.start)
    }
  }
  B <- B.start
  U <- U.start

  if (is.null(W)) {
    W <- matrix(1, N, N)
    W[is.na(as.matrix(d))] <- 0
    flag <- TRUE
  } else {
    flag <- FALSE
  }
  w <- as.dist(W)
  H <- diag(rowSums(W)) - W
  GptV <- ginv(t(V) %*% H %*% V) %*% t(V)

  eta.d <- sum(w*d^2)
  V.tilda <- V %*% B
  dz <- condDist(U, V.tilda, one_n_t)
  sigma <- rep(0, it.max)
  sigma[1] <- eta.d + sum(diag(t(U)%*%H%*%U)) +
    sum(diag(t(V.tilda) %*% H %*% V.tilda)) - 2*sum(w*d*dz)
  gamma <- gamma*eta.d
  if (flag) {
    for (iter in 2:it.max) {
      cz <- -w*d/dz
      cz[which(cz == Inf)] <- 0
      cz <- as.matrix(cz)
      cz_diag <- -rowSums(cz)
      cz <- cz + diag(cz_diag)
      U <- cz %*% U / N
      B <- GptV %*% cz %*% V.tilda
      V.tilda <- V %*% B
      dz <- condDist(U, V.tilda, one_n_t)
      sigma[iter] <- eta.d +
        sum(diag(t(U)%*%H%*%U)) + sum(diag(t(V.tilda) %*% H %*% V.tilda)) -
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
      B <- GptV %*% cz %*% V.tilda
      V.tilda <- V %*% B
      dz <- condDist(U, V.tilda, one_n_t)
      sigma[iter] <- eta.d +
        sum(diag(t(U)%*%H%*%U)) + sum(diag(t(V.tilda) %*% H %*% V.tilda)) -
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

  stress <- sigma[iter]/eta.d

  list(U = U, B = B, stress = stress, sigma = sigma, init = init,
       U.start = U.start, B.start = B.start)
}
