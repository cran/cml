condSmacof <- function (d, V, u.dim, W,
                        method = c('matrix', 'vector'), exact = TRUE,
                        it.max = 1000, gamma = 1e-05,
                        init = c('none', 'eigen', 'user'),
                        U.start, B.start)
{
  d <- as.dist(d)
  V <- as.matrix(V)
  q <- ncol(V)
  if (q > 1) {
    method <- match.arg(method, c('matrix', 'vector'))
  } else if (q == 1) {
    method <- 'matrix'
  }
  N <- nrow(V)
  if (u.dim > (N - 1 - q))
    stop("Max u.dim is N-1-ncol(V)!")

  init <- match.arg(init, c('none', 'eigen', 'user'))
  if (init == 'none') {
    B.start <- diag(rep(1,q))
    U.start <- matrix(runif(N * u.dim, min = -1), ncol = u.dim)
  } else if (init == 'eigen') {
    tmp <- condMDSeigen(d, V, u.dim, method = method)
    B.start <- tmp$B
    U.start <- tmp$U
  } else if (init == 'user') {
    B.start <- as.matrix(B.start)
    U.start <- as.matrix(U.start)
  }

  if (missing(W)) {
    W <- matrix(1, N, N)
  }
  W[is.na(as.matrix(d))] <- 0
  w <- as.dist(W)

  H <- diag(rowSums(W)) - W

  if (any(W!=1)) {
    exact <- TRUE
    Hp <- solve(H + 1) - N^(-2)
  } else {
    if (exact) {
      Hp <- H/(N^2)
    }
  }

  U <- U.start
  eta.d <- sum(w*d^2)
  gamma <- gamma*eta.d
  one_n_t <- t(rep(1,N))

  sigma <- rep(0, it.max)
  if (method == 'matrix') {
    B <- B.start
    V.tilda <- V %*% B
    dz <- condDist(U, V.tilda, one_n_t)
    sigma[1] <- eta.d + sum(diag(t(U)%*%H%*%U)) +
      sum(diag(t(V.tilda) %*% H %*% V.tilda)) - 2*sum(w*d*dz)
    GptV <- mpinv(t(V) %*% H %*% V) %*% t(V)

    if (exact) {
      for (iter in 2:it.max) {
        cz1 <- cz(w, d, dz)
        U <- Hp %*% (cz1 %*% U)
        B <- GptV %*% cz1 %*% V.tilda
        V.tilda <- V %*% B
        dz <- condDist(U, V.tilda, one_n_t)
        sigma[iter] <- eta.d +
          sum(diag(t(U)%*%H%*%U)) + sum(diag(t(V.tilda) %*% H %*% V.tilda)) -
          2*sum(w*d*dz)
        if (sigma[iter - 1] - sigma[iter] < gamma)
          break()
      }
    } else {
      for (iter in 2:it.max) {
        cz1 <- cz(w, d, dz)
        U <- cz1 %*% U / N
        B <- GptV %*% cz1 %*% V.tilda
        V.tilda <- V %*% B
        dz <- condDist(U, V.tilda, one_n_t)
        sigma[iter] <- eta.d +
          sum(diag(t(U)%*%H%*%U)) + sum(diag(t(V.tilda) %*% H %*% V.tilda)) -
          2*sum(w*d*dz)
        if (sigma[iter - 1] - sigma[iter] < gamma)
          break()
      }
    }
  } else if (method == 'vector') {
    b <- diag(B.start)
    tV <- t(V)
    g <- diag(x = tV %*% H %*% V)
    b2 <- b^2
    dz <- condDist2(U, V %*% (b2*tV), one_n_t)
    sigma[1] <- eta.d + sum(diag(t(U)%*%H%*%U)) + sum(b2*g) - 2*sum(w*d*dz)

    if (exact) {
      for (iter in 2:it.max) {
        cz1 <- cz(w, d, dz)
        U <- Hp %*% (cz1 %*% U)
        b <- diag(tV %*% cz1 %*% V)*b/g
        b2 <- b^2
        dz <- condDist2(U, V %*% (b2*tV), one_n_t)
        sigma[iter] <- eta.d + sum(diag(t(U)%*%H%*%U)) + sum(b2*g) - 2*sum(w*d*dz)
        if (sigma[iter - 1] - sigma[iter] < gamma) {
          B <- diag(b)
          break()
        }
      }
    } else {
      for (iter in 2:it.max) {
        cz1 <- cz(w, d, dz)
        U <- cz1 %*% U / N
        b <- diag(tV %*% cz1 %*% V)*b/g
        b2 <- b^2
        dz <- condDist2(U, V %*% (b2*tV), one_n_t)
        sigma[iter] <- eta.d + sum(diag(t(U)%*%H%*%U)) + sum(b2*g) - 2*sum(w*d*dz)
        if (sigma[iter - 1] - sigma[iter] < gamma) {
          B <- diag(b)
          break()
        }
      }
    }
  }

  if (iter == it.max) {
    warning("Max iteration reached!")
  } else {
    sigma <- sigma[1:iter]
  }

  list(U = U, B = B, stress = sigma[iter]/eta.d, sigma = sigma, init = init,
       U.start = U.start, B.start = B.start, method = method, exact = exact)
}
