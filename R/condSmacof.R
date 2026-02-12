condSmacof <- function (d, V, u.dim, W,
                        method = c('matrix', 'vector'), exact = TRUE,
                        it.max = 1000, gamma = 1e-05,
                        init = c('none', 'eigen', 'user'),
                        U.start, B.start)
{
  d <- as.dist(d)
  V <- as.matrix(V)
  Vt <- t(V)
  N <- nrow(V)
  q <- ncol(V)
  if (q > 1) {
    method <- match.arg(method, c('matrix', 'vector'))
  } else if (q == 1) {
    method <- 'vector'
  }
  if (u.dim > (N - 1 - q))
    stop("Max u.dim is N-1-ncol(V)!")

  if (missing(W)) {
    W <- matrix(1, N, N)
    general <- FALSE
  } else {
    W <- (W +t(W))/2
    general <- TRUE
  }

  D <- as.matrix(d)
  D.na <- is.na(D)
  if (any(D.na)) {
    W[D.na] <- 0
    general <- TRUE
  }

  w <- as.dist(W)

  init <- match.arg(init, c('none', 'eigen', 'user'))
  if (init == 'none') {
    B.start <- diag(1, q, q)
    U.start <- matrix(runif(N * u.dim, min = -1), ncol = u.dim)
  } else if (init == 'eigen') {
    tmp <- condMDSeigen(d, V, u.dim, method = method)
    B.start <- matrix(tmp$B, q, q)
    U.start <- matrix(tmp$U, N, u.dim)
  } else if (init == 'user') {
    B.start <- matrix(B.start, q, q)
    U.start <- matrix(U.start, N, u.dim)
  }

  U <- U.start
  B <- B.start

  eta <- sum(w*d^2)
  gamma <- gamma*eta
  one_n_t <- t(rep(1,N))
  sigma <- rep(0, it.max)

  if (general) { # general weights
    H <- diag(rowSums(W)) - W
    Hp <- solve(H + 1) - N^(-2)

    if (method == 'matrix') {
      SpVt <- mpinv(Vt %*% H %*% V) %*% Vt
      V.tilde <- V %*% B
      dz <- condDist(U, V.tilde, one_n_t)
      sigma[1] <- sum(w*(d-dz)^2)

      for (iter in 2:it.max) {
        cz1 <- cz(w, d, dz)
        U <- Hp %*% (cz1 %*% U)
        B <- SpVt %*% cz1 %*% V.tilde
        V.tilde <- V %*% B
        dz <- condDist(U, V.tilde, one_n_t)
        sigma[iter] <- sum(w*(d-dz)^2)
        if (sigma[iter - 1] - sigma[iter] < gamma)
          break()
      }
    } else if (method == 'vector') {
      b <- diag(B)
      s <- diag(x = Vt %*% H %*% V)
      dz <- condDist2(U, V %*% (b^2*Vt), one_n_t)
      sigma[1] <- sum(w*(d-dz)^2)

      for (iter in 2:it.max) {
        cz1 <- cz(w, d, dz)
        U <- Hp %*% (cz1 %*% U)
        b <- diag(Vt %*% cz1 %*% V)*b/s
        dz <- condDist2(U, V %*% (b^2*Vt), one_n_t)
        sigma[iter] <- sum(w*(d-dz)^2)
        if (sigma[iter - 1] - sigma[iter] < gamma) {
          break()
        }
      }
      B <- diag(b, q, q)
    }
  } else { # weights of 1
    Hp <- (diag(1, N, N) - matrix(1/N, N, N))/N
    vs <- colSums(V)
    Vip <- solve(Vt %*% V)
    if (method == 'matrix') {
      Vipvso <- Vip %*% (vs %*% t(vs))
      SpVt <- ((diag(1, q, q) + Vipvso/(N - sum(diag(Vipvso)))) %*% Vip/N) %*% Vt
      V.tilde <- V %*% B
      dz <- condDist(U, V.tilde, one_n_t)
      us <- colSums(U)
      vs.tilde <- colSums(V.tilde)
      sigma[1] <- sum(w*(d-dz)^2)

      if (exact) {
        for (iter in 2:it.max) {
          cz1 <- cz(w, d, dz)
          U <- Hp %*% (cz1 %*% U)
          B <- SpVt %*% cz1 %*% V.tilde
          V.tilde <- V %*% B
          dz <- condDist(U, V.tilde, one_n_t)
          sigma[iter] <- sum(w*(d-dz)^2)
          if (sigma[iter - 1] - sigma[iter] < gamma)
            break()
        }
      } else {
        for (iter in 2:it.max) {
          cz1 <- cz(w, d, dz)
          U <- cz1 %*% U / N
          B <- SpVt %*% cz1 %*% V.tilde
          V.tilde <- V %*% B
          dz <- condDist(U, V.tilde, one_n_t)
          sigma[iter] <- sum(w*(d-dz)^2)
          if (sigma[iter - 1] - sigma[iter] < gamma)
            break()
        }
      }
    } else if (method == 'vector') {
      b <- diag(B)
      s <- N*colSums(V^2) - vs^2
      dz <- condDist2(U, V %*% (b^2*Vt), one_n_t)
      sigma[1] <- sum(w*(d-dz)^2)

      if (exact) {
        for (iter in 2:it.max) {
          cz1 <- cz(w, d, dz)
          U <- Hp %*% (cz1 %*% U)
          b <- diag(Vt %*% cz1 %*% V)*b/s
          dz <- condDist2(U, V %*% (b^2*Vt), one_n_t)
          sigma[iter] <- sum(w*(d-dz)^2)
          if (sigma[iter - 1] - sigma[iter] < gamma) {
            break()
          }
        }
      } else {
        for (iter in 2:it.max) {
          cz1 <- cz(w, d, dz)
          U <- cz1 %*% U / N
          b <- diag(Vt %*% cz1 %*% V)*b/s
          dz <- condDist2(U, V %*% (b^2*Vt), one_n_t)
          sigma[iter] <- sum(w*(d-dz)^2)
          if (sigma[iter - 1] - sigma[iter] < gamma) {
            break()
          }
        }
      }
      B <- diag(b, q, q)
    }
  }

  if (iter == it.max) {
    warning("Max iteration reached!")
  } else {
    sigma <- sigma[1:iter]
  }

  list(U = U, B = B, stress = sigma[iter]/eta, sigma = sigma, init = init,
       U.start = U.start, B.start = B.start, method = method, exact = exact)
}
