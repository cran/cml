condSmacof_mer <- function (d, V, u.dim, W,
                            method = c('matrix', 'vector'), exact = TRUE,
                            it.max = 1000, gamma = 1e-05,
                            init = c('none', 'eigen', 'user'),
                            U.start, B.start, V.tilde.start)
{
  d <- as.dist(d)
  V <- as.matrix(V)
  N <- nrow(V)
  q <- ncol(V)
  if (q > 1) {
    method <- match.arg(method, c('matrix', 'vector'))
  } else if (q == 1) {
    method <- 'vector'
  }

  if (u.dim > (N - 1 - q))
    stop("Max u.dim is N-1-ncol(V)!")

  V.sorted <- sort.int(complete.cases(V), decreasing = T, index.return = T)
  V.sorted.rec.id <- sort.int(V.sorted$ix,decreasing = F,index.return = T)$ix
  N1 <- sum(V.sorted$x)
  N2 <- N - N1
  V1 <- matrix(V[V.sorted$ix[1:N1],], N1, q)
  V1t <- t(V1)
  V2 <- matrix(V[V.sorted$ix[(N1+1):N],], N2, q)
  missing.ids <- which(is.na(V2))
  V2[missing.ids] <- 0

  D <- as.matrix(d)[V.sorted$ix,V.sorted$ix]
  d <- as.dist(D)
  d1 <- as.dist(D[1:N1,1:N1])

  if (missing(W)) {
    W <- matrix(1, N, N)
    general <- FALSE
  } else {
    W <- (W +t(W))/2
    W <- W[V.sorted$ix,V.sorted$ix]
    general <- TRUE
  }

  D.na <- is.na(D)
  if (any(D.na)) {
    W[D.na] <- 0
    general <- TRUE
  }

  w <- as.dist(W)

  init <- match.arg(init, c('none', 'eigen', 'user'))
  if (init == 'none') {
    if (u.dim > (N1 - 1 - q)) {
      tmp <- condSmacof(d1, V1, u.dim, W[1:N1,1:N1], method = method,
                        it.max = it.max, gamma = gamma, init = 'none')
      B.start <- matrix(tmp$B, q, q)
      U1.start <- matrix(tmp$U, N1, u.dim)
    } else {
      B.start <- diag(1, q, q)
      U1.start <- matrix(runif(N1 * u.dim, min = -1), ncol = u.dim)
    }
    U2.start <- matrix(runif(N2 * u.dim, min = -1), ncol = u.dim)
    V2.tilde.start <- V2 %*% B.start
  } else if (init == 'eigen') {
    tmp <- condMDSeigen(d1, V1, u.dim, method=method)
    B.start <- matrix(tmp$B, q, q)
    U1.start <- matrix(tmp$U, N1, u.dim)
    U2.start <- matrix(runif(N2 * u.dim, min = -1), ncol = u.dim)
    V2.tilde.start <- V2 %*% B.start
  } else if (init == 'user') {
    U.start <- matrix(matrix(U.start, N, u.dim)[V.sorted$ix,], N, u.dim)
    B.start <- matrix(B.start, q, q)
    V.tilde.start <- matrix(V.tilde.start, N, q)
    U1.start <- matrix(U.start[1:N1,], N1, u.dim)
    U2.start <- matrix(U.start[(N1+1):N,], N2, u.dim)
    V2.tilde.start <- matrix(V.tilde.start[V.sorted$ix[(N1+1):N],], N2, q)
  }

  B <- B.start
  U1 <- U1.start
  U2 <- U2.start
  U <- rbind(U1, U2)
  V1.tilde <- V1 %*% B
  V2.tilde <- V2.tilde.start
  V.tilde <- rbind(V1.tilde, V2.tilde)

  V.tilde.start <- matrix(V.tilde[V.sorted.rec.id,],,q)

  eta <- sum(w*d^2)
  gamma2 <- gamma*eta
  one_n_t <- t(rep(1,N))
  sigma <- rep(0, it.max)

  if (general) { # general weights ----
    H <- diag(rowSums(W)) - W
    dz <- condDist(U, V.tilde, one_n_t)
    sigma[1] <- sum(w*(d-dz)^2)
    Hp <- solve(H + 1) - N^(-2)
    H12H22p <- H[1:N1,(N1+1):N] %*% solve(H[(N1+1):N,(N1+1):N])
    H21V1 <- H[(N1+1):N,1:N1] %*% V1
    S <- V1t %*% H[1:N1,1:N1] %*% V1
    G <- H21V1 %*% solve(S)
    GV1t <- G %*% V1t
    Kvp <- solve(H[(N1+1):N,(N1+1):N] - G %*% t(H21V1))
    Kb <- S - V1t %*% H12H22p %*% H21V1
    Kbp <- solve(Kb)

    if (method == 'matrix') { #>> matrix ----
      KbpV1t <- Kbp %*% V1t

      for (iter in 2:it.max) {
        cz1 <- cz(w, d, dz)
        U <- Hp %*% (cz1 %*% U)
        B <- KbpV1t %*% ((cz1[1:N1,1:N1] - H12H22p %*% cz1[(N1+1):N,1:N1]) %*% V1.tilde +
                           (cz1[1:N1,(N1+1):N] - H12H22p %*% cz1[(N1+1):N,(N1+1):N]) %*% V2.tilde)
        V2.tilde <- Kvp %*% ((cz1[(N1+1):N,1:N1] - GV1t %*% cz1[1:N1,1:N1]) %*% V1.tilde +
                               (cz1[(N1+1):N,(N1+1):N] - GV1t %*% cz1[1:N1,(N1+1):N]) %*% V2.tilde)
        V1.tilde <- V1 %*% B
        V.tilde <- rbind(V1.tilde, V2.tilde)
        dz <- condDist(U, V.tilde, one_n_t)
        sigma[iter] <- sum(w*(d-dz)^2)
        if (sigma[iter - 1] - sigma[iter] < gamma2)
          break()
      }
    } else if (method == 'vector') { #>> vector ----
      b <- diag(B)

      for (iter in 2:it.max) {
        cz1 <- cz(w, d, dz)
        U <- Hp %*% (cz1 %*% U)
        for (m in 1:q)
          b[m] <- (V1t %*% ((cz1[1:N1,1:N1] - H12H22p %*% cz1[(N1+1):N,1:N1]) %*% V1[,m] * b[m] +
                  (cz1[1:N1,(N1+1):N] - H12H22p %*% cz1[(N1+1):N,(N1+1):N]) %*% V2.tilde[,m]))/Kb[m,m]
        V2.tilde <- Kvp %*% ((cz1[(N1+1):N,1:N1] - GV1t %*% cz1[1:N1,1:N1]) %*% V1.tilde +
                               (cz1[(N1+1):N,(N1+1):N] - GV1t %*% cz1[1:N1,(N1+1):N]) %*% V2.tilde)
        V1.tilde <- sapply(1:q, function(i) V1[,i]*b[i])
        V.tilde <- rbind(V1.tilde, V2.tilde)
        dz <- condDist(U, V.tilde, one_n_t)
        sigma[iter] <- sum(w*(d-dz)^2)
        if (sigma[iter - 1] - sigma[iter] < gamma2)
          break()
      }
      B <- diag(b, q, q)
    }
  } else { # weights of 1 ----
    dz <- condDist(U, V.tilde, one_n_t)
    sigma[1] <- sum(w*(d-dz)^2)
    V1i <- V1t %*% V1
    V1ip <- solve(V1i)
    v1s <- colSums(V1)
    v1so <- v1s %*% t(v1s)
    V1ipv1so <-  V1ip %*% v1so
    trV1ipv1so <- sum(diag(V1ipv1so))
    Sp <- (diag(1, q, q) + V1ipv1so/(N - trV1ipv1so)) %*% V1ip/N
    v1stSp <- t(v1s) %*% Sp
    v1stSpV1t <- v1stSp %*% V1t
    g <- c(1 + v1stSp %*% v1s)
    g.coeff <- g/(N - g*N2)
    Kbp <- (diag(1, q, q) + V1ipv1so/(N1 - trV1ipv1so)) %*% V1ip / N

    if (method == 'matrix') { #>> matrix ----
      KbpV1t <- Kbp %*% V1t

      if (exact) {
        for (iter in 2:it.max) {
          cz1 <- cz(w, d, dz)
          U <- scale(cz1 %*% U / N, center = TRUE, scale = FALSE)
          B <- KbpV1t %*% ((cz1[1:N1,1:N1] +
                              matrix(colSums(matrix(cz1[(N1+1):N,1:N1],N-N1,N1))/N1,
                                     N1, N1, byrow=TRUE)) %*% V1.tilde +
                             (cz1[1:N1,(N1+1):N] +
                                matrix(colSums(matrix(cz1[(N1+1):N,(N1+1):N],N-N1,N-N1))/N1,
                                       N1, N2, byrow=TRUE)) %*% V2.tilde)
          V2.tilde <- (cz1[(N1+1):N,1:N1] +
                         matrix(v1stSpV1t %*% cz1[1:N1,1:N1],
                                N2, N1, byrow=TRUE)) %*% V1.tilde +
            (cz1[(N1+1):N,(N1+1):N] +
               matrix(v1stSpV1t %*% cz1[1:N1,(N1+1):N],
                      N2, N2, byrow=TRUE)) %*% V2.tilde
          V2.tilde <- (V2.tilde + matrix(g.coeff*colSums(V2.tilde),
                                         N2, q, byrow=TRUE))/N
          V1.tilde <- V1 %*% B
          V.tilde <- rbind(V1.tilde, V2.tilde)
          dz <- condDist(U, V.tilde, one_n_t)
          sigma[iter] <- sum(w*(d-dz)^2)
          if (sigma[iter - 1] - sigma[iter] < gamma2)
            break()
        }
      } else {
        for (iter in 2:it.max) {
          cz1 <- cz(w, d, dz)
          U <- cz1 %*% U / N
          B <- KbpV1t %*% ((cz1[1:N1,1:N1] +
                              matrix(colSums(matrix(cz1[(N1+1):N,1:N1],N-N1,N1))/N1,
                                     N1, N1, byrow=TRUE)) %*% V1.tilde +
                             (cz1[1:N1,(N1+1):N] +
                                matrix(colSums(matrix(cz1[(N1+1):N,(N1+1):N],N-N1,N-N1))/N1,
                                       N1, N2, byrow=TRUE)) %*% V2.tilde)
          V2.tilde <- (cz1[(N1+1):N,1:N1] +
                         matrix(v1stSpV1t %*% cz1[1:N1,1:N1],
                                N2, N1, byrow=TRUE)) %*% V1.tilde +
            (cz1[(N1+1):N,(N1+1):N] +
               matrix(v1stSpV1t %*% cz1[1:N1,(N1+1):N],
                      N2, N2, byrow=TRUE)) %*% V2.tilde
          V2.tilde <- (V2.tilde + matrix(g.coeff*colSums(V2.tilde),
                                         N2, q, byrow=TRUE))/N
          V1.tilde <- V1 %*% B
          V.tilde <- rbind(V1.tilde, V2.tilde)
          dz <- condDist(U, V.tilde, one_n_t)
          sigma[iter] <- sum(w*(d-dz)^2)
          if (sigma[iter - 1] - sigma[iter] < gamma2)
            break()
        }
      }
    } else if (method == 'vector') { #>> vector ----
      b <- diag(B)
      Kb <- N*V1i - N/N1*v1so

      if (exact) {
        for (iter in 2:it.max) {
          cz1 <- cz(w, d, dz)
          U <- scale(cz1 %*% U / N, center = TRUE, scale = FALSE)
          for (m in 1:q)
            b[m] <- (V1t %*% ((cz1[1:N1,1:N1] + matrix(colSums(matrix(cz1[(N1+1):N,1:N1],N-N1,N1))/N1,
                                                       N1, N1, byrow=TRUE)) %*% V1[,m] * b[m] +
                                (cz1[1:N1,(N1+1):N] + matrix(colSums(matrix(cz1[(N1+1):N,(N1+1):N],N-N1,N-N1))/N1,
                                                             N1, N2, byrow=TRUE)) %*% V2.tilde[,m]))/Kb[m,m]
          V2.tilde <- (cz1[(N1+1):N,1:N1] +
                         matrix(v1stSpV1t %*% cz1[1:N1,1:N1],
                                N2, N1, byrow=TRUE)) %*% V1.tilde +
            (cz1[(N1+1):N,(N1+1):N] +
               matrix(v1stSpV1t %*% cz1[1:N1,(N1+1):N],
                      N2, N2, byrow=TRUE)) %*% V2.tilde
          V2.tilde <- (V2.tilde + matrix(g.coeff*colSums(V2.tilde),
                                         N2, q, byrow=TRUE))/N
          V1.tilde <- sapply(1:q, function(i) V1[,i]*b[i])
          V.tilde <- rbind(V1.tilde, V2.tilde)
          dz <- condDist(U, V.tilde, one_n_t)
          sigma[iter] <- sum(w*(d-dz)^2)
          if (sigma[iter - 1] - sigma[iter] < gamma2)
            break()
        }
      } else {
        for (iter in 2:it.max) {
          cz1 <- cz(w, d, dz)
          U <- cz1 %*% U / N
          for (m in 1:q)
            b[m] <- (V1t %*% ((cz1[1:N1,1:N1] + matrix(colSums(matrix(cz1[(N1+1):N,1:N1],N-N1,N1))/N1,
                                                       N1, N1, byrow=TRUE)) %*% V1[,m] * b[m] +
                                (cz1[1:N1,(N1+1):N] + matrix(colSums(matrix(cz1[(N1+1):N,(N1+1):N],N-N1,N-N1))/N1,
                                                             N1, N2, byrow=TRUE)) %*% V2.tilde[,m]))/Kb[m,m]
          V2.tilde <- (cz1[(N1+1):N,1:N1] +
                         matrix(v1stSpV1t %*% cz1[1:N1,1:N1],
                                N2, N1, byrow=TRUE)) %*% V1.tilde +
            (cz1[(N1+1):N,(N1+1):N] +
               matrix(v1stSpV1t %*% cz1[1:N1,(N1+1):N],
                      N2, N2, byrow=TRUE)) %*% V2.tilde
          V2.tilde <- (V2.tilde + matrix(g.coeff*colSums(V2.tilde),
                                         N2, q, byrow=TRUE))/N
          V1.tilde <- sapply(1:q, function(i) V1[,i]*b[i])
          V.tilde <- rbind(V1.tilde, V2.tilde)
          dz <- condDist(U, V.tilde, one_n_t)
          sigma[iter] <- sum(w*(d-dz)^2)
          if (sigma[iter - 1] - sigma[iter] < gamma2)
            break()
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

  tmp <- M <- matrix(0, N2, q)
  M[missing.ids] <- 1
  tmp[missing.ids] <- (V2.tilde - V2 %*% B)[missing.ids]
  V.hat <- rbind(V1, V2 + (tmp %*% solve(B))*M)[V.sorted.rec.id,]

  list(U = matrix(U[V.sorted.rec.id,],N,u.dim), B = B,
       stress = sigma[iter]/eta, sigma = sigma, init = init,
       U.start = matrix(rbind(U1.start, U2.start)[V.sorted.rec.id,]),
       B.start = B.start, V.tilde.start = V.tilde.start,
       V.hat = V.hat, method = method, exact = exact)
}

