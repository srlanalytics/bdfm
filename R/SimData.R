#' Simulate Dynamic Factor Model
#'
#' @param r number of rows, i.e. time periods (T is TURE in R, but I use T in
#'   C++)
#' @export
#' @import MASS
SimDFM <- function(r) {
  m <- 3 # number of factors (built into function)
  k <- 10 # number of observed series
  burn <- 200
  tt <- r + burn

  # Transition matrix
  B <- matrix(c(
    .4, 0, 0, .3, 0, 0, .1, 0, 0,
    0, .3, .2, 0, .2, .1, 0, .1, 0,
    0, .2, .4, 0, 0, .2, 0, .1, .1
  ), m, m * 3, byrow = TRUE)

  # Use this to check whether B is stationary
  # A <- rbind(B,cbind(diag(1,6,6),matrix(0,6,3)))
  # tmp <- eigen(A)
  # tmp$values

  # Covariance matrix --- shocks to factors
  sigE <- matrix(c(
    1, 0, 0,
    0, 1, -.5,
    0, -.5, 1
  ), m, m, byrow = TRUE) # Changed from first submission

  V <- mvrnorm(tt, rep(0, m), sigE)

  X <- matrix(0, tt, m)
  for (ti in 4:tt) {
    X[ti, ] <- t(B %*% c(X[ti - 1, ], X[ti - 2, ], X[ti - 3, ])) + V[ti, ]
  }

  H <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, -1, 1, -1, .5, .5, .5, .5, .5, 1, 1, 0, -1, 1, 0, .5, -1, 0, .5, 1), 10, 3, byrow = T)

  R <- diag(1, k, k)
  U <- mvrnorm(r, rep(0, k), R)
  Y_true <- X[(burn + 1):tt, ] %*% t(H)
  Y_obs <- X[(burn + 1):tt, ] %*% t(H) + U
  ind <- rep(c(F, F, F, F, F, T, T, T, T, T), ceiling(r / 10))[1:r]
  Y_obs[ind, 1:3] <- NA

  return(list(
    Y_true = Y_true,
    Y_obs = Y_obs,
    X = X[(burn + 1):tt, ],
    B = B,
    q = sigE,
    H = H,
    r = diag(1, k, k)
  ))
}
