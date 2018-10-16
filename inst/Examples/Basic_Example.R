library(tsbox)

# this is to produce a pdf on the server, rather than plot it. Temporary
# soluthion for Christoph's setting.

gr <- function(...) {
  pdf("plots.pdf", width = 10, height = 7)
  eval(...)
  dev.off()
}

# --- Simulate Data from a state space model ---
# In this example there are 2 states, 4 observables, and 2 lags in the transition equation

# States
library(MASS)
r <- 200 # number of observables (rows)
A <- matrix(c(1, -.5, .1, .7), 2, 2, byrow = T) # transition matrix
Q <- diag(1, 2, 2) # variance of shocks to transition equation
X <- matrix(0, r, 2) # states
E <- matrix(rnorm(2 * r), r, 2) # shocks to transition equation
for (j in 2:r) {
  X[j, ] <- X[j - 1, ] %*% t(A) + E[j, ]
}
# gr(ts.plot(X, col = 1:2))

# Observations
H <- matrix(c(.5, 1, -1, 2, 1, -1, 1, -.5), 4, 2, byrow = T)
R <- diag(1, 4, 4)
Eps <- matrix(rnorm(4 * r), r, 4)
Y <- X %*% t(H) + Eps

# gr(ts.plot(Y))

# ts.plot(Y)

# --- Filtering and Smoothing ---
# One could use steady state filter values here, but in practice we will begin with difuse states

# Initializing
m <- 2 # number of factors in Z
Z <- matrix(0, r, m) # initial states are zero
Zp <- matrix(0, r, m) # store predictions for smoother
p0 <- 10^5 * diag(1, m, m) # large (difuse) initial variance
P0 <- array(0, c(m, m, r)) # stored values for smoother
P0[, , 1] <- p0
P1 <- array(0, c(m, m, r)) # stored values for smoother
P1[, , 1] <- p0

# Filtering
for (j in 2:r) {
  Zp[j, ] <- A %*% Z[j - 1, ] # prediction step
  p1 <- A %*% p0 %*% t(A) + Q # variance of Zp
  p1 <- (p1 + t(p1)) / 2 # reduce rounding error
  P1[, , j] <- p1 # store value for smoother
  yp <- H %*% Zp[j, ] # predictions for observables
  S <- H %*% p1 %*% t(H) + R # variance of yp
  C <- p1 %*% t(H) # covariance of Zp and yp
  Z[j, ] <- t(Zp[j, ]) + t(Y[j, ] - yp) %*% solve(S, t(C)) # update factors
  p0 <- p1 - C %*% solve(S, t(C)) # factor variance
  p0 <- (p0 + t(p0)) / 2 # reduce rounding error
  P0[, , j] <- p0 # store value for smoother
}

gr(ts.plot(cbind(Z[, 1], X[, 1]), col = c("steelblue", "red"))) # plot filtered values

# Smoothing
Zs <- Z
Ps <- P0
for (j in (r - 1):1) {
  g <- P0[, , j] %*% t(solve(P1[, , j + 1], A))
  Zs[j, ] <- Z[j, ] + (Zs[j + 1, ] - Zp[j + 1, ]) %*% t(g)
  ps <- P0[, , j] - g %*% (P1[, , j + 1] - Ps[, , j + 1]) %*% t(g)
  Ps[, , j] <- (ps + t(ps)) / 2
}

gr({
  ts.plot(cbind(Zs[, 1], X[, 1]), col = c("steelblue", "red"), lwd = 2, lty = c(1, 2)) # plot smoothed values
  legend("topleft",
    legend = c("Estimated Values", "True (Unobserved) Values"),
    lty = c(1, 2), lwd = 2, col = c("steelblue", "red")
  )
})
# Dealing with missing values

Ymv <- Y
# ind      <- matrix(rnorm(4*r),r,4) >1 #randomly assigning missing values
# Ymv[ind] <- NA
Ymv[30:100, 1] <- NA
Ymv[70:140, 2] <- NA

# --- Filtering and Smoothing ---
# One could use steady state filter values here, but in practice we will begin with difuse states

# Initializing
m <- 2 # number of factors in Z
Z <- matrix(0, r, m) # initial states are zero
Zp <- matrix(0, r, m) # store predictions for smoother
p0 <- 10^5 * diag(1, m, m) # large (difuse) initial variance
P0 <- array(0, c(m, m, r)) # stored values for smoother
P0[, , 1] <- p0
P1 <- array(0, c(m, m, r)) # stored values for smoother
P1[, , 1] <- p0

# Filtering
for (j in 2:r) {
  Zp[j, ] <- A %*% Z[j - 1, ] # prediction step
  p1 <- A %*% p0 %*% t(A) + Q # variance of Zp
  p1 <- (p1 + t(p1)) / 2 # reduce rounding error
  P1[, , j] <- p1 # store value for smoother
  ind <- is.finite(Ymv[j, ]) # identify non-missing values of Y[j,]
  yp <- H[ind, ] %*% Zp[j, ] # predictions for observables
  S <- H[ind, ] %*% p1 %*% t(H[ind, ]) + R[ind, ind] # variance of yp
  C <- p1 %*% t(H[ind, ]) # covariance of Zp and yp
  Z[j, ] <- t(Zp[j, ]) + t(Y[j, ind] - yp) %*% solve(S, t(C)) # update factors
  p0 <- p1 - C %*% solve(S, t(C)) # factor variance
  p0 <- (p0 + t(p0)) / 2 # reduce rounding error
  P0[, , j] <- p0 # store value for smoother
}

gr({
  ts.plot(cbind(Z[, 1], X[, 1]), col = c("steelblue", "red")) # plot filtered values
})

# Smoothing
Zs <- Z
Ps <- P0
for (j in (r - 1):1) {
  g <- P0[, , j] %*% t(solve(P1[, , j + 1], A))
  Zs[j, ] <- Z[j, ] + (Z[j + 1, ] - Zp[j + 1, ]) %*% t(g)
  ps <- P0[, , j] - g %*% (P1[, , j + 1] - Ps[, , j + 1]) %*% t(g)
  Ps[, , j] <- (ps + t(ps)) / 2
}

gr({
  ts.plot(cbind(Zs[, 1], X[, 1]), col = c("steelblue", "red")) # plot smoothed values
})

# Fill in missing values in Ymv
Y_hat <- X %*% t(H)
Y_fill <- Ymv
Y_fill[!is.finite(Ymv)] <- Y_hat[!is.finite(Ymv)]


gr({
  ts.plot(cbind(Y_fill[, 1], Ymv[, 1], Y[, 1]), col = c("steelblue", "red", "red"), lty = c(1, 1, 2), lwd = c(2, 2, 1)) # plot smoothed values
  legend("topleft",
    legend = c("Estimated Values", "True Observed Values", "True Missing Values"),
    lty = c(1, 1, 2), lwd = c(2, 2, 1), col = c("steelblue", "red", "red")
  )
})
