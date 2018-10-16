
# Helper Functions

#' @export
AnyNA <- function(y) {
  any(is.na(y))
}

#' @export
NumNA <- function(y) {
  length(which(is.na(y)))
}

#' @export
Y_sub <- function(Y) {
  # lots of ways to select a full submatrix of the data... this is one idea
  Ysub <- na.omit(Y)
  obs_old <- nrow(Ysub) * ncol(Ysub)
  miss_rows <- apply(Y, MARGIN = 2, FUN = NumNA)
  dims <- dim(na.omit(Y[, -which(miss_rows == max(miss_rows))]))
  obs_new <- dims[1] * dims[2]
  if (obs_new < obs_old) {
    ind <- apply(Y, MARGIN = 1, FUN = AnyNA)
    Ysub <- Y[!ind, ]
  }
  while(obs_new > obs_old || obs_old == 0){
    obs_old <- obs_new
    miss_rows <- apply(Y, MARGIN = 2, FUN = NumNA)
    Y <- Y[, -which(miss_rows == max(miss_rows))]
    ind <- apply(Y, MARGIN = 1, FUN = AnyNA) # indexes of Y that do not correspond to Ysub
    Ysub <- Y[!ind, ]
    # Ysub      <- na.omit(Y)
    obs_new <- nrow(Ysub) * ncol(Ysub)
    # print(dim(Ysub))
  }
  return(list(Ysub = Ysub, ind = !ind))
}


#' Estimate uniform frequency dynamic factor model
#'
#' @param Y Data in matrix format with time in rows
#' @param m number of factors m
#' @param p number of lags in transition equation
#' @param FC number of periods ahead to forecast
#' @param lam_B prior tightness on B (Bp, the prior for B in the tranistion eq. is built in as zero)
#' @param lam_H prior tightness on H (additive. A value equal to the number of time periods in Data is a pretty strong prior.)
#' @param nu_q prior deg. of freedom for transition equation
#' @param nu_r prior deg. of freedom for observables used to identify the model
#' @param ID Factor Identification. The default is to use principal components.
#' @param reps number of repetitions from which to sample
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib BDFM
BDFM <- function(Y, m, p, FC = 0, Bp = NULL, lam_B = 1, Hp = NULL, lam_H = 1, nu_q = NULL, nu_r = NULL, ID = "PC_full", ITC = T, reps = 1000) {

  # ----------- Preliminaries -----------------

  k <- ncol(Y)
  r <- nrow(Y)
  if (ITC) {
    itc <- colMeans(Y, na.rm = T)
    Y <- Y - matrix(1, r, 1) %x% t(itc) # De-mean data. Data is not automatically stadardized --- that is left up to the user and should be done before estimation if desired.
  } else {
    itc <- rep(0, k)
  }
  if(is.null(nu_q)){
    nu_q = lam_B
  }
  if(is.null(nu_r)){
    nu_r = lam_H
  }

  # Identification is based on the first m variables so the initial guess just uses these variables as factors.

  if (ID == "PC_sub") {
    Ysub <- Y_sub(Y) # submatrix of Y with complete data, i.e. no missing values
    PC <- PrinComp(Ysub$Ysub, m)
    tmp <- matrix(NA, r, m)
    tmp[Ysub$ind, ] <- PC$components
    Y <- cbind(tmp, Y)
    k <- k + m
  } else if (ID == "PC_full") {
    PC <- PrinComp(Y, m)
    Y <- cbind(PC$components, Y)
    k <- k + m
    if (!any(!is.na(PC$components))) {
      stop("Every period contains missing data. Try a different identification method.")
    }
  }

  H <- matrix(0, k, m)
  H[1:m, 1:m] <- diag(1, m, m)
  # for variables not used to normalize
  if (k > m) {
    xx <- as.matrix(Y[, 1:m])
    yy <- as.matrix(Y[, -(1:m)])
    tmp <- t(QuickReg(xx, yy))
    H[-(1:m), ] <- tmp
  }

  Rvec <- rep(1, k) # arbitrary initial guess

  B_in <- matrix(0, m, m * p)
  B_in[1:m, 1:m] <- diag(.1, m, m) # arbitrary initial guess

  q <- diag(1, m, m)

  if (is.null(Bp)) {
    Bp <- matrix(0, m, m * p)
  }
  if (is.null(Hp)) {
    Hp <- matrix(0, k, m)
  }
  if(FC>0){
    tmp <- matrix(NA,FC,k)
    Y   <- rbind(Y, tmp)
    r   <- r + FC
  }

  Parms <- EstDFM(B = B_in, Bp = Bp, lam_B = lam_B, q = q, nu_q = nu_q, H = H, Hp = Hp, lam_H = lam_H, R = Rvec, nu_r = nu_r, Y = Y, reps = reps)

  if (ID %in% c("PC_sub", "PC_full")) {
    B <- Parms$B
    q <- Parms$Q
    H <- as.matrix(Parms$H[-(1:m), ])
    R <- diag(c(Parms$R[-(1:m)]), k, k)

    Est <- DSmooth(B, q, H, R, Y[, -(1:m)])

    Out <- list(
      B = B,
      q = q,
      H = H,
      R = R,
      itc = itc,
      Ys = Est$Ys + matrix(1, r, 1) %x% t(itc),
      Z = Est$Z,
      Qstore = Parms$Qstore,
      Bstore = Parms$Bstore,
      Rstore = Parms$Rstore[-(1:m), ],
      Hstore = Parms$Hstore[-(1:m), , ]
    )
  } else {
    B <- Parms$B
    q <- Parms$Q
    H <- Parms$H
    R <- diag(c(Parms$R), k, k)

    Est <- DSmooth(B, q, H, R, Y)

    Out <- list(
      B = B,
      q = q,
      H = H,
      R = R,
      Ys = Est$Ys + matrix(1, r, 1) %x% t(itc),
      Z = Est$Z,
      Qstore = Parms$Qstore, # lets us look at full distribution
      Bstore = Parms$Bstore,
      Rstore = Parms$Rstore,
      Hstore = Parms$Hstore
    )
  }
  return(Out)
}


#' @export
mlDFM <- function(Y, m, p, tol = 0.01, Loud = FALSE) {
  Y <- as.matrix(Y)
  r <- nrow(Y)
  k <- ncol(Y)
  itc <- colMeans(Y, na.rm = TRUE)
  tmp <- na.omit(Y)
  rtmp <- nrow(tmp)
  sA <- m * (p + 1) # number of factors and size of A matrix

  Ydm <- tmp - matrix(1, rtmp, 1) %x% t(itc)
  scle <- sqrt(colMeans(Ydm^2))
  Ytmp <- Ydm / (matrix(1, rtmp, 1) %x% t(scle)) # scale before taking principle components
  # loadings on principle components and initial guess for H
  PC <- PrinComp(Ytmp, m)
  H <- PC$loadings

  if (sum(sign(H[, 1])) < 0) {
    H[, 1] <- -H[, 1]
  }
  H <- H * (matrix(1, 1, m) %x% scle) # scale H up

  # Arbitrary initial guess for A
  A <- Matrix(0, sA, sA)
  A[1:m, 1:m] <- .1 * Diagonal(m)
  A[(m + 1):sA, 1:(m * p)] <- Diagonal(m * p)

  # Arbirary initial guess for Q
  Q <- Matrix(0, sA, sA)
  Q[1:m, 1:m] <- Diagonal(m)

  # Arbitrary intitial guess for R
  R <- diag(1, k, k)

  count <- 0
  Lik0 <- -100000
  Conv <- 100
  while (Conv > tol | count < 5) {
    Est <- KestExact(A, Q, H, R, Y, itc, m, p)
    A <- Est$A
    Q <- Est$Q
    H <- Est$H
    R <- Est$R
    itc <- Est$itc
    Lik1 <- Est$Lik
    Conv <- 200 * (Lik1 - Lik0) / abs(Lik1 + Lik0)
    Lik0 <- Lik1
    if (Loud) {
      print(Conv)
    }
    count <- count + 1
  }


  ############ Final Estimates ############
  if (m * p == 1) {
    A <- sparseMatrix(i = 1, j = 1, x = c(A[1, 1]), dims = c(1, 1), symmetric = FALSE, triangular = FALSE, giveCsparse = TRUE)
    Q <- sparseMatrix(i = 1, j = 1, x = c(Q[1, 1]), dims = c(1, 1), symmetric = FALSE, triangular = FALSE, giveCsparse = TRUE)
  }
  else {
    A <- A[1:(m * p), 1:(m * p)]
    Q <- Q[1:(m * p), 1:(m * p)]
  }
  Ydm <- Y - matrix(1, r, 1) %x% t(itc)
  HJ <- sparseMatrix(i = rep(1:k, m), j = (1:m) %x% rep(1, k), x = c(H), dims = c(k, m * p), symmetric = FALSE, triangular = FALSE, giveCsparse = TRUE)

  Smth <- Ksmoother(A, Q, HJ, R, Ydm)

  Ys <- Smth$Ys + matrix(1, r, 1) %x% t(itc)
  Yf <- Smth$Yf + matrix(1, r, 1) %x% t(itc)
  Z <- Smth$Z

  return(list(
    Ys = Ys,
    Yf = Yf,
    Lik = Smth$Lik,
    Z = Z,
    itc = itc,
    A = A,
    Q = Q,
    HJ = HJ,
    R = R,
    K = Smth$Kstr,
    PE = Smth$PEstr
  ))
}
