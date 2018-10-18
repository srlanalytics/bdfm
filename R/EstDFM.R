
# ----- Helper Functions ---------

AnyNA <- function(y) {
  any(is.na(y))
}

NumNA <- function(y) {
  length(which(is.na(y)))
}

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

# ----------------------------------------

#' Estimate Bayesian dynamic factor model
#'
#' @param Y Data in matrix format with time in rows
#' @param factors number of factors 
#' @param lags number of lags in transition equation
#' @param forecast number of periods ahead to forecast
#' @param B_prior prior matrix for B in the transition equation. Default is zeros.
#' @param lam_B prior tightness on B 
#' @param H_prior prior matrix for H (loadings) in the observation equation. Default is zeros. 
#' @param lam_H prior tightness on H
#' @param nu_q prior deg. of freedom for transition equation, entered as vector with length equal to the number of factors.  
#' @param nu_r prior deg. of freedom for observables, entered as vector with length equal to the number of observables.
#' @param ID Factor Identification. 'PC_full' is the default (using all observed series), 'PC_sub' finds a submatrix of the data that maximizes the number of observations for a square (no missing values) data set. Use 'PC_sub' when many observations are missing. 
#' @param reps number of repetitions for MCMC sampling
#' @param burn number of iterations to burn in MCMC sampling
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib BDFM
bdfm <- function(Y, factors = 1, lags = 2, forecast = 0, B_prior = NULL, lam_B = 0, H_prior = NULL, lam_H = 0, nu_q = 0, nu_r = NULL, ID = 'PC_full', intercept = TRUE, reps = 1000, burn = 500) {
  
  try_tsbox <- try(library(tsbox), silent = T) #Is tsbox present?
  if(inherits(try_tsbox,"try-error")){ #If not, output warning
    warning("Preserving time series attributes of input data requires the packages tsbox; time series attributes of inputs will be lost.")
    Y <- as.matrix(Y)
    E <- BDFM(Y = Y, m = factors, p = lags, FC = forecast, Bp = B_prior, lam_B = lam_B, Hp = H_prior, lam_H = lam_H, nu_q = nu_q, nu_r = nu_r, ID = ID, ITC = intercept, reps = reps, burn = burn)
    colnames(E$values) <- colnames(Y)
  }else{
    # non time series
    if (!ts_boxable(Y)) {
      if (is.matrix(Y)) {
        E <- BDFM(Y = Y, m = factors, p = lags, FC = forecast, Bp = B_prior, lam_B = lam_B, Hp = H_prior, lam_H = lam_H, nu_q = nu_q, nu_r = nu_r, ID = ID, ITC = intercept, reps = reps, burn = burn)
        colnames(E$values) <- colnames(Y)
      }else{
        stop('Data format not accepted. Data must be ts_boxable() or matrix format.')
      }
    }else{  # time series --- i.e. is ts_boxable

    # convert to mts
    Y.uc  <- unclass(ts_ts(Y))
    Y.tsp <- attr(Y.uc, "tsp")
    attr(Y.uc, "tsp") <- NULL
    
    E <- BDFM(Y = Y.uc, m = factors, p = lags, FC = forecast, Bp = B_prior, lam_B = lam_B, Hp = H_prior, lam_H = lam_H, nu_q = nu_q, nu_r = nu_r, ID = ID, ITC = intercept, reps = reps, burn = burn)
    
    ts(E$values, start = Y.tsp[1], frequency = Y.tsp[3])
    ts(E$factors, start = Y.tsp[1], frequency = Y.tsp[3])
    colnames(E$values) <- colnames(Y)
    
    E$values  <- copy_class(E$values, Y)
    E$factors <- copy_class(E$factors, Y)
    
    
    }
  }
  class(E) <- "bdfm"
  return(E)
}
  

BDFM <- function(Y, m, p, FC, Bp = NULL, lam_B = 0, Hp = NULL, lam_H = 0, nu_q = 0, nu_r = NULL, ID = "PC_full", ITC = T, reps = 1000, burn = 500) {
  
  # ----------- Preliminaries -----------------
  Y <- as.matrix(Y)
  k <- ncol(Y)
  r <- nrow(Y)
  if (ITC) {
    itc <- colMeans(Y, na.rm = T)
    Y <- Y - matrix(1, r, 1) %x% t(itc) # De-mean data. Data is not automatically stadardized --- that is left up to the user and should be done before estimation if desired.
  } else {
    itc <- rep(0, k)
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
      stop("Every period contains missing data. Try setting ID to PC_sub.")
    }
  }
  
  # ----------- Format Priors ------------------
  #enter priors multiplicatively so that 0 is a weak prior and 1 is a strong prior (additive        priors are relative to the number of observations)
  lam_B <- r*lam_B + 1
  nu_q  <- r*nu_q  + lam_B
  lam_H <- r*lam_H + 1 
  if(is.null(nu_r)){
    nu_r = rep(1,k)*lam_H
  }else{
    if(length(nu_r) != k){stop("Length of nu_r must equal the number of observed series")}
    nu_r = rep(1,k)*lam_H + r*nu_r
  }
  # ---------------------------------------------
  

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

  Parms <- EstDFM(B = B_in, Bp = Bp, lam_B = lam_B, q = q, nu_q = nu_q, H = H, Hp = Hp, lam_H = lam_H, R = Rvec, nu_r = nu_r, Y = Y, reps = reps, burn = burn)

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
      values  = Est$Ys,
      factors = Est$Z[,1:m],
      Qstore  = Parms$Qstore,
      Bstore  = Parms$Bstore,
      Rstore  = Parms$Rstore[-(1:m), ],
      Hstore  = Parms$Hstore[-(1:m), , ],
      Kstore  = Est$Kstr,
      PEstore = Est$PEstore
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
      values  = Est$Ys,
      factors = Est$Z[,1:m],
      Qstore  = Parms$Qstore, # lets us look at full distribution
      Bstore  = Parms$Bstore,
      Rstore  = Parms$Rstore,
      Hstore  = Parms$Hstore,
      Kstore  = Est$Kstr,
      PEstore = Est$PEstore
    )
  }
  return(Out)
}

#' ML Estimation for Dynamic Factor Models
#' 
#' Estimate a dynamic factor model by maximum likelihood using the EM Algorithm in Watson and Engle (1983). 
#' 
#' @param Y data with time indexed by row
#' @param factors number of factors to estimate
#' @param lags number of lags in the transition equation
#' @param forecast number of periods ahead to forecast
#' @param tol tollerance for convergence of the EM Algorithm 
#' @param Loud T/F, print change in log likelihood at each iteration
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib BDFM
mlDFM <- function(Y, factors = 1, lags = 2, forecast = 0, tol = 0.01, Loud = FALSE) {
  
  try_tsbox <- try(library(tsbox), silent = T) #Is tsbox present?
  if(inherits(try_tsbox,"try-error")){ #If not, output warning
    warning("Preserving time series attributes of input data requires the packages tsbox; time series attributes of inputs will be lost.")
    Y <- as.matrix(Y)
    E <- MLdfm(Y, m = factors, p = lags, FC = forecast, tol = tol, Loud = Loud)
    colnames(E$values) <- colnames(Y)
  }else{
    # non time series
    if (!ts_boxable(Y)) {
      if (is.matrix(Y)) {
        E <- MLdfm(Y, m = factors, p = lags, FC = forecast, tol = tol, Loud = Loud)
        colnames(E$values) <- colnames(Y)
      }else{
        stop('Data format not accepted. Data must be ts_boxable() or matrix format.')
      }
    }else{  # time series --- i.e. is ts_boxable
      
      # convert to mts
      Y.uc  <- unclass(ts_ts(Y))
      Y.tsp <- attr(Y.uc, "tsp")
      attr(Y.uc, "tsp") <- NULL
      
      E <- MLdfm(Y, m = factors, p = lags, FC = forecast, tol = tol, Loud = Loud)
      
      ts(E$values, start = Y.tsp[1], frequency = Y.tsp[3])
      ts(E$factors, start = Y.tsp[1], frequency = Y.tsp[3])
      colnames(E$values) <- colnames(Y)
      
      E$values   <- copy_class(E$values, Y)
      E$factors <- copy_class(E$factors, Y)
    }
  }
  class(E) <- "mlDFM"
  return(E)
}

MLdfm <- function(Y, m, p, FC = 0, tol = 0.01, Loud = FALSE) {

  Y <- as.matrix(Y)
  r <- nrow(Y)
  k <- ncol(Y)
  itc <- colMeans(Y, na.rm = TRUE)
  Ytmp <- na.omit(Y)
  sA <- m * (p + 1) # number of factors and size of A matrix

  Ytmp <- scale(Ytmp) # scale before taking principle components
  scl  <- attr(Ytmp, "scaled:scale")
  # loadings on principle components and initial guess for H
  PC   <- PrinComp(Ytmp, m)
  H    <- PC$loadings

  if (sum(sign(H[, 1])) < 0) {
    H[, 1] <- -H[, 1]
  }
  H <- H * (matrix(1, 1, m) %x% scl) # scale H up

  # Arbitrary initial guess for A
  A <- Matrix(0, sA, sA)
  A[1:m, 1:m] <- .1 * Diagonal(m)
  A[(m + 1):sA, 1:(m * p)] <- Diagonal(m * p)

  # Arbirary initial guess for Q
  Q <- Matrix(0, sA, sA)
  Q[1:m, 1:m] <- Diagonal(m)

  # Arbitrary intitial guess for R
  R <- diag(1, k, k)
  
  if(FC>0){
    tmp <- matrix(NA,FC,k)
    Y   <- rbind(Y, tmp)
    r   <- r + FC
  }

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
  B    <- A[1:m,1:(m*p)]

  return(list(
    values = Smth$Ys + matrix(1, r, 1) %x% t(itc),
    Lik = Smth$Lik,
    factors = Smth$Z[,1:m],
    B = B,
    Q = Q,
    H = H,
    R = R,
    A = A,
    HJ = HJ,
    itc = itc,
    Kstore  = Smth$Kstr,
    PEstore = Smth$PEstr
  ))
}
#' 
#' # methods
#' #' @export
#' #' @method predict bdfm
#' predict.bdfm <- function(x) {
#'   x$values
#' }
#' 
#' #' @export
#' #' @method print bdfm
#' print.bdfm <- function(x) {
#'   print("a bdfm model: short description")
#' }
#' 
#' #' @export
#' #' @method summary bdfm
#' summary.bdfm <- function(x) {
#'   print("a bdfm model: summary description")
#' }
