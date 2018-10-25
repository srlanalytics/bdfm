BDFM <- function(Y, m, p, FC, Bp = NULL, lam_B = 0, Hp = NULL, lam_H = 0, nu_q = 0, nu_r = NULL, ID = "PC_full", ITC = T, reps = 1000, burn = 500, loud = FALSE) {

  # ----------- Preliminaries -----------------
  Y <- as.matrix(Y)
  k <- ncol(Y)
  r <- nrow(Y)
  n_obs <- sum(is.finite(Y))
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

  Parms <- EstDFM(B = B_in, Bp = Bp, lam_B = lam_B, q = q, nu_q = nu_q, H = H, Hp = Hp, lam_H = lam_H, R = Rvec, nu_r = nu_r, Y = Y, reps = reps, burn = burn, Loud = loud)

  if (ID %in% c("PC_sub", "PC_full")) {
    B <- Parms$B
    q <- Parms$Q
    H <- as.matrix(Parms$H[-(1:m), ])
    R <- diag(c(Parms$R[-(1:m)]))

    Est <- DSmooth(B, q, H, R, Y[, -(1:m), drop = FALSE])

    BIC <- log(n_obs)*(m*p + m^2 + k*m + k) - 2*Est$Lik

    Out <- list(
      B = B,
      q = q,
      H = H,
      R = R,
      itc = itc,
      values  = Est$Ys + matrix(1, r, 1) %x% t(itc),
      factors = Est$Z[,1:m],
      Qstore  = Parms$Qstore,
      Bstore  = Parms$Bstore,
      Rstore  = Parms$Rstore[-(1:m), ],
      Hstore  = Parms$Hstore[-(1:m), , ],
      Kstore  = Est$Kstr,
      PEstore = Est$PEstore,
      Lik     = Est$Lik,
      BIC     = BIC
    )
  } else {
    B <- Parms$B
    q <- Parms$Q
    H <- Parms$H
    R <- diag(c(Parms$R), k, k)

    Est <- DSmooth(B, q, H, R, Y)

    BIC <- log(n_obs)*(m*p + m^2 + k*m + k) - 2*Est$Lik

    Out <- list(
      B = B,
      q = q,
      H = H,
      R = R,
      values  = Est$Ys + matrix(1, r, 1) %x% t(itc),
      factors = Est$Z[,1:m],
      Qstore  = Parms$Qstore, # lets us look at full distribution
      Bstore  = Parms$Bstore,
      Rstore  = Parms$Rstore,
      Hstore  = Parms$Hstore,
      Kstore  = Est$Kstr,
      PEstore = Est$PEstore,
      Lik     = Est$Lik,
      BIC     = BIC
    )
  }
  return(Out)
}


#' MCMC Routine for Bayesian Dynamic Factor Models
#'
#' \code{CppBDFM} is the core C++ function for estimating a linear-Gaussian Bayesain dynamic factor model by MCMC methods using Durbin and Koopman's disturbance smoother. This function may be called directly by advanced users. The only dependencies are the  Armadillo (\url{http://arma.sourceforge.net/}) linear algebra library for C++ and the packages needed for interfacing with R (\code{\link{Rcpp}} and \code{\link{RcppArmadillo}}).
#'
#' @param B  initial guess for B in transition equation
#' @param Bp prior for B
#' @param lam_B prior tightness for B (additive)
#' @param q initial guess for q in the transition equation
#' @param nu_q prior "degrees of freedom" for inverse-Whishart prior for q (additive, prior scale is fixed so that increasing nu_q shrinks the variance towards zero)
#' @param H initial guess for H in the trasition equation
#' @param Hp prior for H
#' @param lam_H prior tightness for H (additive)
#' @param R initial guess for diagonal elements of R in the transition equation, entered as a vector
#' @param nu_r prior deg. of freedom for elements of R, entered as a vector (additive, prior scale is fixed so that increasing nu_r[j] shrinks the variance of shocks to series j towards zero)
#' @param Y Input data. Data must be scaled and centered prior to estimation if desired.
#' @param reps number of repetitions for MCMC sampling
#' @param burn number of iterations to burn in MCMC sampling
#' @param Loud print status of function during evaluation.
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib BDFM
CppBDFM <- function(B, Bp, lam_B, q, nu_q, H, Hp, lam_H, R, nu_r, Y, reps, burn, Loud){
  OUT <- EstDFM(B = B, Bp = Bp, lam_B = lam_B, q = q, nu_q = nu_q, H = H, Hp = Hp, lam_H = lam_H, R = R, nu_r = nu_r, Y = Y, reps = reps, burn = burn, Loud = Loud)
  return(OUT)
}