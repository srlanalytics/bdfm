bdfm <- function(Y, m, p, Bp, lam_B, Hp, lam_H, nu_q, nu_r, ID, store_idx, freq, LD, reps, burn, verbose, orthogonal_shocks) {

  # Preliminaries
  Y <- as.matrix(Y)
  k <- ncol(Y)
  r <- nrow(Y)
  n_obs <- sum(is.finite(Y))

  if (is.null(freq)) { # uniform frequency
    freq <- rep(1, k)
  }
  if (is.null(LD)) { # don't worry about aggregating differenced low freq. data
    LD <- rep(0, k)
  }

  # If data is mixed frequency number of lags needed may be bigger than p
  nlags <- rep(0, k)
  for (j in 1:k) {
    if (LD[j] == 0) {
      nlags[j] <- freq[j]
    } else if (LD[j] == 1) {
      nlags[j] <- 2 * freq[j] - 1
    } else {
      stop("Invalid diffs argument. Values must be 0 for level data or 1 for differenced data. Second differences not supported.")
    }
  }
  pp <- max(c(nlags, p))

  Jb <- Diagonal(m * p)
  if (pp > p) {
    Jb <- cbind(Jb, Matrix(0, m * p, m * (pp - p)))
  }

  # Identification is based on the first m variables so the initial guess just uses these variables as factors.

  if (is.numeric(ID)) {
    if(length(ID)<m){
      stop("Number of factors is too great for selected identification routine. Try fewer factors or 'pc_wide'")
    }
    PC <- PrinComp(Y[, ID], m)
    Y <- cbind(PC$components, Y)
    k <- k + m
    if (!is.null(nu_r)) {
      nu_r <- c(rep(0, m), nu_r)
    }
    freq <- c(rep(1, m), freq)
    LD <- c(rep(0, m), LD)
  } else if (ID == "pc_wide") {
    PC <- PrinComp(Y, m)
    Y <- cbind(PC$components, Y)
    k <- k + m
    if (!is.null(nu_r)) {
      nu_r <- c(rep(0, m), nu_r)
    }
    freq <- c(rep(1, m), freq)
    LD <- c(rep(0, m), LD)
    if (!any(!is.na(PC$components))) {
      stop("Every period contains missing data. Try setting ID to pc_long.")
    }
  }else if (ID == "pc_long") {
    long <- apply(Y,2,function(e) sum(is.finite(e)))
    long <- long>=median(long)
    if(sum(long)<m){
      stop("Number of factors is too great for selected identification routine. Try fewer factors or 'pc_wide'")
    }
    PC <- PrinComp(Y[,long, drop = FALSE], m)
    Y <- cbind(PC$components, Y)
    k <- k + m
    if (!is.null(nu_r)) {
      nu_r <- c(rep(0, m), nu_r)
    }
    freq <- c(rep(1, m), freq)
    LD <- c(rep(0, m), LD)
    if (!any(!is.na(PC$components))) {
      stop("Every period contains missing data. Try setting ID to pc_sub.")
    }
  } else if (ID != "name") {
    warning(paste(ID, "is not a valid identification string or index vector, defaulting to pc_long"))
    ID <- "pc_long"
    long <- apply(Y,2,function(e) sum(is.finite(e)))
    long <- long>=median(long)
    if(sum(long)<m){
      stop("Number of factors is too great for selected identification routine. Try fewer factors or 'pc_wide'")
    }
    PC <- PrinComp(Y[,long, drop = FALSE], m)
    Y <- cbind(PC$components, Y)
    k <- k + m
    if (!is.null(nu_r)) {
      nu_r <- c(rep(0, m), nu_r)
    }
    freq <- c(rep(1, m), freq)
    LD <- c(rep(0, m), LD)
    if (!any(!is.na(PC$components))) {
      stop("Every period contains missing data. Try setting ID to pc_sub.")
    }
  }

  # Format Priors
  # enter priors multiplicatively so that 0 is a weak prior and 1 is a strong
  # prior (additive priors are relative to the number of observations)
  lam_B <- r * lam_B + 1
  nu_q <- r * nu_q + 1
  lam_H <- r * lam_H + 1
  if (is.null(nu_r)) {
    nu_r <- rep(1, k)
  } else {
    if (length(nu_r) != k) {
      stop("Length of nu_r must equal the number of observed series")
    }
    nu_r <- r * nu_r + rep(1, k)
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
  if (is.null(store_idx)) {
    store_idx <- 0
    store_Y <- FALSE
  } else {
    store_Y <- TRUE
    if (ID %in% c("pc_wide", "pc_long") || is.numeric(ID)) {
      store_idx <- store_idx + m - 1 # -1 due to zero indexing in C++
    } else {
      store_idx <- store_idx - 1
    }
  }

  Parms <- EstDFM(B = B_in, Bp = Bp, Jb = Jb, lam_B = lam_B, q = q, nu_q = nu_q, H = H, Hp = Hp, lam_H = lam_H, R = Rvec, nu_r = nu_r, Y = Y, freq = freq, LD = LD, store_Y = store_Y, store_idx = store_idx, reps = reps, burn = burn, verbose = verbose)

  if (ID %in% c("pc_wide", "pc_long") || is.numeric(ID)) {
    B <- Parms$B
    q <- Parms$Q
    H <- as.matrix(Parms$H[-(1:m), ])
    R <- diag(c(Parms$R[-(1:m)]), nrow = k-m, ncol = k-m)

    Y <- Y[, -(1:m), drop = FALSE] #drop components used to identify model
    
    if(orthogonal_shocks){ #if we want to return a model with orthogonal shocks, rotate the parameters
      id <- Identify(H,q)
      H  <- H%*%id[[1]]
      B  <- id[[2]]%*%B%*%(diag(1,p,p)%x%id[[1]])
      q  <- id[[2]]%*%q%*%t(id[[2]])
    }

    Est <- DSmooth(
      B = B, Jb = Jb, q = q, H = H, R = R,
      Y = Y, freq = freq[-(1:m)], LD = LD[-(1:m)]
    )

    # stopifnot(!all(is.na(Est$Ys)))

    #Format output a bit
    rownames(H) <- colnames(Y)
    R <- diag(R)
    names(R) <- colnames(Y)

    BIC <- log(n_obs) * (m * p + m^2 + k * m + k) - 2 * Est$Lik

    Out <- list(
      B = B,
      q = q,
      H = H,
      R = R,
      Jb = Jb,
      values = Est$Ys, # + matrix(1, r, 1) %x% t(itc),
      factors = Est$Z[, 1:m],
      unsmoothed_factors = Est$Zz[, 1:m],
      predicted_factors  = Est$Zp[, 1:m],
      Qstore = Parms$Qstore,
      Bstore = Parms$Bstore,
      Hstore = Parms$Hstore[-(1:m), , , drop = FALSE],
      Rstore = Parms$Rstore[-(1:m), , drop = FALSE],
      Kstore = Est$Kstr,
      PEstore = Est$PEstr,
      Lik = Est$Lik,
      BIC = BIC,
      Ystore = Parms$Ystore,
      Ymedian = Parms$Y_median
    )
  } else {
    B <- Parms$B
    q <- Parms$Q
    H <- Parms$H
    R <- diag(c(Parms$R), k, k)
    
    if(orthogonal_shocks){ #if we want to return a model with orthogonal shocks, rotate the parameters
      id <- Identify(H,q)
      H  <- H%*%id[[1]]
      B  <- id[[2]]%*%B%*%(diag(1,p,p)%x%id[[1]])
      q  <- id[[2]]%*%q%*%t(id[[2]])
    }

    Est <- DSmooth(B = B, Jb = Jb, q = q, H = H, R = R, Y = Y, freq = freq, LD = LD)

    #Format output a bit
    rownames(H) <- colnames(Y)
    R <- diag(R)
    names(R) <- colnames(Y)

    BIC <- log(n_obs) * (m * p + m^2 + k * m + k) - 2 * Est$Lik

    Out <- list(
      B = B,
      q = q,
      H = H,
      R = R,
      Jb = Jb,
      values = Est$Ys, # + matrix(1, r, 1) %x% t(itc),
      factors = Est$Z[, 1:m],
      unsmoothed_factors = Est$Zz[, 1:m],
      predicted_factors  = Est$Zp[, 1:m],
      Qstore = Parms$Qstore, # lets us look at full distribution
      Bstore = Parms$Bstore,
      Hstore = Parms$Hstore,
      Rstore = Parms$Rstore,
      Kstore = Est$Kstr,
      PEstore = Est$PEstr,
      Lik = Est$Lik,
      BIC = BIC,
      Ystore = Parms$Ystore,
      Ymedian = Parms$Y_median
    )
  }
  return(Out)
}


#' MCMC Routine for Bayesian Dynamic Factor Models
#'
#' \code{Cppbdfm} is the core C++ function for estimating a linear-Gaussian
#' Bayesain dynamic factor model by MCMC methods using Durbin and Koopman's
#' disturbance smoother. This function may be called directly by advanced
#' users. The only dependencies are the  Armadillo
#' (\url{http://arma.sourceforge.net/}) linear algebra library for C++ and the
#' packages needed for interfacing with R (\code{\link{Rcpp}} and
#' \code{\link{RcppArmadillo}}).
#'
#' @param B  initial guess for B in transition equation
#' @param Bp prior for B
#' @param Jb Helper matrix for transition equation, identity matrix if uniform frequency
#' @param lam_B prior tightness for B (additive)
#' @param q initial guess for q in the transition equation
#' @param nu_q prior "degrees of freedom" for inverse-Whishart prior for q (additive, prior scale is fixed so that increasing nu_q shrinks the variance towards zero)
#' @param H initial guess for H in the trasition equation
#' @param Hp prior for H
#' @param lam_H prior tightness for H (additive)
#' @param R initial guess for diagonal elements of R in the transition equation, entered as a vector
#' @param nu_r prior deg. of freedom for elements of R, entered as a vector (additive, prior scale is fixed so that increasing nu_r[j] shrinks the variance of shocks to series j towards zero)
#' @param Y Input data. Data must be scaled and centered prior to estimation if desired.
#' @param freq vector, number of high frequency periods in an observation
#' @param LD vector, 0 for level data and 1 for differenced data
#' @param Ystore T/F, should the distribution of Y be stored
#' @param store_idx, if Ystore is TRUE, index of which observed series to store. Note C++ uses zero indexing (i.e. subtract 1 from the R index value)
#' @param reps number of repetitions for MCMC sampling
#' @param burn number of iterations to burn in MCMC sampling
#' @param verbose print status of function during evaluation.
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib bdfm
Cppbdfm <- function(B, Bp, Jb, lam_B, q, nu_q, H, Hp, lam_H, R, nu_r, Y, freq, LD, Ystore = FALSE, store_idx = 0, reps = 1000, burn = 500, verbose = FALSE) {
  OUT <- EstDFM(B = B, Bp = Bp, Jb = Jb, lam_B = lam_B, q = q, nu_q = nu_q, H = H, Hp = Hp, lam_H = lam_H, R = R, nu_r = nu_r, Y = Y, freq = freq, LD = LD, reps = reps, burn = burn, verbose = verbose)
  return(OUT)
}
