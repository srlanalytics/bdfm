#' Estimate Bayesian dynamic factor model
#'
#' @param Y Data in matrix format with time in rows
#' @param factors number of factors
#' @param lags number of lags in transition equation
#' @param forecast number of periods ahead to forecast
#' @param method character, method to be used
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom stats dnorm na.omit ts var
#' @importFrom utils head tail
#' @importFrom stats dnorm na.omit ts var
#' @importFrom utils head tail
#' @useDynLib BDFM
dfm <- function(Y, factors = 1, lags = 2, forecast = 0, method = "Bayesian") {

  # non time series
  if (!ts_boxable(Y) && is.matrix(Y)) {
      E <- BDFM(Y = Y, m = factors, p = lags, FC = forecast, Bp = B_prior, lam_B = lam_B, Hp = H_prior, lam_H = lam_H, nu_q = nu_q, nu_r = nu_r, ID = ID, ITC = intercept, reps = reps, burn = burn)
      colnames(E$values) <- colnames(Y)
      return(E)
  }else{
    # time series
    stopifnot(ts_boxable(Y))
    # convert to mts
    Y.uc  <- unclass(ts_ts(Y))
    Y.tsp <- attr(Y.uc, "tsp")
    attr(Y.uc, "tsp") <- NULL

    if(method == "Bayesian"){
      B_prior   <- getOption('B_prior', default = NULL)
      lam_B     <- getOption('lam_B', default = 0)
      H_prior   <- getOption('H_prior', default = NULL)
      lam_H     <- getOption('lam_H', default = 0)
      nu_q      <- getOption('nu_q', default = 0)
      nu_r      <- getOption('nu_r', default = NULL)
      ID        <- getOption('ID', default = "PC_full")
      intercept <- getOption('intercept', default = T)
      reps      <- getOption('reps', default = 1000)
      burn      <- getOption('burn', default = 500)
      E <- BDFM(Y = Y.uc, m = factors, p = lags, FC = forecast, Bp = B_prior, lam_B = lam_B, Hp = H_prior, lam_H = lam_H, nu_q = nu_q, nu_r = nu_r, ID = ID, ITC = intercept, reps = reps, burn = burn)
    }else if(method == "ML"){
      tol       <- getOption('tol', default = 0.01)
      Loud      <- getOption('Loud', default = FALSE)
      E <- MLdfm(Y, m = factors, p = lags, FC = forecast, tol = tol, Loud = Loud)
    }else if(method == "PC"){
      B_prior   <- getOption('B_prior', default = NULL)
      lam_B     <- getOption('lam_B', default = 0)
      H_prior   <- getOption('H_prior', default = NULL)
      lam_H     <- getOption('lam_H', default = 0)
      nu_q      <- getOption('nu_q', default = 0)
      nu_r      <- getOption('nu_r', default = NULL)
      ID        <- getOption('ID', default = "PC_full")
      intercept <- getOption('intercept', default = T)
      reps      <- getOption('reps', default = 1000)
      burn      <- getOption('burn', default = 500)
      E <- PCdfm(Y.uc, m = factors, p = lags, FC = forecast, Bp = B_prior, lam_B = lam_B, Hp = H_prior, lam_H = lam_H, nu_q = nu_q, nu_r = nu_r, ID = ID, ITC = intercept, reps = reps, burn = burn)
    }else{
      stop("method must be either Bayesian, ML, or PC")
    }

    ts(E$values, start = Y.tsp[1], frequency = Y.tsp[3])  #make predicted values a ts timeseries
    ts(E$factors, start = Y.tsp[1], frequency = Y.tsp[3]) #make estimated factors a ts timeseries
    colnames(E$values) <- colnames(Y) #apply column names from Y

    E$values  <- copy_class(E$values, Y)  #put predicted values back into original class
    E$factors <- copy_class(E$factors, Y) #put estimated factors back into original class
  }
  class(E) <- "dfm"
  return(E)
}







# methods
#' @export
#' @method predict dfm
predict.dfm <- function(object, ...) {
  object$values
}

#' @export
#' @method print dfm
print.dfm <- function(x, ...) {
  cat("Call: \n Bayesian dynamic factor model with", nrow(x$B), "factor(s) and", ncol(x$B)/nrow(x$B), "lag(s).")
  cat("\n \n")
  cat("Log Likelihood:", x$Lik)
  cat("\n \n")
  cat("BIC:", x$BIC)
}

#' @export
#' @method summary dfm
summary.dfm <- function(object, ...) {
  cat("Call: \n Bayesian dynamic factor model with", nrow(object$B), "factor(s) and", ncol(object$B)/nrow(object$B), "lag(s).")
  cat("\n \n")
  cat("Log Likelihood:", object$Lik)
  cat("\n \n")
  cat("BIC:", object$BIC)
  cat("\n \n")
  cat("Posterior medians for transition equation: \n")
  cat("\n Coefficients B: \n")
  print(object$B)
  cat("\n Covariance Q: \n")
  print(object$q)
  cat("\n \n")
  cat("Posterior medians for observation equation: \n")
  cat("\n Coefficients H: \n")
  H <- data.frame(object$H)
  row.names(H) <- colnames(object$values)
  colnames(H) <- as.character(seq(1,ncol(H)))
  print(H)
  cat("\n shocks R: \n")
  r <- data.frame(diag(object$R))
  row.names(r) <- colnames(object$values)
  colnames(r) <- "Variance of Shocks"
  print(r)

}
