#' Estimate dynamic factor model
#'
#' @param Y data in matrix format with time in rows
#' @param factors number of factors
#' @param lags number of lags in transition equation
#' @param forecast number of periods ahead to forecast
#' @param method character, method to be used
#' @param B_prior prior matrix for B in the transition equation. Default is
#'   zeros.
#' @param lam_B prior tightness on B
#' @param H_prior prior matrix for H (loadings) in the observation equation.
#'  Default is zeros.
#' @param lam_H prior tightness on H
#' @param nu_q prior deg. of freedom for transition equation, entered as vector
#'   with length equal to the number of factors.
#' @param nu_r prior deg. of freedom for observables, entered as vector with
#'   length equal to the number of observables.
#' @param identification factor identification. 'PC_full' is the default (using
#'   all observed series), 'PC_sub' finds a submatrix of the data that maximizes
#'   the number of observations for a square (no missing values) data set. Use
#'   'PC_sub' when many observations are missing.
#' @param intercept logical, should an icept be included?
#' @param store_idx, if estimation is Bayesian, index of input data to store the full posterior distribution of predicted values.
#' @param reps number of repetitions for MCMC sampling
#' @param burn number of iterations to burn in MCMC sampling
#' @param loud print status of function during evalutation. If ML, print
#'   difference in likelihood at each iteration of the EM algorithm.
#' @param EM_tolerance tolerance for convergence of EM algorithm (method `ml` only). Convergence
#'   criteria is calculated as 200 * (Lik1 - Lik0) / abs(Lik1 + Lik0) where Lik1
#'   is the log likelihood from this iteration and Lik0 is the likelihood from
#'   the previous iteration.
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom stats dnorm na.omit ts var
#' @importFrom utils head tail
#' @examples
#' fdeaths0 <- fdeaths
#' fdeaths0[length(fdeaths0)] <- NA
#' dta <- cbind(fdeaths0, mdeaths)
#'
#' library(bdfm)
#' m <- dfm(dta, forecast = 2)
#' predict(m)
#'
#' @useDynLib bdfm
dfm <- function(Y, factors = 1, lags = 2, forecast = 0,
                method = c("bayesian", "ml", "pc"),
                B_prior = NULL, lam_B = 0, H_prior = NULL, lam_H = 0, nu_q = 0,
                nu_r = NULL, identification = "PC_full", intercept = TRUE,
                store_idx = NULL, reps = 1000, burn = 500, loud = FALSE,
                EM_tolerance = 0.01) {

  method <- match.arg(method)  # checks and picks the first if unspecified

  tsobjs <- c("zoo", "xts", "tslist", "tbl_ts", "timeSeries",
    "tbl_time", "tbl_df", "data.table", "data.frame", "dts")

  if (!requireNamespace("tsbox") & any(class(Y) %in% tsobjs)) {
    stop('"tsbox" is needed to support non ts-time-series. To install: \n\n  install.packages("tsbox")', call. = FALSE)
  }

  # non time series
  if (!any(class(Y) %in% c(tsobjs, "ts", "mts")) && is.matrix(Y)) {

    ans <- dfm_core(
      Y = Y, m = factors, p = lags, FC = forecast, method = method, Bp = B_prior,
      lam_B = lam_B, Hp = H_prior, lam_H = lam_H, nu_q = nu_q, nu_r = nu_r,
      ID = identification, ITC = intercept, store_idx = store_idx, reps = reps,
      burn = burn, loud = loud, tol = EM_tolerance
    )
    colnames(ans$values) <- colnames(Y)

  } else {

    # no need for tsbox if Y is ts or mts
    if (inherits(Y, "ts")) {
      Y.uc  <- unclass(Y)
    } else {
      stopifnot(requireNamespace("tsbox"))
      # time series
      stopifnot(tsbox::ts_boxable(Y))
      # convert to mts
      Y.uc  <- unclass(tsbox::ts_ts(Y))
    }

    Y.tsp <- attr(Y.uc, "tsp")
    attr(Y.uc, "tsp") <- NULL
    ans <- dfm_core(
      Y = Y.uc, m = factors, p = lags, FC = forecast, method = method, Bp = B_prior,
      lam_B = lam_B, Hp = H_prior, lam_H = lam_H, nu_q = nu_q, nu_r = nu_r,
      ID = identification, ITC = intercept, store_idx = store_idx, reps = reps,
      burn = burn, loud = loud, tol = EM_tolerance
    )

    # make values a ts timeseries
    ans$values <- ts(ans$values, start = Y.tsp[1], frequency = Y.tsp[3])
    ans$factors <- ts(ans$factors, start = Y.tsp[1], frequency = Y.tsp[3])
    if(!is.null(store_idx)){
      ans$Ymedian <- ts(ans$Ymedian, start = Y.tsp[1], frequency = Y.tsp[3])
    }
    # apply column names from Y
    colnames(ans$values) <- colnames(Y)

    # put values back into original class
    if (!inherits(Y, "ts")) {
      ans$values  <- tsbox::copy_class(ans$values, Y)
      ans$factors <- tsbox::copy_class(ans$factors, Y, preserve.mode = FALSE)
      if(!is.null(store_idx)){
         ans$Ymedian <- tsbox::copy_class(ans$Ymedian, Y)
      }
    }
  }
  class(ans) <- "dfm"
  return(ans)
}

dfm_core <- function(Y, m, p, FC, Bp, lam_B, Hp, lam_H, nu_q, nu_r, ID, ITC,
                     store_idx, reps, burn, tol, loud, method) {
  if (method == "bayesian") {
    bdfm(
      Y = Y, m = m, p = p, FC = FC, Bp = Bp,
      lam_B = lam_B, Hp = Hp, lam_H = lam_H, nu_q = nu_q, nu_r = nu_r,
      ID = ID, ITC = ITC, store_idx = store_idx, reps = reps, burn = burn,
      loud = loud
    )
  } else if (method == "ml") {
    MLdfm(
      Y = Y, m = m, p = p, FC = FC, tol = tol,
      loud = loud
    )
  } else if (method == "pc") {
    PCdfm(
      Y, m = m, p = p, FC = FC, Bp = Bp,
      lam_B = lam_B, Hp = Hp, lam_H = lam_H, nu_q = nu_q, nu_r = nu_r,
      ID = ID, ITC = ITC, reps = reps, burn = burn
    )
  }
}


#' extractor function for factors
#' @param x obeject of class `"dfm"`
#' @export
factors <- function(x) {
  stopifnot(inherits(x, "dfm"))
  x$factors
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
  cat(
    "Call: \n Bayesian dynamic factor model with",
    nrow(x$B), "factor(s) and", ncol(x$B)/nrow(x$B), "lag(s)."
  )
  cat("\n \n")
  cat("Log Likelihood:", x$Lik)
  cat("\n \n")
  cat("BIC:", x$BIC)
}

#' @export
#' @method summary dfm
summary.dfm <- function(object, ...) {
  cat(
    "Call: \n Bayesian dynamic factor model with", nrow(object$B),
    "factor(s) and", ncol(object$B)/nrow(object$B), "lag(s)."
  )
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
