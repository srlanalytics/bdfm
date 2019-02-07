#' Estimate dynamic factor model
#'
#' @param Y data in matrix format with time in rows
#' @param factors number of factors
#' @param lags number of lags in transition equation
#' @param forecast number of periods ahead to forecast
#' @param method character, method to be used
#' @param frequency_mix numeric, number of high frequency periods in observation if data is mixed frequency
#' @param differences, numeric, 0 for levels, 1 for first differences, only specified for mixed frequency models
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
#' summary(m)
#' @useDynLib bdfm
dfm <- function(data, factors = 1, lags = 3, forecasts = 0,
                method = c("bayesian", "ml", "pc"), scale = TRUE, logs = NULL, diffs = NULL,
                frequency_mix = "auto", as_differenced = NULL, trans_prior = NULL,
                trans_shrink = 0, trans_df = 0, obs_prior = NULL, obs_shrink = 0,
                obs_df = NULL, identification = "PC_full",
                store_idx = NULL, reps = 1000, burn = 500, loud = FALSE,
                EM_tolerance = 0.01) {
  method <- match.arg(method) # checks and picks the first if unspecified

  if (!is.null(frequency_mix) && method != "bayesian") {
    stop("Mixed freqeuncy models are only supported for Bayesian estimation")
  }

  tsobjs <- c(
    "zoo", "xts", "tslist", "tbl_ts", "timeSeries",
    "tbl_time", "tbl_df", "data.table", "data.frame", "dts"
  )

  if (!requireNamespace("tsbox") & any(class(Y) %in% tsobjs)) {
    stop('"tsbox" is needed to support non ts-time-series. To install: \n\n  install.packages("tsbox")', call. = FALSE)
  }

  # non time series
  if (!any(class(data) %in% c(tsobjs, "ts", "mts")) && is.matrix(data)) {
    ans <- dfm_core(
      Y = data, m = factors, p = lags, FC = forecast, method = method,
      scale = scale, logs = logs, diffs = diffs, freq = frequency_mix,
      asD = as_differenced, Bp = trans_prior, lam_B = trans_shrink, nu_q = trans_df,
      Hp = obs_prior, lam_H = obs_shrink, nu_r = obs_df,
      ID = identification, store_idx = store_idx, reps = reps,
      burn = burn, loud = loud, tol = EM_tolerance
    )
    colnames(ans$values) <- colnames(data)
    ans$dates <- NULL
  } else {

    # no need for tsbox if Y is ts or mts
    if (inherits(data, "ts")) {
      Y.uc <- unclass(data)
    } else {
      stopifnot(requireNamespace("tsbox"))
      # time series
      stopifnot(tsbox::ts_boxable(data))
      # convert to mts
      Y.uc <- unclass(tsbox::ts_ts(data))
    }

    Y.tsp <- attr(Y.uc, "tsp")
    attr(Y.uc, "tsp") <- NULL
    ans <- dfm_core(
      Y = Y.uc, m = factors, p = lags, FC = forecast, method = method,
      scale = scale, logs = logs, diffs = diffs, freq = frequency_mix,
      asD = as_differenced, Bp = trans_prior, lam_B = trans_shrink, nu_q = trans_df,
      Hp = obs_prior, lam_H = obs_shrink, nu_r = obs_df,
      ID = identification, store_idx = store_idx, reps = reps,
      burn = burn, loud = loud, tol = EM_tolerance
    )

    # make values a ts timeseries
    ans$values <- ts(ans$values, start = Y.tsp[1], frequency = Y.tsp[3])
    ans$factors <- ts(ans$factors, start = Y.tsp[1], frequency = Y.tsp[3])
    if (!is.null(store_idx)) {
      ans$Ymedian <- ts(ans$Ymedian, start = Y.tsp[1], frequency = Y.tsp[3])
    }
    # apply column names from Y
    colnames(ans$values) <- colnames(data)

    # return a date vector
    ans$dates <- tsbox::ts_regular(tsbox::ts_df(ans$values))[, 1]

    # put values back into original class
    if (!inherits(Y, "ts")) {
      ans$values <- tsbox::copy_class(ans$values, data)
      ans$factors <- tsbox::copy_class(ans$factors, data, preserve.mode = FALSE)
      if (!is.null(store_idx)) {
        ans$Ymedian <- tsbox::copy_class(ans$Ymedian, data)
      }
    }
  }
  class(ans) <- "dfm"
  return(ans)
}

dfm_core <- function(Y, m, p, FC, method, scale, logs, diffs, freq, asD,
                     Bp, lam_B, nu_q, Hp, lam_H, nu_r, ID,
                     store_idx, reps, burn, loud, tol) {

  #-------Data processing-------------------------

  # frequency
  if (freq == "auto") {
    freq <- apply(Y, MARGIN = 2, FUN = get_freq)
  } else if (!is.integer(freq) || length(freq) != ncol(Y)) {
    stop("Argument 'freq' must be 'auto' or integer valued with 
         length equal to the number data series")
  }

  # logs
  if (!is.null(logs)) {
    if (is.character(logs)) {
      logs <- unlist(sapply(logs, FUN = grep, colnames(Y)))
    } else if (!is.numeric(logs)) {
      stop("Argument 'logs' must be either a character (string) vector or numeric index values")
    }
    Y[, logs] <- log(Y[, logs])
  }

  # differences
  if (!is.null(diffs)) {
    Y_lev <- Y
    if (is.character(diffs)) {
      diffs <- unlist(sapply(diffs, FUN = grep, colnames(Y)))
    } else if (!is.numeric(diffs)) {
      stop("Argument 'diffs' must be either a character (string) vector or numeric index values")
    }
    Y[, diffs] <- sapply(diffs, mf_diff, fq = freq, Y = Y)
  }

  # specify which series are differenced for mixed frequency estimation
  LD <- rep(0, ncol(Y))
  if (!is.null(asD)) {
    if (is.character(asD)) {
      preD <- unlist(sapply(asD, FUN = grep, colnames(Y)))
    } else if (!is.numeric(asD)) {
      stop("Argument 'asD' must be either a character (string) vector or numeric index values")
    }
  } else {
    LD[diffs] <- 1 # in bdfm 1 indicates differenced data, 0 level data
  }


  if (scale) {
    Y <- bdfm::scale(Y)
  }

  if (method == "bayesian") {
    est <- bdfm(
      Y = Y, m = m, p = p, FC = FC, Bp = Bp,
      lam_B = lam_B, Hp = Hp, lam_H = lam_H, nu_q = nu_q, nu_r = nu_r,
      ID = ID, store_idx = store_idx, freq = freq, LD = LD, reps = reps,
      burn = burn, loud = loud
    )
  } else if (method == "ml") {
    est <- MLdfm(
      Y = Y, m = m, p = p, FC = FC, tol = tol,
      loud = loud
    )
  } else if (method == "pc") {
    est <- PCdfm(
      Y,
      m = m, p = p, FC = FC, Bp = Bp,
      lam_B = lam_B, Hp = Hp, lam_H = lam_H, nu_q = nu_q, nu_r = nu_r,
      ID = ID, reps = reps, burn = burn
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
    nrow(x$B), "factor(s) and", ncol(x$B) / nrow(x$B), "lag(s)."
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
    "factor(s) and", ncol(object$B) / nrow(object$B), "lag(s)."
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
  colnames(H) <- as.character(seq(1, ncol(H)))
  print(H)
  cat("\n shocks R: \n")
  r <- data.frame(diag(object$R))
  row.names(r) <- colnames(object$values)
  colnames(r) <- "Variance of Shocks"
  print(r)
}
