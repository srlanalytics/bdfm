#' Estimate dynamic factor model
#'
#' @param data data in matrix format with time in rows
#' @param factors number of factors
#' @param lags number of lags in transition equation
#' @param forecasts number of periods ahead to forecasts
#' @param method character, method to be used
#' @param scale scale data before estimation (True/False)?
#' @param logs names or index values of series which should be entered in log levels
#' @param diffs names or index values of series which should be differenced
#' @param outlier_threshold drop observations more than x standard deviations from the mean
#' @param frequency_mix 'auto' or numeric, number of high frequency periods in observation if data is mixed frequency
#' @param pre_differenced, names or index values of low freqeuncy series entered in differences (not necessary if already specified in 'diffs')
#' @param trans_prior prior matrix for B in the transition equation. Default is
#'   zeros.
#' @param trans_shrink prior tightness on B matrix in trasition equation
#' @param trans_df prior deg. of freedom for transition equation
#' @param obs_prior prior matrix for H (loadings) in the observation equation
#'  Default is zeros.
#' @param obs_shrink prior tightness on H (loadings) in the observation equation
#' @param obs_df prior deg. of freedom for observables, entered as vector with
#'   length equal to the number of observables.
#' @param identification factor identification. 'pc_long' is the default and finds series with the most
#' observations over time. 'pc_full' uses all observed series, 'pc_sub' finds a submatrix of the data that maximizes
#'   the number of observations for a square (no missing values) data set. Users may also enter a
#'   numeric vector for specified series.
#' @param store_idx, if estimation is Bayesian, index of input data to store the full posterior distribution of predicted values.
#' @param reps number of repetitions for MCMC sampling
#' @param burn number of iterations to burn in MCMC sampling
#' @param verbose print status of function during evalutation. If ML, print
#'   difference in likelihood at each iteration of the EM algorithm. Default is `TRUE` in interactive mode, `FALSE` otherwise.
#' @param tol tolerance for convergence of EM algorithm (method `ml` only). Convergence
#'   criteria is calculated as 200 * (Lik1 - Lik0) / abs(Lik1 + Lik0) where Lik1
#'   is the log likelihood from this iteration and Lik0 is the likelihood from
#'   the previous iteration.
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom stats dnorm na.omit ts var approx frequency is.ts loess median model.matrix na.exclude predict setNames start
#' @importFrom utils head tail
#' @examples
#' fdeaths0 <- fdeaths
#' fdeaths0[length(fdeaths0)] <- NA
#' dta <- cbind(fdeaths0, mdeaths)
#'
#' library(bdfm)
#' m <- dfm(dta, forecasts = 2)
#' summary(m)
#' @useDynLib bdfm
dfm <- function(data, factors = 1, lags = "auto", forecasts = "auto",
                method = c("bayesian", "ml", "pc"), scale = TRUE, logs = NULL, diffs = NULL,
                outlier_threshold = 4, frequency_mix = "auto", pre_differenced = NULL,
                trans_prior = NULL, trans_shrink = 0, trans_df = 0, obs_prior = NULL, obs_shrink = 0,
                obs_df = NULL, identification = "pc_long",
                store_idx = NULL, reps = 1000, burn = 500, verbose = interactive(),
                tol = 0.01) {

  call <- match.call

  method <- match.arg(method) # checks and picks the first if unspecified

  if (!is.null(frequency_mix) && method != "bayesian") {
    stop("Mixed freqeuncy models are only supported for Bayesian estimation")
  }

  # check need for tsbox
  tsobjs <- c(
    "zoo", "xts", "tslist", "tbl_ts", "timeSeries",
    "tbl_time", "tbl_df", "data.table", "data.frame", "dts"
  )
  if (any(class(data) %in% tsobjs) && !requireNamespace("tsbox")) {
    stop('"tsbox" is needed to support non ts-time-series. To install: \n\n  install.packages("tsbox")', call. = FALSE)
  }

  # non time series
  if (!any(class(data) %in% c(tsobjs, "ts", "mts")) && is.matrix(data)) {
    ans <- dfm_core(
      Y = data, m = factors, p = lags, FC = forecasts, method = method,
      scale = scale, logs = logs, diffs = diffs, freq = frequency_mix,
      preD = pre_differenced, Bp = trans_prior, lam_B = trans_shrink, trans_df = trans_df,
      Hp = obs_prior, lam_H = obs_shrink, obs_df = obs_df,
      ID = identification, store_idx = store_idx, reps = reps,
      burn = burn, verbose = verbose, tol = tol
    )
    colnames(ans$values) <- colnames(data)
    ans$dates <- NULL
  } else {
    # no need for tsbox if Y is ts or mts
    if (inherits(data, "ts")) {
      Y.uc <- as.matrix(unclass(data))
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
      Y = Y.uc, m = factors, p = lags, FC = forecasts, method = method,
      scale = scale, logs = logs, diffs = diffs, outlier_threshold = outlier_threshold, freq = frequency_mix,
      preD = pre_differenced, Bp = trans_prior, lam_B = trans_shrink, trans_df = trans_df,
      Hp = obs_prior, lam_H = obs_shrink, obs_df = obs_df,
      ID = identification, store_idx = store_idx, reps = reps,
      burn = burn, verbose = verbose, tol = tol
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
    # ans$dates <- tsbox::ts_regular(tsbox::ts_df(ans$values))[, 1]

    # put values back into original class
    if (!inherits(data, "ts")) {
      ans$values <- tsbox::copy_class(ans$values, data)
      ans$factors <- tsbox::copy_class(ans$factors, data, preserve.mode = FALSE)
      if (!is.null(store_idx)) {
        ans$Ymedian <- tsbox::copy_class(ans$Ymedian, data)
      }
    }
  }

  ans$call <- match.call()
  class(ans) <- "dfm"
  return(ans)
}

# m <- 3
# p <- "auto"
# freq <- "auto"
# Bp <- NULL
# preD <- 1
# lam_B = 0
# nu_q = 0
# Hp = NULL
# lam_H = 0
# obs_df = NULL
# ID = "pc_long"
# store_idx = 2
# reps = 1000
# burn = 500
# verbose = T
# tol = .01
# FC = 3
# logs = c( 2,  4,  5,  8,  9, 10, 11, 12, 15, 16, 17, 21, 22)
# diffs = c(2, 4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24)
# outlier_threshold <- 4
# scale = TRUE

dfm_core <- function(Y, m, p, FC, method, scale, logs, outlier_threshold, diffs, freq, preD,
                     Bp, lam_B, trans_df, Hp, lam_H, obs_df, ID,
                     store_idx, reps, burn, verbose, tol) {



  #-------Data processing-------------------------

  k <- NCOL(Y) # number of series

  # frequency
  if (freq == "auto") {
    freq <- apply(Y, MARGIN = 2, FUN = get_freq)
  } else if (!is.integer(freq) || length(freq) != ncol(Y)) {
    stop("Argument 'freq' must be 'auto' or integer valued with
         length equal to the number data series")
  }

  if (p == "auto"){
    p <- max(freq)
  }

  if (FC == "auto"){
    FC <- max(freq)
  }

  # add forecast periods
  if (FC > 0) {
    tmp <- matrix(NA, FC, k)
    Y <- rbind(Y, tmp)
  }

  # logs
  if (!is.null(logs)) {
    logs <- standardize_index(logs, Y)
    Y[, logs] <- log(Y[, logs])
  }

  # differences
  if (!is.null(diffs)) {
    Y_lev <- Y
    diffs <- standardize_index(diffs, Y)
    Y[, diffs] <- sapply(diffs, mf_diff, fq = freq, Y = Y)
  }

  # if (!is.null(trans_df)) {
  #   trans_df <- standardize_index(trans_df, Y)
  # }

  if (!is.null(obs_df)) {
    obs_df <- standardize_numeric(obs_df, Y)
  }

  # specify which series are differenced for mixed frequency estimation

  LD <- rep(0, NCOL(Y))
  if (!is.null(preD)) {
    preD <- standardize_index(preD, Y)
  }
  LD[unique(c(preD, diffs))] <- 1 # in bdfm 1 indicates differenced data, 0 level data

  # drop outliers
  Y[abs(scale(Y)) > outlier_threshold] <- NA

  if (scale) {
    Y <- 100*scale(Y)
    y_scale  <- attr(Y, "scaled:scale")
    y_center <- attr(Y, "scaled:center")
  }

  if (!is.null(store_idx)){
    if(length(store_idx)>1){
      stop("Length of 'store_idx' cannot be greater than 1")
    }
    store_idx <- standardize_index(store_idx, Y)
  }

  if(!ID%in%c("pc_full", "pc_sub", "pc_long", "name")){
    ID <- standardize_index(ID, Y)
  }

  if (method == "bayesian") {
    est <- bdfm(
      Y = Y, m = m, p = p, Bp = Bp,
      lam_B = lam_B, Hp = Hp, lam_H = lam_H, nu_q = trans_df, nu_r = obs_df,
      ID = ID, store_idx = store_idx, freq = freq, LD = LD, reps = reps,
      burn = burn, verbose = verbose
    )
  } else if (method == "ml") {
    est <- MLdfm(
      Y = Y, m = m, p = p, tol = tol,
      verbose = verbose
    )
  } else if (method == "pc") {
    est <- PCdfm(
      Y, m = m, p = p, Bp = Bp,
      lam_B = lam_B, Hp = Hp, lam_H = lam_H, nu_q = trans_df, nu_r = obs_df,
      ID = ID, reps = reps, burn = burn
    )
  }

  # undo scaling
  if(scale){
    est$values <- (matrix(1, nrow(est$values), 1) %x% t(y_scale)) * (est$values / 100) + (matrix(1, nrow(est$values), 1) %x% t(y_center))
    if(!is.null(store_idx) && method == "bayesian"){
      est$Ystore <- est$Ystore*(y_scale[store_idx]/100) + y_center[store_idx]
      est$Ymedain <- est$Ymedian*(y_scale[store_idx]/100) + y_center[store_idx]
    }
  }

  # undo differences
  if (!is.null(diffs)) {
    est$values[,diffs] <- sapply(diffs, FUN = level, fq = freq, Y_lev = Y_lev, vals = est$values)
    if(!is.null(store_idx) && method == "bayesian" && store_idx%in%diffs){
      est$Ymedain <- level_simple(est$Ymedain, y_lev = Y_lev[,store_idx], fq = freq[store_idx])
      est$Ystore  <- apply(est$Ystore, MARGIN = 2, FUN = level_simple, y_lev = Y_lev[,store_idx], fq = freq[store_idx])
    }
  }

  # undo logs
  if (!is.null(logs)) {
    est$values[,logs] <- exp(est$values[,logs])
    if(!is.null(store_idx) && method == "bayesian" && store_idx%in%logs){
      est$Ymedain <- exp(est$Ymedain)
      est$Ystore  <- exp(est$Ystore)
    }
  }

  return(est)
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

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  cat("Bayesian dynamic factor model with",
    nrow(x$B), "factor(s) and", ncol(x$B) / nrow(x$B), "lag(s)."
  )
  cat("\n")
  cat("Log Likelihood:", x$Lik, " ")
  cat("BIC:", x$BIC)
  cat("\n")
  invisible(x)
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
