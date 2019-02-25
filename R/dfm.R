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
      scale = scale, logs = logs, diffs = diffs, outlier_threshold = outlier_threshold, freq = frequency_mix,
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
      data_unclassed <- as.matrix(unclass(data))
    } else {
      stopifnot(requireNamespace("tsbox"))
      # time series
      stopifnot(tsbox::ts_boxable(data))
      # convert to mts
      data_unclassed <- unclass(tsbox::ts_ts(data))
    }
    data_tsp <- attr(data_unclassed, "tsp")
    attr(data_unclassed, "tsp") <- NULL
    ans <- dfm_core(
      Y = data_unclassed, m = factors, p = lags, FC = forecasts, method = method,
      scale = scale, logs = logs, diffs = diffs, outlier_threshold = outlier_threshold, freq = frequency_mix,
      preD = pre_differenced, Bp = trans_prior, lam_B = trans_shrink, trans_df = trans_df,
      Hp = obs_prior, lam_H = obs_shrink, obs_df = obs_df,
      ID = identification, store_idx = store_idx, reps = reps,
      burn = burn, verbose = verbose, tol = tol
    )

    # make values a ts timeseries
    ans$values <- ts(ans$values, start = data_tsp[1], frequency = data_tsp[3])
    ans$factors <- ts(ans$factors, start = data_tsp[1], frequency = data_tsp[3])
    if (!is.null(store_idx)) {
      ans$Ymedian <- ts(ans$Ymedian, start = data_tsp[1], frequency = data_tsp[3])
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



