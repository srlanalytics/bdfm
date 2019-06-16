#' Estimate a Dynamic Factor Model
#'
#' Estimates a Bayesian or non-Bayesian dynamic factor Model. With the default
#' options, `dfm` calls  automatic procedures that works well in many
#' circumstances.
#'
#' **Specifying series**: Individual series can be specified either by *names*
#' (recommended) or index values. An index value refers to the position of the
#' series in `data`.
#'
#' **Specifying parameters for specific series**: Parameters for individual
#' series can be specified using a *named* vector (recommended) or using a
#' unnamed vector of the same length as as the number of series in `data`.
#'
#' @param data one or multiple time series. The data to be used for estimation.
#'   This can be entered as a `"ts"` object of as a matrix. If
#'   [tsbox](https://tsbox.help) is installed, any ts-boxable time series can be
#'   supplied (`ts`, `xts`, `zoo`, `data.frame`, `data.table`, `tbl`, `tbl_ts`,
#'   `tbl_time`, or `timeSeries`)
#' @param factors integer. The number of unobserved factors to be estimated. A
#'   larger number of factors leads to a more complex model.
#' @param lags integer. The number of lags in the transition equation. If
#'  `"auto"` (default), the number is equal to highest frequency in `data`.
#' @param forecasts integer. Number of periods ahead to forecasts.
#' @param method character. Method to be used; one of `"bayesian"`, `"ml"` or
#'   `"pc". See details.
#' @param scale logical. Should data be scaled before estimation? `TRUE`
#'   (default) resolves some numerical problems during estimation. `FALSE`
#'   ensures that the coefficient estimates are interpretable.
#' @param logs names or index values (see details). Series of which the
#'   logarithm is taken (can be combined with `diffs`). If `"auto"` (default)
#'   this is done for all series that are differentiated and have no values < 0.
#' @param diffs names or index values (see details). Series to be
#'   differentiated. If `"auto"` (default), a modified Durbin-Watson test is
#'   performed.
#' @param outlier_threshold integer. Observations more than `outlier_threshold`
#'   standard deviations from the series mean are removed. This is useful to
#'   increase the stability of the estimation.
#' @param frequency_mix integer or `"auto"`. Number of high frequency periods
#'   in a low frequency period. If `"auto"` (default), this is inferred from the
#'   time series.
#' @param pre_differenced names or index values (see details). series entered in
#'   differences (If series are specified in 'diffs', this is not needed.)
#' @param trans_prior prior matrix for B in the transition equation. Default is
#'   zeros. To use a random walk prior with, for example,
#'   m factors and p lags, set `trans_prior = cbind(diag(1,m,m), matrix(0,m,m*(p-1)))`.
#' @param trans_shrink prior tightness on B matrix in transition equation. Used to shrink
#'   forecast values towards the prior `trans_prior`,
#'   which may help reduce parameter uncertainty.
#' @param trans_df prior degrees of freedom for transition equation. 
#'   Shrinking shocks to the trasition equation will increase the magnitude
#'   of shocks to the observation equation dampening updates from observed series (method `bayesian` only). 
#' @param obs_prior prior matrix for H (loadings) in the observation equation
#'  Default is zeros.
#' @param obs_shrink prior tightness on H (loadings) in the observation equation; a greater value will shrink estimates of loadings
#'  more aggressively towards the prior 'obs_prior'. When the prior is zero (the default value), this is an alternative
#'  (and typically more stable) approach to dampening the impact of updates from observed series. 
#' @param obs_df named vector (see details). prior degrees of freedom
#'   for gamma distribution in the observation equation. This is useful to give specific series a larger weight,
#'   e.g. 1. (default 0, method `bayesian` only).
#' @param identification  names or index values (see details), or character. Factor identification. `"pc_long"`
#'   (default) finds series with the most observations over time, on which it uses principal components. '"pc_full"'
#'   uses all observed series, `"pc_sub"` finds a submatrix of the data that
#'   maximizes the number of observations for a square (no missing values) data
#'   set. Identification can also be done manually, by supplying names or index
#'   values.
#' @param keep_posterior names or index values (see details). Series of which to keep the full posterior distribution of predicted values (method `bayesian` only). This is useful for forecasting.
#' @param interpolate logical. If data is mixed frequency, should low frequency be interpolated?
#' @param orthogonal_shocks return a rotation of the model with orthogonal shocks and factors. This is useful ....
#' @param reps number of repetitions for MCMC sampling
#' @param burn number of iterations to burn in MCMC sampling
#' @param verbose print status of function during evaluation. If ML, print
#'   difference in likelihood at each iteration of the EM algorithm. Default is
#'  `TRUE` in interactive mode, `FALSE` otherwise, so it does not appear, e.g., in `reprex::reprex()`.
#'
#' @param tol tolerance for convergence of EM algorithm (method `ml` only).
#' @seealso `vignette("dfm")`, for a more comprehensive intro to the package.
#' @seealso [*Practical Implementation of Factor Models*](http://srlquantitative.com/docs/Factor_Models.pdf) for a comprehensive overview of dynamic factor models.
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom stats dnorm na.omit ts var approx frequency is.ts loess median model.matrix na.exclude predict setNames start
#' @importFrom utils head tail
#' @examples
#'
#' dta <- cbind(fdeaths, mdeaths)
#'
#' m0 <- dfm(dta, forecast = 2) # estimation with 2 period forecast
#' predict(m0)                  # series with imputations and forecasts
#' summary(m0)                  # summary of the model
#' factors(m0)                  # estimated factor
#'
#' # informative priors: giving 'fdeaths' a higher weight
#' m1 <- dfm(dta, obs_df = c("fdeaths" = 1))
#' summary(m1)
#'
#' \dontrun{
#' # Forecasting U.S. GDP
#' m1 <- dfm(econ_us,
#'   pre_differenced = "A191RL1Q225SBEA",
#'   keep_posterior = "A191RL1Q225SBEA"
#' )
#'
#' # interpolating low frequency series
#' dta_mixed <- econ_us[, c(1, 3)]
#' predict(dfm(dta_mixed))
#' predict(dfm(dta_mixed, interpolate = TRUE))
#' }
#' @useDynLib bdfm
dfm <- function(data,
                factors = 1,
                lags = "auto",
                forecasts = 0,
                method = c("bayesian", "ml", "pc"),
                scale = TRUE,
                logs = "auto",
                diffs = "auto",
                outlier_threshold = 4,
                frequency_mix = "auto",
                pre_differenced = NULL,
                trans_prior = NULL,
                trans_shrink = 0,
                trans_df = 0,
                obs_prior = NULL,
                obs_shrink = 0,
                obs_df = NULL,
                identification = "pc_long",
                keep_posterior = NULL,
                interpolate = FALSE,
                orthogonal_shocks = FALSE,
                reps = 1000,
                burn = 500,
                verbose = interactive(),
                tol = 0.01
                ) {

  call <- match.call

  method <- match.arg(method) # checks and picks the first if unspecified

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
    all_NA <- apply(X = data, MARGIN = 2, FUN = AllNA)
    if(any(all_NA)){
      stop(paste(colnames(data)[all_NA], collapse=", "), " contain no observations. Remove these series before estimation.")
    }
    ans <- dfm_core(
      Y = data, m = factors, p = lags, FC = forecasts, method = method,
      scale = scale, logs = logs, diffs = diffs, outlier_threshold = outlier_threshold, freq = frequency_mix,
      preD = pre_differenced, Bp = trans_prior, lam_B = trans_shrink, trans_df = trans_df,
      Hp = obs_prior, lam_H = obs_shrink, obs_df = obs_df,
      ID = identification, keep_posterior = keep_posterior, reps = reps,
      burn = burn, verbose = verbose, tol = tol, interpolate = interpolate,
      orthogonal_shocks = orthogonal_shocks
    )
    colnames(ans$values) <- colnames(data)
    ans$dates <- NULL
  } else {
    # no requirement for tsbox if data is ts or mts
    if (inherits(data, "ts")) {
      data_tsp <- attr(data, "tsp")
      data_unclassed <- as.matrix(unclass(data))
    } else {
      # all other time series classes are handled by tsbox
      stopifnot(requireNamespace("tsbox"))
      stopifnot(tsbox::ts_boxable(data))
      data_tsp <- attr(tsbox::ts_ts(data), "tsp")
      data_unclassed <- unclass(tsbox::ts_ts(data))
    }

    attr(data_unclassed, "tsp") <- NULL

    all_NA <- apply(X = data_unclassed, MARGIN = 2, FUN = AllNA)
    if(any(all_NA)){
      stop(paste(colnames(data_unclassed)[all_NA], collapse=", "), " contain no observations. Remove these series before estimation.")
    }

    ans <- dfm_core(
      Y = as.matrix(data_unclassed), m = factors, p = lags, FC = forecasts, method = method,
      scale = scale, logs = logs, diffs = diffs, outlier_threshold = outlier_threshold, freq = frequency_mix,
      preD = pre_differenced, Bp = trans_prior, lam_B = trans_shrink, trans_df = trans_df,
      Hp = obs_prior, lam_H = obs_shrink, obs_df = obs_df,
      ID = identification, keep_posterior = keep_posterior, reps = reps,
      burn = burn, verbose = verbose, tol = tol, interpolate = interpolate,
      orthogonal_shocks = orthogonal_shocks
    )

    # re-apply time series properties and colnames from input
    ans$values <- ts(ans$values, start = data_tsp[1], frequency = data_tsp[3])
    ans$adjusted <- ts(ans$adjusted, start = data_tsp[1], frequency = data_tsp[3])
    ans$factors <- ts(ans$factors, start = data_tsp[1], frequency = data_tsp[3])
    if (!is.null(keep_posterior)) {
      ans$Ymedian <- ts(ans$Ymedian, start = data_tsp[1], frequency = data_tsp[3])
      ans$idx_update <- ts(ans$idx_update, start = data_tsp[1], frequency = data_tsp[3])
    }
    colnames(ans$values) <- colnames(data)

    # put values back into original class (other than ts)
    if (!inherits(data, "ts")) {
      ans$values <- tsbox::copy_class(ans$values, data)
      ans$adjusted <- tsbox::copy_class(ans$adjusted, data)
      ans$factors <- tsbox::copy_class(ans$factors, data, preserve.mode = FALSE)
      if (!is.null(keep_posterior)) {
        ans$Ymedian <- tsbox::copy_class(ans$Ymedian, data)
        ans$idx_update <- tsbox::copy_class(ans$idx_update, data)
      }
    }
  }

  ans$method <- method
  ans$call <- match.call()
  class(ans) <- "dfm"
  return(ans)
}



