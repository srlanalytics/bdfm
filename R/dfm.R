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
#'   This can be entered as a `"ts"` object or as a matrix. If
#'   [tsbox](https://tsbox.help) is installed, any ts-boxable time series can be
#'   supplied (`ts`, `xts`, `zoo`, `data.frame`, `data.table`, `tbl`, `tbl_ts`,
#'   `tbl_time`, or `timeSeries`)
#' @param factors integer. The number of unobserved factors to be estimated. A
#'   larger number of factors leads to a more complex model. Denoted as 'm'
#'   in the documentation.
#' @param lags integer. The number of lags in the transition equation. If
#'  `"auto"` (default), the number is equal to highest frequency in `data`.
#'  Denoted as 'p' in the documentation.
#' @param forecasts integer. Number of periods ahead to forecasts.
#' @param method character. Method to be used; one of `"bayesian"`, `"ml"` or
#'   `"pc"`. See details.
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
#'   differences (If series are specified in `diffs`, this is not needed.)
#' @param trans_prior m x mp (m: factors, p: lags) prior matrix for B (the transition matrix) in the transition equation. Default is
#'   zeros. E.g., to use a random walk prior with m factors and p lags, set
#'   `trans_prior = cbind(diag(1,m,m), matrix(0,m,m*(p-1)))`.
#' @param trans_shrink numeric. Prior tightness on B matrix in transition equation where a value of zero is
#'   used to denote an improper (flat) prior (i.e. no shrinkage). Use
#'   to shrink forecast values towards the prior `trans_prior`,
#'   which may help reduce parameter uncertainty in estimation.
#' @param trans_df numeric. Prior degrees of freedom for inverse-Wishart distribution of shocks
#'   in the transition equation, where 0 implies no shrinkage.
#'   Shrinking shocks to the transition equation will increase the magnitude
#'   of shocks to the observation equation dampening updates from observed
#'   series. High values of `trans_df` can lead to
#'   instability in simulations.
#' @param obs_prior k x m (k: observed series, m: factors) prior matrix for H (loadings) in the observation equation
#'   Default is zeros.
#' @param obs_shrink numeric. Prior tightness on H (loadings) in the observation
#'   equation where a value of zero is
#'   used to denote an improper (flat) prior (i.e. no shrinkage). A greater value will shrink estimates of loadings more
#'   aggressively towards the prior `obs_prior`. When the prior is zero (the
#'   default value), this is an alternative (and typically more stable) approach
#'   to dampening the impact of updates from observed series.
#' @param obs_df named vector (see details). prior degrees of freedom
#'   for inverse chi-squared distribution in the observation equation. This is useful to give
#'   specific series a larger weight, e.g. 1. (default 0).
#' @param identification names or index values (see details), or character.
#'   Factor identification. `"pc_long"` (default) identifies on principal components
#'   from series with at least the median number of observations. `"pc_wide"`
#'   identifies on principal components using all series, where rows of the observations
#'   matrix containing missing data are omitted. `"name"` uses Stock and Watson's "naming
#'   factors" identification, i.e. identifying on the first m series provided where m is the
#'   number of factors. Identification can also be done manually, by supplying names or index
#'   values from which identifying series are derived via principal components.
#' @param keep_posterior names or index values (see details). Series of which to
#'   keep the full posterior distribution of predicted values (method
#'   `"bayesian"` only). This is useful for forecasting as the posterior median forecast
#'   value tends to me more accurate than forecasts using the posterior median parameter
#'   estimates, and allows for the evaluation of forecast accuracy.
#' @param interpolate logical. Should output return intra-frequency estimates of low
#'   frequency observables? Put differently, if the model includes monthly and quarterly data,
#'   should output include
#'   estimates of quarterly data every month (where quarterly refers to an aggregate of the
#'   current and previous 2 months; for `interpolate = TRUE`) or just at the end of the quarter (months 3, 6, 9, and 12;
#'   for `interpolate = FALSE`, default)?
#' @param orthogonal_shocks logical. Return a rotation of the model with orthogonal
#'   shocks and factors. This is used to isolate the impact of each factor on observables,
#'   allowing for a clean interpretation of how shocks (which, if `TRUE` are not correlated)
#'   impact observed series.
#' @param reps integer. Number of repetitions for MCMC sampling
#' @param burn integer. Number of iterations to burn in MCMC sampling
#' @param verbose logical. Print status of function during evaluation. Default is
#'  `TRUE` in interactive mode, `FALSE` otherwise, so it does not appear, e.g.,
#'  in `reprex::reprex()`.
#' @param tol numeric. Tolerance for convergence of EM algorithm (method `"ml"`
#'   only). The default value is 0.01 which corresponds to the convergence
#'   criteria used in Doz, Giannone, and Reichlin (2012).
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
                verbose = interactive() && !isTRUE(getOption("knitr.in.progress")),
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
      data_ts <- data
    } else {
      # all other time series classes are handled by tsbox
      stopifnot(requireNamespace("tsbox"))
      stopifnot(tsbox::ts_boxable(data))
      data_ts <- tsbox::ts_ts(data)
    }

    data_tsp <- attr(data_ts, "tsp")
    data_unclassed <- as.matrix(unclass(data_ts))
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
    ans$factor_update <- lapply(ans$factor_update, function(e) ts(e, start = data_tsp[1], frequency = data_tsp[3]))
    if (!is.null(keep_posterior)) {
      ans$Ymedian <- ts(ans$Ymedian, start = data_tsp[1], frequency = data_tsp[3])
      ans$Ystore <- ts(ans$Ystore, start = data_tsp[1], frequency = data_tsp[3])
      ans$idx_update <- ts(ans$idx_update, start = data_tsp[1], frequency = data_tsp[3])
    }
    colnames(ans$values) <- colnames(data_ts)

    # put values back into original class (other than ts)
    if (!inherits(data, "ts")) {
      ans$values <- tsbox::copy_class(ans$values, data)
      ans$adjusted <- tsbox::copy_class(ans$adjusted, data)
      ans$factors <- tsbox::copy_class(ans$factors, data, preserve.mode = FALSE)
      ans$factor_update <- lapply(ans$factor_update, function(e) tsbox::copy_class(e, data))
      if (!is.null(keep_posterior)) {
        ans$Ymedian <- tsbox::copy_class(ans$Ymedian, data)
        ans$Ystore <- tsbox::copy_class(ans$Ystore, data)
        ans$idx_update <- tsbox::copy_class(ans$idx_update, data)
      }
    }
  }

  ans$method <- method
  ans$call <- match.call()
  class(ans) <- "dfm"
  return(ans)
}



