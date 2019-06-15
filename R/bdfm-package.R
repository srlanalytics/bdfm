#' Bayesain and Maximum Likelihood Estimation of Dynamic Factor Models
#'
#' Estimates dynamic factor models by simulation using the Durbin and Koopman
#' (2012) disturbance smoother. Maximum likelihood estimation via Watson and
#' Engle (1983) and 2-step estimation via principal components is also
#' supported. Input data may be mixed frequency, noisy, have missing values, or
#' ragged edges with different start or end dates.
#'
#' The best way to start is to have a look at the vignette:
#'
#'   \code{vignette("dfm")}
#'
#' @name bdfm-package
#' @aliases bdfm
#' @docType package
#' @keywords package
#' @seealso \code{\link{dfm}} for the core function and more information on
#'   package usage.
NULL



#' US Economic Data
#'
#' US Economic Data Series from [FRED](https://fred.stlouisfed.org). The latest
#' version of the data can be retrieved by running the example code. You need
#' FRED API Key, which can be optained
#' [here](https://research.stlouisfed.org/docs/api/api_key.html)
#'
#' @source https://fred.stlouisfed.org
#'
#' @docType data
#'
#' @format Time series of class \code{"ts"}.
#'
#' @source FRED
#'
#' @name econ_us
#' @keywords datasets
#' @examples
#'
#' \dontrun{
#' library(tsbox)
#' library(fredr)
#' library(bdfm)
#'
#' # Use Sys.setenv(FRED_API_KEY = 'XXXX') to set FRED API Key.
#' # Optain here: https://research.stlouisfed.org/docs/api/api_key.html
#'
#' series_q <- c(
#'   'A191RL1Q225SBEA',       # 01 Real GDP, seasonally adjusted, quarterly, annualized % change
#'   'W068RCQ027SBEA'         # 02 Governemnt expenditures
#' )
#'
#' series_m <- c(
#'   'USSLIND',               # 03 Federal Reserve leading index, monthly, percent
#'   'PCEDG',                 # 04 persional consumption: durable goods, monthly, level
#'   'PCEND',                 # 05 persional consumption: non-durable goods, monthly, level
#'   'UMCSENT',               # 06 Consumer Sentiment, monthly, delayed 1 month for free data
#'   'UNRATE',                # 07 Unemployment, monthly
#'   'JTSJOL',                # 08 Job openenings, total non-farm
#'   'INDPRO',                # 09 Industrial Production Index, monthly, level
#'   'CSUSHPINSA',            # 10 Case-Shiller home price index, monthly, two month lag, level
#'   'HSN1F',                 # 11 New 1 family houses sold, level
#'   'TSIFRGHT',              # 12 Freight transportation index, monthly, 2-3 month lag, level
#'   'FRGSHPUSM649NCIS',      # 13 CASS freight index, level, not SA
#'   'CAPUTLG2211S',          # 14 Electricity usage, % capacity, monthly
#'   'IPG2211S',              # 15 Electricity, industrial production index, monthly, level
#'   'DGORDER',               # 16 New Orders, durable manufacturing goods, monthly, level
#'   'AMTMNO',                # 17 New Orderes, all manufacuring industries, level
#'   'MNFCTRIRSA',            # 18 Manufacturers inventories:sales ratio
#'   'RETAILIRSA',            # 19 Retail inventories:sales ratio
#'   'WHLSLRIRSA',            # 20 Wholesalers, inventories:sales ratio
#'   'CPILFESL'               # 21 CPI
#' )
#'
#' series_wd <- c(
#'   'ICSA',                  # 22 Initial claims, SA, weekly
#'   'TWEXB',                 # 23 exchange rate index, weekly
#'   'T10Y3M'                 # 24 10Y to 3M treasury spread, daily
#' )
#'
#' multi_fredr <- function(x, by = 0, frequency = NULL) {
#'   z <- lapply(x, fredr, observation_start = as.Date("1980-01-01"), frequency = frequency)
#'   lapply(z, ts_lag, by = by)
#' }
#' econ_us <- ts_ts(do.call(rbind, c(
#'   multi_fredr(series_q, by = "2 month"),
#'   multi_fredr(series_m),
#'   multi_fredr(series_wd, frequency = "m")
#'   )
#' ))
#'
#' # A lot of our data is already seasonally adjusted. We actually only need to SA
#' # three series (though whether it's necessary to adjust consumer sentiment is
#' # debatable).
#' series_sa <- c('UMCSENT', 'CSUSHPINSA', 'FRGSHPUSM649NCIS')
#' econ_us_sa <- do.call(cbind, lapply(
#'   ts_tslist(econ_us[, series_sa]),
#'   function(e) seas_we(e, lags = 3)$values
#' ))
#' window(econ_us[, series_sa], end = end(econ_us_sa)) <- econ_us_sa
#' }
#'
NULL

