#' Seasonal adujustment, usign Watson and Engle (1983) Algorithm
#'
#' to produce maximum likelihood estimates of model parameters
#'
#' @param data data
#' @param lags integer, number of lags in transition equation
#' @param ... further arguments, passed to [seas_factors()]
#' @param ahead number of periods ahead to produce adjustments and forecasts
#' @param transformation string, transform data before calculating seasonality ("auto", "log", or "none")
#' @param span numeric, span for loess() reg used to detrend non-stationary series. NULL implies do not detrend.
#' @param tol numeric, tolerance for convergence of likelihood function (times 100)
#' @param verbose logical, whether to output convergence of iterations
#' @return an object of class `"seas_we"`. Use `predict` to extract the seasonally adjusted values.
#' @export
#' @useDynLib bdfm
#' @examples
#' m <- seas_we(mdeaths)
#' adj_we <- predict(m)
#' library(seasonal)
#' adj_x13 <- predict(seas(mdeaths))
#' tsbox::ts_plot(mdeaths, adj_x13, adj_we)
seas_we <- function(data, lags = 1, ..., ahead = 0, transformation = "auto", span = 0.5, tol = 0.01, verbose = FALSE) {

  if (!requireNamespace("tsbox") & !is.ts(data)) {
    stop('"tsbox" is needed to support non ts-time-series. To install: \n\n  install.packages("tsbox")', call. = FALSE)
  }

  data.orig <- data
  
  if (requireNamespace("tsbox")) {
    data <- tsbox::ts_ts(data)
    dates <- tsbox::ts_df(data)$time #this drops NA values
    if(ahead>0){
      days <- unique(diff(dates))
      if(length(days)==1){
        dates <- c(dates, seq.Date(from = tail(dates,1), by = days, length.out = ahead+1)[-1])
      }else if(all(c(1,3)%in%days)){
        extra_dates <- seq.Date(from = tail(dates,1), by = "day", length.out = ahead+1)[-1]
        extra_dates <- extra_dates[!weekdays(extra_dates)%in%c('Saturday', 'Sunday')]
        dates <- c(dates,extra_dates)
      }else if(any(c(28:31)%in%days)){
        dates <- c(dates, seq.Date(from = tail(dates,1), by = "month", length.out = ahead+1)[-1])
      }else if(any(c(64:66%in%days))){
        dates <- c(dates, seq.Date(from = tail(dates,1), by = "quarter", length.out = ahead+1)[-1])
      }else if(any(c(356,366)%in%days)){
        dates <- c(dates, seq.Date(from = tail(dates,1), by = "year", length.out = ahead+1)[-1])
      }else{
        stop("frequency not supported for forecasts")
      }
    }
  }

  if (NCOL(data) != 1) stop("applicable to single time series only")

  # dates extraction without tsbox (month and quarterly only)
  if (!requireNamespace("tsbox")) {
    fr <- frequency(data)
    
    if (!(fr %in% c(4, 12))) {
      stop('"tsbox" is needed to support non-standard frequencies. To install: \n\n  install.packages("tsbox")', call. = FALSE)
    }
    if (fr == 12) {
      dates <- seq(
        as.Date(paste(start(data)[1], start(data)[2], 1, sep = "-")),
        length.out = length(data) + ahead,
        by = "month"
      )
    }
    if (fr == 4) {
      dates <- seq(
        as.Date(paste(start(data)[1], ((start(data)[2] - 1) * 3) + 1, 1, sep = "-")),
        length.out = length(data) + ahead,
        by = "quarter")
    }
  }

  N <- seas_factors(dates, ...)
  #N <- seas_factors(dates, effect = "weekdays", holiday = NULL)

  scale.factor <- 1000
  
  # shorthand notation
  y <- c(data)
  p <- lags
  
  if(!transformation%in%c("auto", "log", "none")){
    warning(paste("Transformation must be one of auto, log, or none.", "Transformation", transformation,  "not supported, defaulting to auto"))
    transformation <- "auto"
  }
  
  if(transformation == "auto"){
    tmp <- diff(y, differences = 1)
    tmp_var  <- var(tmp, na.rm = TRUE)/length(y)
    tmp_mean <- mean(tmp, na.rm = TRUE)
    #T/F - take logs if mean diff is significantly different (one s.d.) from zero AND
    #      fewer than 5% of observations are less than zero (i.e. due to errors in data)
    # Properly we should use the second diff, but the financial crises messes that up.
    take_log <- tmp_mean/sqrt(tmp_var)>1 && sum(y<0)/length(y) < 0.05
    if(take_log){
      y[y<0] <- NA
      y      <- log(y)
    }
  }else if(transformation == "log"){
    take_log <- TRUE
    y[y<0] <- NA
    y      <- log(y)
  }else if(transformation == "none"){
    take_log <- FALSE
  }

  # Note that data is scaled up. This avoids small matrix determinents in the
  # calculations which can make results inaccurate.
  y <- scale.factor * y
  
  if(is.numeric(span)){
    dtrend    <- loess(y ~ seq(1,length(y)), span = span, na.action = na.exclude)
    y         <- y-predict(dtrend)
  }else if(is.null(span)){
    y <- y - mean(y, na.rm = TRUE)
  }else{
    warning("span should be either numeric or NULL. Non-numeric value entered; 
            data will not be detrended.")
    y <- y - mean(y, na.rm = TRUE)
  }
  
  row_N <- nrow(N) # number of observations (rows)
  if(row_N>length(y)){
    y <- c(y, rep(NA, row_N-length(y)))
  }
  y <- as.matrix(y) # required for C++ stuff
  B <- matrix(0, 1, p) # empty parameter matrix to be filled in (transition equation)
  B[1,1] <- .1
  q <- var(y, na.rm = T) / 2 # arbitrary initial guess for q
  r <- q # arbitrary initial guess for r
  M <- t(QuickReg(N, y)) # abitrary initial guess for M (OLS using the C++ function QuickReg)

  Lik0 <- -1e10
  count <- 0
  conv <- 100
  while (conv > tol | count < 6) { # loop until the likelihood function converges
    est <- KSeas(B, q, M, r, y, N) # calculate estimates in C++
    B <- est$B
    q <- est$q
    r <- est$r
    M <- est$M
    Z0 <- est$Z0
    Lik1 <- est$Lik
    conv <- 100 * (Lik1 - Lik0) / abs(Lik1 + Lik0)
    Lik0 <- Lik1
    if (verbose) {
      message(conv)
    }
    count <- count + 1
  }

  if(take_log){
    sa   <- exp(N %*% t(M) / scale.factor)
    y_sa <- c(data.orig)/sa[1:length(data.orig)]
  }else{
    sa   <- N %*% t(M) / scale.factor
    y_sa <- c(data.orig) - sa[1:length(data.orig)]
  }
  
  # make values a ts timeseries
  tsp <- tsp(data)
  y_sa <- ts(y_sa, start = tsp[1], frequency = tsp[3])
  sa <- ts(sa, start = tsp[1], frequency = tsp[3])

  # put values back into original class
  if (!inherits(data.orig, "ts")) {
    y_sa  <- tsbox::copy_class(y_sa, data.orig)
    sa <- tsbox::copy_class(sa, data.orig)
  }

  z <- list(
    values = y_sa, # seasonally adjusted Y
    factor = sa,  # seasonal adjustments
    M = M,
    take_log = take_log #T/F, log transformation?
  )

  class(z) <- "seas_we"
  z

}



# methods
#' @export
#' @method predict seas_we
predict.seas_we <- function(object, ...) {
  object$values
}

# methods
#' @export
#' @method print seas_we
print.seas_we <- function(x, ...) {
  cat("A seas_we object....\n")
}



