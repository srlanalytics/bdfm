#' Seasonal adujustment, usign Watson and Engle (1983) Algorithm
#'
#' to produce maximum likelihood estimates of model parameters
#'
#' @param data data
#' @param lags integer, number of lags in transition equation
#' @param ... further arguments, passed to [seas_factors()]
#' @param tol numeric, tolerance for convergence of likelihood function (times 100)
#' @param verbose logical, whether to output convergence of iterations
#' @return an object of class `"seas_we"`. Use `predict` to extract the seasonally adjusted values.
#' @export
#' @useDynLib BDFM
#' @examples
#' m <- seas_we(mdeaths)
#' adj_we <- predict(m)
#' library(seasonal)
#' adj_x13 <- predict(seas(mdeaths))
#' tsbox::ts_plot(mdeaths, adj_x13, adj_we)
seas_we <- function(data, lags = 1, ..., tol = 0.01, verbose = FALSE) {

  if (!requireNamespace("tsbox") & !is.ts(data)) {
    stop('"tsbox" is needed to support non ts-time-series. To install: \n\n  install.packages("tsbox")', call. = FALSE)
  }

  data.orig <- data

  if (requireNamespace("tsbox")) {
    data <- tsbox::ts_ts(data)
    dates <- tsbox::ts_df(data)$time
  }

  if (NCOL(data) != 1) stop("applicable to single time series only")

  # dates extraction without tsbox (month and quarterly only)
  if (!requireNamespace("tsbox")) {
    fr <- frequency(data)
    if (NCOL(data) != 1) stop("applicable to single time series only")

    if (!(fr %in% c(4, 12))) {
      stop('"tsbox" is needed to support non-standard frequencies. To install: \n\n  install.packages("tsbox")', call. = FALSE)
    }
    if (fr == 12) {
      dates <- seq(
        as.Date(paste(start(data)[1], start(data)[2], 1, sep = "-")),
        length.out = length(data),
        by = "month"
      )
    }
    if (fr == 4) {
      dates <- seq(
        as.Date(paste(start(data)[1], ((start(data)[2] - 1) * 3) + 1, 1, sep = "-")),
        length.out = length(data),
        by = "quarter")
    }
  }


  N <- seas_factors(dates, ...)

  # Note that data is scaled up by 100. This avoids small matrix determinents in the
  # calculations which can make results inaccurate.

  scale.factor <- 100

  # shorthand notation
  y <- scale.factor * c(data)
  p <- lags

  row_N <- nrow(N) # number of observations (rows)
  y <- as.matrix(y) # required for C++ stuff
  if(row_N>nrow(y)){
    y <- c(y, rep(NA, row_N-nrow(y)))
  }
  m <- ncol(N) # number of seasonal factors (columns of N)
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

  y_sa <- (y - N %*% t(M)) / scale.factor
  sa <- N %*% t(M) / scale.factor


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
    M = M
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



