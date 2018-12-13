#' Seasonal Factors at Daily, Monthly or Quarterly Frequency
#'
#' Generate seasonal factors at daily, monthly or quarterly frequency
#'
#' @param dates vector of class `Date`
#' @param effect character, one of "months", "quarters", "workdays", "weekdays"
#' @param holiday vector of class `Date`, holidays to adjust for, e.g. `easter
#' @param start integer, shift start of holiday
#' @param end integer, shift end of holiday
#' @param dist_fun distribution function for `effect` `"months"` and `"quarters"`. E.g., `dist_norm`.
#' @export
#' @examples
#' from <- as.Date("2000-01-01")
#' to <- as.Date("2005-03-01")
#' dates.d <- seq(from, to, by = "day")
#' dates.m <- seq(from, to, by = "month")
#' dates.q <- seq(from, to, by = "quarter")
#' seas_factors(dates.d)
#' seas_factors(dates.m)
#' seas_factors(dates.q)
#' seas_factors(dates.d, dist_fun = dist_norm)
seas_factors <- function(dates,
                         effect = c("months", "workdays"),
                         holiday = easter,
                         start = -1, end = 0,
                         dist_fun = NULL
                         ) {
  effect <- match.arg(
    effect,
    choices = c("months", "quarters", "workdays", "weekdays"),
    several.ok = TRUE
  )

  if ("months" %in% effect && "quarter" %in% effect) {
    stop("cannot use 'months' and 'quarter' simultanously", call. = FALSE)
  }
  if ("weekdays" %in% effect && "workdays" %in% effect) {
    stop("cannot use 'weekdays' and 'workdays' simultanously", call. = FALSE)
  }

  stopifnot(inherits(dates, "Date"))

  # z <- as.matrix(rep(1, length(dates)))
  # colnames(z) <- "(intercept)"

  z <- NULL

  # all freqs
  if ("months" %in% effect) {
    z <- cbind(z, dummy_matrix(tolower(months(dates, abbreviate = TRUE))))
  }
  if ("quarters" %in% effect) {
    z <- cbind(z, dummy_matrix(quarters(dates)))
  }
  # apply distribution function, if provided
  if (!is.null(dist_fun)) {
    stopifnot(is.function(dist_fun))
    z <- apply(z, 2, apply_to_segments, fun = dist_fun)
  }

  if ("workdays" %in% effect) {
    z <- cbind(z, workdays(dates))
  }
  if ("weekdays" %in% effect) {
    z <- cbind(z, dummy_matrix(tolower(weekdays(dates, abbreviate = TRUE))))
  }

  if (!is.null(holiday)) {
    z <- cbind(z, genhol_dates(holiday, dates, start = start, end = end))
  }

  z

}


#' @export
#' @param x numeric vector to be distributed
#' @name seas_factors
dist_unif <- function(x) rep(1 / length(x), length(x))

#' @export
#' @name seas_factors
dist_norm <- function(x) {
  z <- dnorm(seq(-2, 2, length.out = length(x)))
  z / sum(z)
}


