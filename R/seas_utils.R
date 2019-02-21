# applies a function on each continous segment of 1s.
# @example
#
# dist_unif <- function(x) rep(1 / length(x), length(x))
# dist_norm <- function(x) {
#   z <- dnorm(seq(-2, 2, length.out = length(x)))
#   z / sum(z)
# }
# x <- rep(c(1, 0, 0, 1, 0, 0), each = 3)
# apply_to_segments(x, dist_norm)
apply_to_segments <- function(x, fun) {
  stopifnot(identical(sort(unique(x)), c(0, 1)))

  dx <- diff(x)
  idx.start <- which(dx == 1) + 1
  # add 1 at start if vector starts with 1
  if (x[1] == 1) {
    idx.start <- c(1, idx.start)
  }

  idx.end <- which(dx == -1)
  # add length(x) at end if vector ends with 1
  if (x[length(x)] == 1) {
    idx.end <- c(idx.end, length(x))
  }

  stopifnot(identical(length(idx.start), length(idx.end)))

  z <- x
  for (i in seq(idx.start)) {
    z[idx.start[i]:idx.end[i]] <- fun(z[idx.start[i]:idx.end[i]])
  }
  z
}


daily_seq <- function(dates) {
  freq <- which_freq(dates)

  if (freq == "day") return(dates)

  if (freq == "month") {
    start_of_period <- function(x) {
      as.Date(paste(
        as.POSIXlt(x)$year + 1900L,
        as.POSIXlt(x)$mon + 1L,
        1,
        sep = "-"
      ))
    }
    end_of_period <- function(x) {
      seq(start_of_period(x), length.out = 2, by = "month")[2] - 1
    }
  }

  if (freq == "quarter") {
    start_of_period <- function(x) {
      as.Date(paste(
        as.POSIXlt(x)$year + 1900L,
        floor(as.POSIXlt(x)$mon / 3) + 1L,
        1,
        sep = "-"
      ))
    }
    end_of_period <- function(x) {
      seq(start_of_period(x), length.out = 2, by = "quarter")[2] - 1
    }
  }

  seq.Date(from = start_of_period(dates[1]), to = end_of_period(dates[length(dates)]), by = "day")
}


year_month <- function(dates) {
  dates.POSIXlt <- as.POSIXlt(dates)
  paste(
    dates.POSIXlt$year + 1900L,
    sprintf("%02d", dates.POSIXlt$mon + 1L),
    sep = "-"
  )
}

year_quarter <- function(dates) {
  dates.POSIXlt <- as.POSIXlt(dates)
  paste(
    dates.POSIXlt$year + 1900L,
    sprintf("%02d", floor(dates.POSIXlt$mon / 3) + 1L),
    sep = "-"
  )
}

which_freq <- function(dates) {
  stopifnot(inherits(dates, "Date"))
  udiff <- unique(diff(dates))
  if (identical(udiff, 1)) {
    freq <- "day"
  } else if (all(udiff %in% 28:31)) {
    freq <- "month"
  } else if (all(udiff %in% 90:92)) {
    freq <- "quarter"
  } else {
    stop("unknown frequency. 'dates' must be daily, monthly or quarterly", call. = FALSE)
  }
  freq
}

# relevant timestamp function, depending on frequency: day, month, or quarter
timestamp_fun <- function(freq) {
  if (freq == "day") timestamp_fun <- function(x) x
  if (freq == "month") timestamp_fun <- year_month
  if (freq == "quarter") timestamp_fun <- year_quarter
  timestamp_fun
}


# Working Days per Month, Quarter
#
# Could be stored permanently as they do not depend on data
#
# @examples
# from <- as.Date("2000-01-01")
# to <- as.Date("2005-03-01")
# dates.d <- seq(from, to, by = "day")
# dates.m <- seq(from, to, by = "month")
# dates.q <- seq(from, to, by = "quarter")
# workdays(dates.d)
# workdays(dates.m)
# workdays(dates.q)
workdays <- function(dates) {
  freq <- which_freq(dates)
  timestamp <- timestamp_fun(freq)

  sq <- daily_seq(dates)
  is.weekday <- as.integer(!(weekdays(sq) %in% c("Saturday", "Sunday")))

  z <- as.matrix(tapply(is.weekday, timestamp(sq), sum))
  colnames(z) <- "workdays"
  z
}


dummy_matrix <- function(x) {
  # factors with same level order as x
  xf <- factor(x, levels = unique(x))
  z <- model.matrix(~ xf + 0)
  attr(z, "assign") <- NULL
  attr(z, "contrasts") <- NULL
  colnames(z) <- levels(xf)
  # z[,-1]
  z
}



# Date based version of seasonal::genhol
#
# that uses a sequence of dates, rather than tsp info
#
# should gain same arguments as 'genhol' in seasonal, start, end
genhol_dates <- function(x, dates, start = 0, end = 0) {
  stopifnot(inherits(x, "Date"))
  stopifnot(inherits(dates, "Date"))
  freq <- which_freq(dates)
  timestamp <- timestamp_fun(freq)
  sq <- daily_seq(dates)


  idx.holiday <- which(sq %in% x)
  st <- idx.holiday + start
  en <- idx.holiday + end
  idx.holiday.extended <- unlist(Map(function(st, en) st:en, st, en))

  is.holiday <- logical(length(sq))
  is.holiday[idx.holiday.extended] <- TRUE

  z <- as.matrix(tapply(is.holiday, timestamp(sq), sum))
  colnames(z) <- "holiday"
  z
}
