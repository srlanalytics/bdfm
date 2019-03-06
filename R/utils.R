AnyNA <- function(y) {
  any(is.na(y))
}

NumNA <- function(y) {
  length(which(is.na(y)))
}

Y_sub <- function(Y) {
  # lots of ways to select a full submatrix of the data... this is one idea
  Ysub <- na.omit(Y)
  obs_old <- nrow(Ysub) * ncol(Ysub)
  miss_rows <- apply(Y, MARGIN = 2, FUN = NumNA)
  dims <- dim(na.omit(Y[, -which(miss_rows == max(miss_rows))]))
  obs_new <- dims[1] * dims[2]
  if (obs_new < obs_old) {
    ind <- apply(Y, MARGIN = 1, FUN = AnyNA)
    Ysub <- Y[!ind, ]
  }
  while (obs_new > obs_old || obs_old == 0) {
    obs_old <- obs_new
    miss_rows <- apply(Y, MARGIN = 2, FUN = NumNA)
    Y <- Y[, -which(miss_rows == max(miss_rows))]
    ind <- apply(Y, MARGIN = 1, FUN = AnyNA) # indexes of Y that do not correspond to Ysub
    Ysub <- Y[!ind, ]
    # Ysub      <- na.omit(Y)
    obs_new <- nrow(Ysub) * ncol(Ysub)
    # print(dim(Ysub))
  }
  return(list(Ysub = Ysub, ind = !ind))
}

# #' scale
# #'
# #' Scale data for factor model estimation and store scaling values
# #'
# #' @param x a numeric matrix(like object)
# #' @param center either a logical value or a numeric vector of length equal to the number of columns of x
# #' @param scale either a logical value or a numeric vector of length equal to the number of columns of x
# #' @export
# #' @useDynLib bdfm
# bdfm_scale <- function(x, center = TRUE, scale = TRUE) {
#   x <- 100 * (scale(x, center = center, scale = scale))
#   options(scaled_scale = attr(x, "scaled:scale"))
#   options(scaled_center = attr(x, "scaled:center"))
#   return(x)
# }
#
# #' unscale
# #'
# #' Unscale data using values from bdfm::scale
# #'
# #' @param x a numeric matrix(like object)
# #' @param idx An index of which elements of scaled_center and scaled_scale to use
# #' @export
# #' @useDynLib bdfm
# unscale <- function(x, idx = NULL) {
#   scaled_center <- getOption("scaled_center")
#   scaled_scale <- getOption("scaled_scale")
#   if (is.null(idx)) {
#     idx <- 1:ncol(x)
#   }
#   if (is.null(dim(x))) {
#     if (length(idx) != 1) {
#       stop("Index values idx must have the same length as the number of columns of x")
#     }
#     if (!is.null(scaled_scale)) {
#       x <- scaled_scale[idx] * x / 100
#     }
#     if (!is.null(scaled_center)) {
#       x <- x + scaled_center[idx]
#     }
#   } else {
#     if (length(idx) != ncol(x)) {
#       stop("Index values idx must have the same length as the number of columns of x")
#     }
#     if (!is.null(scaled_scale)) {
#       x <- (matrix(1, nrow(x), 1) %x% t(scaled_scale[idx])) * (x / 100)
#     }
#     if (!is.null(scaled_center)) {
#       x <- x + (matrix(1, nrow(x), 1) %x% t(scaled_center[idx]))
#     }
#   }
#   return(x)
# }

# convert string to numeric index value
standardize_index <- function(x, Y) {
  argname <- deparse(substitute(x))
  if (is.character(x)) {
    out <- unlist(sapply(x, FUN = grep, colnames(Y)))
  } else if (is.numeric(x)) {
    out <- x
  } else {
    stop("Argument '", argname,
      "' must be either a character (string) vector or numeric index values",
      call. = FALSE
    )
  }
  out
}

standardize_numeric <- function(x, Y) {
  argname <- deparse(substitute(x))
  if (!is.null(names(x))) {
    if (any(!(names(x) %in% colnames(Y)))) {
      stop("Names of '", argname, "' must correspond to colnames of 'data'", call. = FALSE)
    }
    out <- setNames(numeric(NCOL(Y)), colnames(Y))
    out[names(x)] <- x
  } else {
    out <- x
  }
  if (length(out) != NCOL(Y)) {
    stop("Lenght of '", argname, "' must correspond to number of columns of 'data'", call. = FALSE)
  }
  out
}

#util to get gain for current period only
one_through_m <- function(x,m){
  return(t(x[1:m,,drop = FALSE]))
}

#hadamard function for mapply to get forecast updates
hadamard <- function(a,b){
  return(a*b)
}

#matrix multiplication for lapply to get forecast updates
multiply_each <- function(x, loading, y_scale = 1){
  return(y_scale*(x%*%loading))
}

factor_updates <- function(idx, gain, PE){ #get factor updates from gain and prediction error
  g  <- gain[[idx]]
  pe <- PE[[idx]]
  fc_ud <- g*(matrix(1,1,NCOL(g))%x%pe)
  return(fc_ud)
}


name_updates <- function(idx, Y, fc_update){
  y <- Y[idx,]
  fc_ud <- fc_update[[idx]]
  if(any(is.finite(y))){
    row.names(fc_ud) <- names(y)[is.finite(y)]
  }
  return(fc_ud)
}

scale_updates <- function(idx, fc_update, y_scale){
  y <- Y[idx,]
  fc_ud <- fc_update[[idx]]
  if(any(is.finite(y))){
    ind <- is.finite(y)
    fc_ud <- y_scale[ind]*fc_ud/100
  }
  return(fc_ud)
}

# auto detect frequency
get_freq <- function(y) {
  out <- median(diff(which(!is.na(y))))
}

# diff mixed frequency data keeping the same length for observed series
mf_diff <- function(ind, fq, Y) {
  Y   <- as.matrix(Y)
  out <- c(rep(NA, fq[ind]), diff(Y[, ind], lag = fq[ind]))
}

# x <- c(NA, 2, 4, NA, NA, 2, NA)
# na_appox(x)
na_appox <- function(x) {
  idx <- seq_along(x)

  x_here <- x[!is.na(x)]
  idx_here <- idx[!is.na(x)]

  # rule 2: NA outside of interpolation should be filled with last known value
  ap <- approx(idx_here, x_here, idx, rule = 2)

  ap$y
}

# convert differenced data back to levels
level <- function(ind, fq, Y_lev, vals) {
  y_lev <- Y_lev[,ind]
  #identify which periods should have observations
  #for high frequency/uniform frequency observations this will be every period
  rmdr    <- median(which(is.finite(y_lev))%%fq[ind])
  indx    <- seq(1, length(y_lev))%%fq[ind] == rmdr
  y       <- vals[indx, ind]
  cs      <- cumsum(y)
  appox   <- na_appox(y_lev[indx] - cs)
  y_lev[indx]   <- appox + cs #return at same frequency as input
  return(y_lev)
}

# convert differenced data back to levels
level_simple <- function(val, y_lev, fq) {
  #identify which periods should have observations
  #for high frequency/uniform frequency observations this will be every period
  rmdr    <- median(which(is.finite(y_lev))%%fq)
  indx    <- seq(1, length(y_lev))%%fq == rmdr
  y       <- val[indx]
  cs      <- cumsum(y)
  appox   <- na_appox(y_lev[indx] - cs)
  y_lev[indx]   <- appox + cs #return at same frequency as input
  return(y_lev)
}


# act <- c(1, 1, 2, 3, NA, NA, NA, 6, 7, 8, NA)
# dd <- c(NA, 1, 2, 1, 1, 1, 5, 2, 1, 2, 3)

# spline(act, n = length(act))

# cs <- c(0, cumsum(dd[-1]))

# offset <- act-cs


# offset <- c(1, 0, -1, -1, -2.5, -4, -5.5, -7, -7, -8, -8)


# cs+offset


# approx(offset)





# Dates_to_cpp <- function(dates) {
#   Dates_cpp <- matrix(0, length(dates), 3)
#   Dates_cpp[, 1] <- as.numeric(format(dates, "%Y"))
#   Dates_cpp[, 2] <- as.numeric(format(dates, "%m"))
#   Dates_cpp[, 3] <- as.numeric(format(dates, "%d"))
#   colnames(Dates_cpp) <- c("Year", "Month", "day")
#   return(Dates_cpp)
# }

# Dates_to_R <- function(Dates_cpp) {
#   dts <- paste0(as.character(Dates_cpp[, 1]), "-", as.character(Dates_cpp[, 2]), "-", as.character(Dates_cpp[, 3]))
#   dts <- as.Date(dts)
#   return(dts)
# }

# #' End Next Month Date
# #'
# #' Get date for the end of the next month
# #'
# #' @param last_date date in the month of interest
# #' @export
# End_Next_Month <- function(last_date) {
#   r <- length(last_date)
#   ymd <- Dates_to_cpp(last_date) + matrix(1, r, 1) %x% matrix(c(0, 1, 0), 1, 3)
#   if (any(ymd[, 2] == 13)) {
#     indx <- ymd[, 2] == 13
#     ymd[indx, 1] <- ymd[indx, 1] + 1
#     ymd[indx, 2] <- 1
#   }
#   end_next_month <- Dates_to_R(end_of_month(ymd))
#   return(end_next_month)
# }

# #' End of Current Month Date
# #'
# #' Get date for the end of the current month
# #'
# #' @param last_date date in the month of interest
# #' @export
# #' @useDynLib bdfm
# End_This_Month <- function(last_date) {
#   end_this_month <- Dates_to_R(end_of_month(Dates_to_cpp(last_date)))
#   return(end_this_month)
# }

# #' End of Last Month Date
# #'
# #' Get date for the end of the previous month
# #'
# #' @param last_date date in the month of interest
# #' @export
# #' @useDynLib bdfm
# End_Last_Month <- function(last_date) {
#   r <- length(last_date)
#   ymd <- Dates_to_cpp(last_date) - matrix(1, r, 1) %x% matrix(c(0, 1, 0), 1, 3)
#   if (any(ymd[, 2] == 0)) {
#     indx <- ymd[, 2] == 0
#     ymd[indx, 1] <- ymd[indx, 1] - 1
#     ymd[indx, 2] <- 12
#   }
#   end_last_month <- Dates_to_R(end_of_month(ymd))
#   return(end_last_month)
# }

# #' End of Month Date
# #'
# #' Get date for the last day of month shift months ago (-) or ahead (+)
# #'
# #' @param last_date date in the month of interest
# #' @param shift Months ahead (positive value) or behind (negative value)
# #' @export
# #' @useDynLib bdfm
# End_of_Month <- function(last_date, shift = 0) {
#   r <- length(last_date)
#   ymd <- Dates_to_cpp(last_date) + matrix(1, r, 1) %x% matrix(c(0, shift, 0), 1, 3)
#   if (any(ymd[, 2] > 12)) {
#     indx <- ymd[, 2] > 12
#     exx <- ymd[indx, 2]
#     yrs <- floor(exx / 12)
#     mths <- exx - 12 * yrs
#     ymd[indx, 1] <- ymd[indx, 1] + yrs
#     ymd[indx, 2] <- mths
#   }
#   if (any(ymd[, 2] < 1)) {
#     indx <- ymd[, 2] < 1
#     exx <- -ymd[indx, 2]
#     yrs <- floor(exx / 12)
#     mths <- 12 - (exx - 12 * yrs)
#     ymd[indx, 1] <- ymd[indx, 1] - yrs - 1
#     ymd[indx, 2] <- mths
#   }

#   end_next_month <- Dates_to_R(end_of_month(ymd))
#   return(end_next_month)
# }
