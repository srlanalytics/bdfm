AnyNA <- function(y) {
  any(!is.finite(y))
}

AllNA <- function(y) {
  all(!is.finite(y))
}

NumNA <- function(y) {
  length(which(is.na(y)))
}

# convert string to numeric index value
standardize_index <- function(x, Y) {
  argname <- deparse(substitute(x))
  if (is.character(x)) {
    out <- unlist(sapply(x, FUN = grep, colnames(Y)))
  } else if (is.numeric(x)) {
    out <- x
  } else if (is.logical(x)){
    if(length(x)!=NCOL(Y)){
      stop("If argument '", argname,
           "' is logical, its length must equal the number of columns in the input data",
           call. = FALSE
      )
    }
    out <- which(x)
  } else {
    stop("Argument '", argname,
      "' must be either a character (string) vector, a logical vector, or numeric index values",
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

#return only values that correspond to the end of the quarter (or whatever the low frequency periods are)
drop_intermediates <- function(ind, freq, Y_raw, vals){
  y_raw <- Y_raw[,ind]
  y     <- vals[,ind]
  rmdr  <- median(which(is.finite(y_raw))%%freq[ind])
  indx    <- seq(1, length(y_raw))%%freq[ind] == rmdr
  y[!indx]<- NA
  return(y)
}


