# Automatic selection for logs and differences

log_diff <- function(y, fast = FALSE){
  if(!fast){
    fq <- get_freq(y)
    rmdr  <- median(which(is.finite(y))%%fq)
    indx  <- seq(1, length(y))%%fq == rmdr
    y <- y[indx]
  }
  # 1st test: test AR(1) coefficient
  # rm outlier here reffers to how many iterations of removing outliers should we do. The threshold is set at 5 sd.
  ts_reg <- UVreg(x = y[-length(y)] - mean(y, na.rm = TRUE), y = y[-1] - mean(y, na.rm = TRUE), rm_outlier = 2)
  
  # 2nd test: does the data grow over time? Normally the first test would get this, but often with super noisy
  # daily data it misses. This is a sort of fall-back.
  lin_reg <- UVreg(x = seq(length(y)), y = y - mean(y, na.rm = TRUE), rm_outlier = 2)
  
  take_diff <- !(ts_reg$B + 3*ts_reg$sd)<1 || (abs(lin_reg$B) - 2*lin_reg$sd)>0
  take_log  <- take_diff && !any(y<0, na.rm = TRUE)  #(ts_reg$B - ts_reg$sd)>1 && !any(y<0, na.rm = TRUE)
  
  log_diff <- c(take_log, take_diff)
  names(log_diff) <- c("Take Logs", "Take Diffs")
  
  return(log_diff)
}


#' Automatically determine if time series data should be differenced or log differenced
#'
#' @param X data with time indexed by rows
#' @export
#' @examples
#' should_log_diff(EuStockMarkets)
#' @useDynLib bdfm
should_log_diff <- function(X){
  out <- apply(X = as.matrix(X), MARGIN = 2, FUN = log_diff)
  rownames(out) <- c("log", "diff")
  return(out)
}


