# Automatic selection for logs and differences

log_diff <- function(y){
  fq <- get_freq(y)
  rmdr  <- median(which(is.finite(y))%%fq)
  indx  <- seq(1, length(y))%%fq == rmdr
  y <- y[indx]
  ts_reg <- UVreg(x = y[-length(y)] - mean(y, na.rm = TRUE), y = y[-1] - mean(y, na.rm = TRUE), itc = FALSE)
  
  take_log  <- (ts_reg$B - ts_reg$sd)>1 && !any(y<0, na.rm = TRUE)
  take_diff <- !(ts_reg$B + 2*ts_reg$sd)<1
  
  return(c(take_log, take_diff))
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


