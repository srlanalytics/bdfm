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
  while(obs_new > obs_old || obs_old == 0){
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

Dates_to_cpp <- function(dates){
  Dates_cpp <- matrix(0,length(dates),3)
  Dates_cpp[,1] <- as.numeric(format(dates, "%Y"))
  Dates_cpp[,2] <- as.numeric(format(dates, "%m"))
  Dates_cpp[,3] <- as.numeric(format(dates, "%d"))
  colnames(Dates_cpp) <- c("Year", "Month", "day")
  return(Dates_cpp)
}

Dates_to_R <- function(Dates_cpp){
  dts <- paste0(as.character(Dates_cpp[,1]),"-",as.character(Dates_cpp[,2]),"-",as.character(Dates_cpp[,3]))
  dts <- as.Date(dts)
  return(dts)
}

#' End Next Month Date
#'
#' Get date for the end of the next month
#'
#' @param last_date date in the month of interest
#' @export
End_Next_Month <- function(last_date){
  r   <- length(last_date)
  ymd <- Dates_to_cpp(last_date) + matrix(1,r,1)%x%matrix(c(0,1,0),1,3)
  if(any(ymd[,2] == 13)){
    indx   <- ymd[,2]==13
    ymd[indx,1] = ymd[indx,1]+1
    ymd[indx,2] = 1
  }
  end_next_month <- Dates_to_R(end_of_month(ymd))
  return(end_next_month)
}

#' End of Current Month Date
#'
#' Get date for the end of the current month
#'
#' @param last_date date in the month of interest
#' @export
#' @useDynLib bdfm
End_This_Month <- function(last_date){
  end_this_month <- Dates_to_R(end_of_month(Dates_to_cpp(last_date)))
  return(end_this_month)
}

#' End of Last Month Date
#'
#' Get date for the end of the previous month
#'
#' @param last_date date in the month of interest
#' @export
#' @useDynLib bdfm
End_Last_Month <- function(last_date){
  r   <- length(last_date)
  ymd <- Dates_to_cpp(last_date) - matrix(1,r,1)%x%matrix(c(0,1,0),1,3)
  if(any(ymd[,2] == 0)){
    indx   <- ymd[,2]==0
    ymd[indx,1] = ymd[indx,1]-1
    ymd[indx,2] = 12
  }
  end_last_month <- Dates_to_R(end_of_month(ymd))
  return(end_last_month)
}

#'End of Month Date
#'
#' Get date for the last day of month shift months ago (-) or ahead (+)
#'
#' @param last_date date in the month of interest
#' @param shift Months ahead (positive value) or behind (negative value)
#' @export
#' @useDynLib bdfm
End_of_Month <- function(last_date, shift = 0){
  r   <- length(last_date)
  ymd <- Dates_to_cpp(last_date) + matrix(1,r,1)%x%matrix(c(0,shift,0),1,3)
  if(any(ymd[,2] > 12)){
    indx   <- ymd[,2]>12
    exx    <- ymd[indx,2]
    yrs    <- floor(exx/12)
    mths   <- exx-12*yrs
    ymd[indx,1] = ymd[indx,1]+yrs
    ymd[indx,2] = mths
  }
  if(any(ymd[,2] < 1)){
    indx   <- ymd[,2]<1
    exx    <- -ymd[indx,2]
    yrs    <- floor(exx/12)
    mths   <- 12-(exx-12*yrs)
    ymd[indx,1] = ymd[indx,1]-yrs-1
    ymd[indx,2] = mths
  }

  end_next_month <- Dates_to_R(end_of_month(ymd))
  return(end_next_month)
}