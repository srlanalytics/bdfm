

#--- creating matrix of predetermined variables ---
# This example includes a few possiblities for predetermined variables.
# It should be simple to include others following the same format.

DayOfWeek <- function(Days, dates) {
  # day_of_week <- seq.Date(from = as.Date(Begin), to = as.Date(End), by = "day")
  day_of_week <- weekdays(dates)
  day_dummy <- matrix(0, length(day_of_week), length(Days))
  for (j in 1:length(Days)) {
    day_dummy[day_of_week %in% Days[j], j] <- 1
  }
  return(day_dummy)
}
sum_yearmonth <- function(yearmonth, daily_seq, x){
  ind <- which(format.Date(daily_seq, "%Y-%m")%in%yearmonth)
  out <- sum(x[ind])
  return(out)
}

#' Seasonal Factors at Daily Frequency
#'
#' Generate predetermined seasonal factors at daily frequency
#'
#' @param dates dates that factors will be generated for (date vector)
#' @param predetermined seasonal adjustment factors to be generated, for example c('June', 'July')
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib BDFM
Predetermined.d <- function(dates, predetermined) {
  mn <- length(predetermined)
  nn <- matrix(0, length(dates), mn)
  indx <- 1

  # Constant intercept term
  if ("int" %in% predetermined) {
    nn[, indx] <- 1
    indx <- indx + 1
  }

  # Day of the week effect for daily data
  Days <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
  if (any(Days %in% predetermined)) {
    Days <- Days[Days %in% predetermined]
    m_tmp <- length(Days)
    nn[, indx:(indx + m_tmp - 1)] <- DayOfWeek(Days, dates)
    indx <- indx + m_tmp
  }


  #--------- Seasons -----------------------

  # Need to re-do this as it doesn't work for beginning and end of sample!!

  # Solution: add 60 days to either end of the date vector to make sure seasonal effects for beginning and end of the sample are included.

  dates_ex_long <- seq.Date(from = as.Date(head(dates, 1) - 60), to = as.Date(tail(dates, 1) + 60), by = "d")

  if ("spring" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "04-16") # spring dates
    if (!length(H_dates) == 0) {
      spn <- -46:45
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  if ("summer" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "07-16") # summer dates
    if (!length(H_dates) == 0) {
      spn <- -45:46
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  if ("fall" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "10-16") # fall dates
    if (!length(H_dates) == 0) {
      spn <- -45:45
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  if ("winter" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "01-15") # winter dates
    if (!length(H_dates) == 0) {
      spn <- -45:44
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  # ------- Months -------------------------------------

  if ("January" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "01-15") # January dates
    if (!length(H_dates) == 0) {
      spn <- -20:20
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  # jim <- cbind.data.frame(dates,nn)

  if ("February" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "02-14") # February dates
    if (!length(H_dates) == 0) {
      spn <- -19:19
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  if ("March" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "03-15") # March dates
    if (!length(H_dates) == 0) {
      spn <- -20:20
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  if ("April" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "04-15") # April dates
    if (!length(H_dates) == 0) {
      spn <- -20:20
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  if ("May" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "05-15") # May dates
    if (!length(H_dates) == 0) {
      spn <- -20:20
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  if ("June" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "06-15") # June dates
    if (!length(H_dates) == 0) {
      spn <- -20:20
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  if ("July" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "07-15") # July dates
    if (!length(H_dates) == 0) {
      spn <- -20:20
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  if ("August" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "08-15") # August dates
    if (!length(H_dates) == 0) {
      spn <- -20:20
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  if ("September" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "09-15") # September dates
    if (!length(H_dates) == 0) {
      spn <- -20:20
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  if ("October" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "10-15") # October dates
    if (!length(H_dates) == 0) {
      spn <- -20:20
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  if ("November" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "11-15") # November dates
    if (!length(H_dates) == 0) {
      spn <- -20:20
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  if ("December" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "12-15") # February dates
    if (!length(H_dates) == 0) {
      spn <- -20:20
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .4) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }

  # Holiday effects


  # Christmas is the only holiday I've put in so far but any other holiday can use the same format.
  if ("Christmas" %in% predetermined) {
    month_day <- format.Date(dates_ex_long, "%m-%d")
    H_dates <- which(month_day == "12-25") # Holiday dates
    if (!length(H_dates) == 0) {
      spn <- -54:0 # for Christmas go back to November 1st
      nm <- dnorm(spn / max(abs(spn)), mean = 0, sd = .25) # effect uses the shape of a normal pdf
      nm <- 10 * nm / sum(nm)
      H_span <- sort(unlist(lapply(H_dates, function(bob) {
        bob + spn
      }))) # indexes to span (Holiday span)
      H_span <- H_span[H_span > 0 & H_span <= length(dates_ex_long)]
      efx <- rep(nm, length(H_dates))
      H_dates <- dates_ex_long[H_span]
      ind_efx <- which(H_dates %in% dates)
      ind_date <- which(dates %in% H_dates)
      nn[ind_date, indx] <- efx[ind_efx]
    }
    indx <- indx + 1
  }
  if("Easter"%in%predetermined || "CNY"%in%predetermined || "Diwali"%in%predetermined){
    load(system.file("data/Holiday.RData", package = "BDFM"))
  }

  if("CNY"%in%predetermined){ #Chinese New Year
    #month_day   <- format.Date(dates_ex_long, "%m-%d")
    H_dates     <- which(dates_ex_long%in%(holiday$cny+7))  #Mid of CNY
    if(!length(H_dates)==0){
      spn         <- -7:7
      nm          <- dnorm(spn/max(abs(spn)), mean = 0, sd = .5) #effect uses the shape of a normal pdf
      nm          <- -10*nm/sum(nm)
      H_span      <- sort(unlist(lapply(H_dates,function(bob){bob+spn}))) #indexes to span (Holiday span)
      H_span      <- H_span[H_span>0 & H_span<=length(dates_ex_long)]
      efx         <- rep(nm,length(H_dates))
      H_dates     <- dates_ex_long[H_span]
      ind_efx     <- which(H_dates%in%dates)
      ind_date    <- which(dates%in%H_dates)
      nn[ind_date, indx]  <- efx[ind_efx]
    }
    R_dates     <- which(dates_ex_long%in%(holiday$cny+40))  #Mid of CNY
    if(!length(H_dates)==0){
      spn         <- -25:25
      nm          <- dnorm(spn/max(abs(spn)), mean = 0, sd = .5) #effect uses the shape of a normal pdf
      nm          <- 10*nm/sum(nm)
      R_span      <- sort(unlist(lapply(R_dates,function(bob){bob+spn}))) #indexes to span (Recovery span)
      R_span      <- R_span[R_span>0 & R_span<=length(dates_ex_long)]
      efx         <- rep(nm,length(R_dates))
      R_dates     <- dates_ex_long[R_span]
      ind_efx     <- which(R_dates%in%dates)
      ind_date    <- which(dates%in%R_dates)
      nn[ind_date, indx]  <- efx[ind_efx]
    }
    indx <- indx+1
  }

  if("Easter"%in%predetermined){
    #month_day   <- format.Date(dates_ex_long, "%m-%d")
    H_dates     <- which(dates_ex_long%in%holiday$easter)  #Index of dates
    if(!length(H_dates)==0){
      spn         <- -5:10
      nm          <- dnorm(spn/max(abs(spn)), mean = 0, sd = .4) #effect uses the shape of a normal pdf
      nm          <- 10*nm/sum(nm)
      H_span      <- sort(unlist(lapply(H_dates,function(bob){bob+spn}))) #indexes to span (Holiday span)
      H_span      <- H_span[H_span>0 & H_span<=length(dates_ex_long)]
      efx         <- rep(nm,length(H_dates))
      H_dates     <- dates_ex_long[H_span]
      ind_efx     <- which(H_dates%in%dates)
      ind_date    <- which(dates%in%H_dates)
      nn[ind_date, indx]  <- efx[ind_efx]
    }
    indx <- indx+1
  }

  if("Diwali"%in%predetermined){
    #month_day   <- format.Date(dates_ex_long, "%m-%d")
    H_dates     <- which(dates_ex_long%in%holiday$diwali)  #Index of dates
    if(!length(H_dates)==0){
      spn         <- -10:10
      nm          <- dnorm(spn/max(abs(spn)), mean = 0, sd = .5) #effect uses the shape of a normal pdf
      nm          <- 10*nm/sum(nm)
      H_span      <- sort(unlist(lapply(H_dates,function(bob){bob+spn}))) #indexes to span (Holiday span)
      H_span      <- H_span[H_span>0 & H_span<=length(dates_ex_long)]
      efx         <- rep(nm,length(H_dates))
      H_dates     <- dates_ex_long[H_span]
      ind_efx     <- which(H_dates%in%dates)
      ind_date    <- which(dates%in%H_dates)
      nn[ind_date, indx]  <- efx[ind_efx]
    }
    indx <- indx+1
  }
  rownames(nn) <- as.character(dates)
  return(nn)
}

#' Seasonal Factors at Monthly Frequency
#'
#' Generate predetermined seasonal factors at monthly frequency
#'
#' @param dates dates that factors will be generated for (date vector)
#' @param predetermined seasonal adjustment factors to be generated, for example c('June', 'July')
#' @export
#' @import data.table
#' @useDynLib BDFM
#' @examples
#' library(BDFM)
#' Spain <- read.csv(
#'   system.file("Examples/Spain_IP.csv", package = "BDFM"),
#'   stringsAsFactors = FALSE
#' )
#' y <- Spain$IP
#' dates <- as.Date(Spain$date)
#' predetermined = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December", "Easter", "trading_days")
#' N <- Predetermined.m(dates, predetermined)
Predetermined.m <- function(dates, predetermined) {
  mn <- length(predetermined)
  N <- matrix(0, length(dates), mn)
  indx <- 1

  # Constant intercept term
  if ("int" %in% predetermined) {
    N[, indx] <- 1
    indx <- indx + 1
  }
  if ("spring" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "03", indx] <- .5
    N[Month == "04", indx] <- 1
    N[Month == "05", indx] <- .5
    indx <- indx + 1
  }
  if ("summer" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "06", indx] <- .5
    N[Month == "07", indx] <- 1
    N[Month == "08", indx] <- .5
    indx <- indx + 1
  }
  if ("fall" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "09", indx] <- .5
    N[Month == "10", indx] <- 1
    N[Month == "11", indx] <- .5
    indx <- indx + 1
  }
  if ("winter" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "12", indx] <- .5
    N[Month == "01", indx] <- 1
    N[Month == "02", indx] <- .5
    indx <- indx + 1
  }
  # if("trading_days"%in%predetermined){
  #   daily_seq <- seq.Date(from = as.Date(paste(format.Date(dates[1], "%Y-%m"), "01", sep = "-")), to = End_This_Month(tail(dates,1)), by = "day")
  #   weekday_count <- rep(1, length(daily_seq))
  #   weekday_count[weekdays(daily_seq)%in%c("Saturday", "Sunday")] <- 0
  #   yearmonth <- unique(format.Date(daily_seq, "%Y-%m"))
  #   N[,indx] <- sapply(X = yearmonth, FUN = sum_yearmonth, daily_seq = daily_seq, x = weekday_count)
  #   indx <- indx+1
  # }

  # browser()
  #Fast data.table code, but doesn't seem to work when building package:

  if("trading_days"%in%predetermined){
    daily_sequence <- seq.Date(from = as.Date(paste(year(dates[1]), month(dates[1]), "01", sep = "-")), to = End_This_Month(tail(dates,1)), by = "day")
    weekday_count <- rep(1, length(daily_sequence))
    weekday_count[weekdays(daily_sequence)%in%c("Saturday", "Sunday")] <- 0
    tmp_dt <- data.table(dates = daily_sequence, count = weekday_count)
    out <- tmp_dt[, lapply(.SD,sum), by = .(year(dates), month(dates))]
    N[,indx] <- out$count
    indx <- indx+1
  }
  if ("January" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "01", indx] <- 1
    indx <- indx + 1
  }
  if ("February" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "02", indx] <- 1
    indx <- indx + 1
  }
  if ("March" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "03", indx] <- 1
    indx <- indx + 1
  }
  if ("April" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "04", indx] <- 1
    indx <- indx + 1
  }
  if ("May" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "05", indx] <- 1
    indx <- indx + 1
  }
  if ("June" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "06", indx] <- 1
    indx <- indx + 1
  }
  if ("July" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "07", indx] <- 1
    indx <- indx + 1
  }
  if ("August" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "08", indx] <- 1
    indx <- indx + 1
  }
  if ("September" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "09", indx] <- 1
    indx <- indx + 1
  }
  if ("October" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "10", indx] <- 1
    indx <- indx + 1
  }
  if ("November" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "11", indx] <- 1
    indx <- indx + 1
  }
  if ("December" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "12", indx] <- 1
    indx <- indx + 1
  }
  if ("Christmas" %in% predetermined) {
    Month <- format.Date(dates, "%m")
    N[Month == "12", indx] <- 1
    indx <- indx + 1
  }
  if("Easter"%in%predetermined || "CNY"%in%predetermined || "Diwali"%in%predetermined){
    load(system.file("data/holiday.RData", package = "BDFM"))
  }
  if("Easter"%in%predetermined){
    Month     <- format.Date(dates, "%Y-%m")
    HMonth    <- format.Date(holiday$easter, "%Y-%m")
    H_dates   <- which(Month%in%HMonth)
    N[H_dates, indx] <- 1
    indx <- indx+1
  }
  if("CNY"%in%predetermined){
    Month     <- format.Date(dates, "%Y-%m")
    HMonth    <- format.Date(holiday$cny+7, "%Y-%m") #end of the holiday
    RMonth    <- format.Date(End_Next_Month(holiday$cny+7), "%Y-%m") #recovery period
    H_dates   <- which(Month%in%HMonth)
    R_dates   <- which(Month%in%RMonth)
    N[H_dates, indx] <- -1
    N[R_dates, indx] <- 1
    indx <- indx+1
  }
  if("Diwali"%in%predetermined){
    Month     <- format.Date(dates, "%Y-%m")
    HMonth    <- format.Date(holiday$diwali, "%Y-%m")
    H_dates   <- which(Month%in%HMonth)
    N[H_dates, indx] <- 1
    indx <- indx+1
  }
  return(N)
}

#' Watson and Engle (1983) Algorithm
#'
#' to produce maximum likelihood estimates of model parameters
#'
#' @param y    data
#' @param N    matrix of predetermined factors
#' @param lags    lags in transition equation
#' @param tol  tolerance for convergence of likelihood function (*100)
#' @param Loud T/F whether to output convergence of iterations
#' @export
#' @useDynLib BDFM
SeasAdj_WE <- function(y, N, lags = 1, tol = 0.01, Loud = FALSE) {
  p <- lags #shorthand notation
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

  Lik0 <- -(1e10)
  count <- 0
  Conv <- 100
  while (Conv > tol | count < 6) { # loop until the likelihood function converges
    Est <- KSeas(B, q, M, r, y, N) # calculate estimates in C++
    B <- Est$B
    q <- Est$q
    r <- Est$r
    M <- Est$M
    Z0 <- Est$Z0
    Lik1 <- Est$Lik
    Conv <- 100 * (Lik1 - Lik0) / abs(Lik1 + Lik0)
    Lik0 <- Lik1
    if (Loud) {
      print(Conv)
    }

    count <- count + 1
  }

  y_SA <- y - N %*% t(M)
  sa <- N %*% t(M)

  return(list(
    y_SA = y_SA, # seasonally adjusted Y
    sa = sa, # seasonal adjustments
    M = M
  )) # adjustment parameter
}
