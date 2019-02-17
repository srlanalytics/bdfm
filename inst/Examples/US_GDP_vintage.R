#Third approach to mixed frequncy data:

# This script follows the "standard" approach to mixed frequency data
# of Mariano and Murasawa (2003)

library(bdfm) #estimation routines
library(tsbox) #used to format time series
library(jsonlite) #for calling the FRED API
library(data.table) #used to format data

# Use Sys.setenv(FRED_API_KEY = 'XXXX') to set FRED API Key.
# Optain here:https://research.stlouisfed.org/docs/api/api_key.html

your_key = Sys.getenv("FRED_API_KEY")

# ------------------- Helper functions ------------------------------------------------------

# A handy little utility to call data from the FRED API.
get_fred_lite <- function(series_id, series_name = NULL, observation_start = "1776-07-04", observation_end = "9999-12-31",frequency = NULL, vintage_dates = NULL){
  if(!is.null(vintage_dates)){
    observation_end = max(as.Date(vintage_dates))
    output_type = 2
  }else{
    output_type = 1
  }
  fred_call <- paste0("https://api.stlouisfed.org/fred/series/observations?series_id=",
                      series_id,
                      "&observation_start=",
                      observation_start,
                      "&observation_end=",
                      observation_end,
                      "&output_type=",
                      output_type,
                      "&api_key=",
                      your_key,
                      "&file_type=json")

  if(!is.null(frequency)){
    fred_call <- paste0(fred_call,"&frequency=",frequency)
  }
  if(!is.null(vintage_dates)){
    fred_call <- paste0(fred_call, "&vintage_dates=", vintage_dates)
  }
  X <- fromJSON(fred_call)
  if("value"%in%colnames(X$observations)){
    X$observations <- X$observations[, c("date", "value")]
  }else if(length(grep(series_id, colnames(X$observations)))!=0 ){
    X$observations = data.frame(date = X$observations[,"date"],
                                value = X$observations[,grep(series_id, colnames(X$observations))])
  }
  X$observations[X$observations[,"value"] == ".","value"] <- NA
  if(is.null(series_name)){
    series_name = series_id
  }
  X$observations <- data.frame(date = as.Date(X$observations[,"date"]),
                               series_name = rep(series_name,nrow(X$observations)),
                               value = as.numeric(X$observations[,"value"]) )
  return(X)
}

# Run the actual API call and
Get_Data <- function(series_name, vintage_dates = NULL){
  if(is.null(vintage_dates)){
    if(series_name%in%c('T10Y3M', 'ICSA', 'TWEXB')){
      frq  <- 'm' #Aggregate these up to monthly
      Data <- get_fred_lite(series_id = series_name, observation_start = "1980-01-01", frequency = frq)
    }else{
      Data <- get_fred_lite(series_id = series_name, observation_start = "1980-01-01")
    }
  }else{
    if(series_name%in%c('T10Y3M', 'ICSA', 'TWEXB')){
      frq  <- 'm' #Aggregate these up to monthly
      Data <- get_fred_lite(series_id = series_name, observation_start = "1980-01-01", frequency = frq, vintage_dates = vintage_dates)
    }else{
      Data <- get_fred_lite(series_id = series_name, observation_start = "1980-01-01", vintage_dates = vintage_dates)
    }
  }
  return(Data$observations)
}

# Go from long to mixed frequency wide format. We can also use as_of for backtesting.
# This function requires data.table
call_data <- function(series_names, dt = Data, as_of = NULL, date = "date", series_id = "series_id", value = "value", pub_date = NULL){

  dt <- data.table(dt)

  colnames(dt)[colnames(dt)==date]        <- "date"
  colnames(dt)[colnames(dt)==series_id]   <- "series_id"
  colnames(dt)[colnames(dt)==value]       <- "value"
  if(is.null(pub_date)){
    dt[,pub_date := as.Date(NA)]
  }else{
    colnames(dt)[colnames(dt)==pub_date]    <- "pub_date"
  }

  unique_names <- unique(dt$series_id)
  missing_series <- series_names[!series_names%in%unique_names]

  if(length(missing_series)>0){
    warning(paste("The following series are missing:", paste(missing_series, collapse = ", ")))
  }

  if(is.null(as_of)) as_of <- Sys.time() #Model from current time (default) or set as_of for backtesting

  Out <- dcast(dt[series_id%in%series_names & (pub_date <= as.Date(as_of) | (is.na(pub_date) & as.Date(date) <= as.Date(as_of)) ) ], date ~ series_id,    value.var = "value")

  return(Out)
}


# ------------------------ Step 1: Import Data from FRED API ------------------------

unique_names <- c('A191RL1Q225SBEA', #1 Real GDP, seasonally adjusted, quarterly, annualized % change
                  'W068RCQ027SBEA', #2 Governemnt expenditures
                  'USSLIND', #3 Federal Reserve leading index, monthly, percent
                  'PCEDG', #4 persional consumption: durable goods, monthly, level
                  'PCEND', #5 persional consumption: non-durable goods, monthly, level
                  'UMCSENT', #6 Consumer Sentiment, monthly, delayed 1 month for free data
                  'UNRATE', #7 Unemployment, monthly
                  'JTSJOL', #8 Job openenings, total non-farm
                  'INDPRO', #9 Industrial Production Index, monthly, level
                  'CSUSHPINSA', #10 Case-Shiller home price index, monthly, two month lag, level
                  'HSN1F',     #11 New 1 family houses sold, level
                  'TSIFRGHT', #12 Freight transportation index, monthly, 2-3 month lag, level
                  'FRGSHPUSM649NCIS', #13 CASS freight index, level, not SA
                  'CAPUTLG2211S', #14 Electricity usage, % capacity, monthly
                  'IPG2211S', #15 Electricity, industrial production index, monthly, level
                  'DGORDER', #16 New Orders, durable manufacturing goods, monthly, level
                  'AMTMNO', #17 New Orderes, all manufacuring industries, level
                  'MNFCTRIRSA', #18 Manufacturers inventories:sales ratio
                  'RETAILIRSA', #19 Retail inventories:sales ratio
                  'WHLSLRIRSA', #20 Wholesalers, inventories:sales ratio
                  'CPILFESL', #21 CPI
                  'ICSA', #22 Initial claims, SA, weekly
                  'TWEXB', #23 exchange rate index, weekly
                  'T10Y3M') #24 10Y to 3M treasury spread, daily

RawData <- lapply(unique_names, FUN = Get_Data)
cat("Retrieved Fred Data")
Data    <- rbindlist(RawData)
dt      <- call_data(unique_names, dt = Data, series_id = "series_name") #mixed frequency format
dt      <- setcolorder(dt, c("date", unique_names))


#Quarterly data MUST be indexed at the end of the quarter, not the first month of the quarter!
dt[, A191RL1Q225SBEA := shift(A191RL1Q225SBEA, n = 2, type = "lag")]
dt[, W068RCQ027SBEA := shift(W068RCQ027SBEA, n = 2, type = "lag")]

#It's worth looking at dt. This panel is far from square and includes quarterly and monthly data.
#This is, there are loads of missing observations. For our DFM, that's no problem!

# ------------------------------- Seasonal Adjustment ---------------------------------------------

# A lot of our data is already seasonally adjusted. We actually only need to SA three series
# (though whether it's necessary to adjust consumer sentiment is debatable).

do_SA   <- c('UMCSENT', 'CSUSHPINSA', 'FRGSHPUSM649NCIS')
ind_SA  <- unlist(sapply(do_SA, FUN = grep, unique_names))
for(j in ind_SA){
  tmp_ts  <- ts_ts(dt[,c(1, j+1), with=F]) #convert to ts for sesaonal adjustment
  sa      <- seas_we(tmp_ts, lags = 3, verbose = T)
  tmp_dt  <- ts_data.table(sa$values)
  dt[dt$date%in%tmp_dt$time, j+1] <- tmp_dt$value
}

# Work with matrices from now on to make life easier
Y <- as.matrix(dt[,-1])

# Y.pkg <- `attr<-`(unclass(econ_us), "tsp", NULL)
# all.equal(Y.pkg, Y[-470,])  # TRUE (2019-02-17)

k <- ncol(Y) #number of series

# --------- Take logs where necessary ------------------------------------------------------------
#Specify which variables we should not take logs of
ind_log    <- 1:k
no_log     <- c("A191RL1Q225SBEA","USSLIND", "UMCSENT", "UNRATE", "FRGSHPUSM649NCIS", "CAPUTLG2211S", "MNFCTRIRSA", "RETAILIRSA", "WHLSLRIRSA", "TWEXB", "T10Y3M")
ind_log    <- ind_log[!ind_log%in%unlist(sapply(no_log, FUN = grep, colnames(Y)))]
Y[,ind_log] <- log(Y[,ind_log])

# ----Take differences where necessary, again ensuring correct mixed frequency specification---------
#Difference monthly data
ind_diffM   <- 1:k
no_diff     <- c("A191RL1Q225SBEA", 'W068RCQ027SBEA', "USSLIND")
ind_diffM   <- ind_diffM[!ind_diffM%in%unlist(sapply(no_diff, FUN = grep, colnames(Y)))]
Y[-1, ind_diffM] <- diff(Y[,ind_diffM], lag = 1) #take differences at 1 lag
#Difference Quarterly Data
Y[-(1:3),'W068RCQ027SBEA' ] <- diff(Y[,'W068RCQ027SBEA'], lag = 3)

#Drop first two quarters for consistency
Y  <- Y[-(1:6), ]
dt <- dt[-(1:6), ] #keep dimensions the same everywhere

# ----- Drop outliers (optional but gets rid of some wierd stuff) ------------------------
for(j in 1:k){
  itc <- mean(Y[,j],na.rm = TRUE)
  SD  <- sqrt(mean((Y[,j]-itc)^2,na.rm = TRUE))
  Err <- abs(Y[,j]- itc)/SD #eliminate observations for which shocks are greater than 5 s.d
  Y[Err>4,j] <- NA #drop obs more than 4 s.d. from mean
}

# ------- Scale the Data (optional but usually a good idea) --------------------------------
Y <- scale(Y) #note that this is bdfm::scale() so values are stored for unscaling

# ------- Estimate a DFM and target GDP ----------------------------------------------------
#Specify frequencies
freq <- rep(1,k) #All but the first 2 series are monthly --- 3 months in a quarter
freq[1:2] <- 3 #Specifying which are quarterly (the first two)
#Specify differences vs. levels
differences <- rep(1,k) #every quarterly series is differenced (GDP is already in
# percent change, i.e. log differences); this only matters for quarterly data.


# Target GDP using our prior deg. of freedom on shocks to obs. of GDP
# (GDP is the first series, hence the index 1)
nu_r = rep(0,k)
nu_r[1] <- 1
#Estimate the mixed frequency dfm
est <- dfm(Y, factors = 3, lags = 3, frequency_mix = freq, differences = differences, nu_r = nu_r,
           identification = "PC_full", store_idx = 1, reps = 1000, burn = 500, loud = TRUE)

#Are we drawing from a stationary distribution?
ts.plot(est$Qstore[1,1,])
ts.plot(est$Hstore[1,1,])

summary(est)

est <- dfm(Y, factors = 3, lags = 3, method = "pc", frequency_mix = freq, differences = differences, nu_r = nu_r,
           identification = "PC_sub", store_idx = 1, reps = 1000, burn = 500, loud = TRUE)


yy <- Y[,c(3,4,5,6,7,9)]

est2 <- dfm(yy)

print(est2)


# -------- Convert results back to original units -------------------------------------

Ymedian  <- unscale(est$Ymedian, idx = "A191RL1Q225SBEA")

# Looking at results at high frequency (monthly)

plot(dt$date, Ymedian, type = 'l', col = 'red', xlab = "year", ylab = "GDP")
points(dt$date, dt$A191RL1Q225SBEA, col = 'steelblue')

# Looking at results at low frequency (quarterly)
Qind <- format(dt$date, "%m")%in%c('03','06','09','12')
plot(dt$date[Qind], dt$A191RL1Q225SBEA[Qind], type = 'l', col = 'steelblue', xlab = 'year', ylab = 'GDP')
lines(dt$date[Qind], Ymedian[Qind], col = 'red', lty = 2)
# Our prediction for the next GDP release is three periods ahead of the last release
prediction <- Ymedian[max(which(is.finite(Y[,1])))+3]



