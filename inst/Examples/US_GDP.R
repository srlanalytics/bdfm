
library(bdfm) #estimation routines
library(tsbox) #used to format time series
library(jsonlite) #for calling the FRED API
library(data.table) #used to format data 

your_key = "9c39f8a4fbcc38c22475bc7f26988367"

# ------------------- Helper functions ------------------------------------------------------

# A handy little utility to call data from the FRED API. 
get_fred_lite <- function(series_id, series_name = NULL, observation_start = "1776-07-04", observation_end = "9999-12-31",frequency = NULL){
  fred_call <- paste0("https://api.stlouisfed.org/fred/series/observations?series_id=",
                      series_id,
                      "&observation_start=",
                      observation_start,
                      "&observation_end=",
                      observation_end,
                      "&output_type=1",
                      "&api_key=",
                      your_key,
                      "&file_type=json")
  if(!is.null(frequency)){
    fred_call <- paste0(fred_call,"&frequency=",frequency)
  }
  X <- fromJSON(fred_call)$observations[,c("date", "value")]
  X[X$value == ".","value"] <- NA 
  X <- data.frame(date = as.Date(X$date), value = as.numeric(X$value))
  if(is.null(series_name)){
    series_name = series_id
  }
  X <- cbind.data.frame(X$date, rep(series_name,nrow(X)), X$value)
  colnames(X) <- c("date", "series_name", "value")
  return(X)
}

# Run the actual API call and 
Get_Data <- function(series_name){
  if(series_name%in%c('T10Y3M', 'ICSA', 'TWEXB')){
    frq  <- 'm' #Aggregate these up to monthly
    Data <- get_fred_lite(series_id = series_name, observation_start = "1980-01-01", frequency = frq)
  }else{
    Data <- get_fred_lite(series_id = series_name, observation_start = "1980-01-01")
  }
  return(Data)
} 

# Go from long to mixed frequency wide format. We can also use as_of for backtesting. 
# This function requires data.table
call_data <- function(..., dt = Data, as_of = NULL, date = "date", series_id = "series_id", value = "value", pub_date = NULL){
  
  dt <- data.table(dt)

  colnames(dt)[colnames(dt)==date]        <- "date"
  colnames(dt)[colnames(dt)==series_id]   <- "series_id"
  colnames(dt)[colnames(dt)==value]       <- "value"
  if(is.null(pub_date)){
    dt[,pub_date := as.Date(date)]
  }else{
    colnames(dt)[colnames(dt)==pub_date]    <- "pub_date"
  }
  
  setkey(dt, series_id, date)
  
  # for convenience, so we can write call_data("dfsd", "dsfsdf")
  l <- list(...)
  is.char <- (sapply(l, class) == "character")
  if (any(!is.char)){
    stop("some elements are not character vectors")
  }
  oi <- do.call("c", l) #of interest
  
  unique_names <- unique(dt[,series_id])
  missing_series <- oi[!oi%in%unique_names]
  
  if(length(missing_series)>0){
    warning(paste("The following series are missing:", paste(missing_series, collapse = ", ")))
  }
  
  if(is.null(as_of)) as_of <- Sys.time() #Model from current time (default) or set as_of for backtesting
  
  Out <- dcast(dt[series_id%in%oi & (pub_date <= as.Date(as_of) | (is.na(pub_date) & as.Date(date) <= as.Date(as_of)) ) ], date ~ series_id,    value.var = "value")
  
  return(Out)
}

#Aggregate monthly data allowing us to correctly specify our factor model with bdfm
#Just use a simple moving average for monthly data
moving_ave <- function(Y,n){filter(Y, filter = rep(1/n,n), method = "convolution", sides = 1)}

diffMQ <- function(dt){
  
  #get frequency of each series. In this example 3 is quarterly (3 high frequency observations
  #per low frequency observation) and 1 is monthly. 
  freq <- rep(0,ncol(dt)-1)
  for(j in 2:ncol(dt)){
    freq[j-1] <- median(diff(which(is.finite(unlist(dt[,j,with = F]))))) 
  }
  
  
  
}

# -----------------------------------------------------------------------------------

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

# -------- Aggregate monthly data for a correctly specified mixed frequency model ----------------

# The first two series are quarterly, everything else is monthly, so omit them in our aggregations
Y[, 3:24] <- moving_ave(Y[, 3:24], 3)
k <- ncol(Y) #number of series 

# --------- Take logs where necessary ------------------------------------------------------------
#Specify which variables we should not take logs of
ind_log    <- 1:k
no_log     <- c("A191RL1Q225SBEA","USSLIND", "UMCSENT", "UNRATE", "FRGSHPUSM649NCIS", "CAPUTLG2211S", "MNFCTRIRSA", "RETAILIRSA", "WHLSLRIRSA", "TWEXB", "T10Y3M")
ind_log    <- ind_log[!ind_log%in%unlist(sapply(no_log, FUN = grep, colnames(Y)))]
Y[,ind_log] <- log(Y[,ind_log])

# ----Take differences where necessary, again ensuring correct mixed frequency specification---------
ind_diff   <- 1:k
no_diff    <- c("A191RL1Q225SBEA", "USSLIND")
ind_diff   <- ind_diff[!ind_diff%in%unlist(sapply(no_diff, FUN = grep, colnames(Y)))]
Y[-(1:3), ind_diff] <- diff(Y[,ind_diff], lag = 3)

#Drop first two quarters due to aggregation and differencing 

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
Y <- 100*scale(Y)

# ------- Estimate a DFM and target GDP ----------------------------------------------------
nu_r = rep(0,k)
nu_r[1] <- 1 #target GDP using our prior deg. of freedom on shocks to obs. of GDP
#(GDP is the first series, hence the index 1)
est <- dfm(Y, factors = 3, lags = 3, nu_r = nu_r, store_idx = 1, reps = 1000, burn = 500, loud = TRUE)

# #mixed frequency call
# frequency_mix <- rep(1,ncol(Y)) #1 indicates high frequency data --- 1 period per observation
# frequency_mix[1:2] <- 3 #4 for quarterly data in this case --- 3 months in a quarter
# 
# differences <- rep(0,ncol(Y))
# differences[1:2] <- 1 # the first two low frequency series are (log) differenced
# 
# est <- dfm(Y, factors = 3, lags = 3, frequency_mix = frequency_mix, differences = differences,
#            nu_r = nu_r, store_idx = 1, reps = 1000, burn = 500, loud = TRUE)


#Are we drawing from a stationary distribution?
ts.plot(est$Qstore[1,1,]) 
ts.plot(est$Hstore[1,1,]) 

# -------- Convert results back to original units -------------------------------------

Ymedian <- attr(Y, "scaled:scale")[1]*est$Ymedian/100 + attr(Y, "scaled:center")[1]

# Looking at results at high frequency (monthly)

plot(dt$date, Ymedian, type = 'l', col = 'red', xlab = "year", ylab = "GDP")
points(dt$date, dt$A191RL1Q225SBEA, col = 'steelblue')

# Looking at results at low frequency (quarterly)
Qind <- format(dt$date, "%m")%in%c('03','06','09','12')
plot(dt$date[Qind], dt$A191RL1Q225SBEA[Qind], type = 'l', col = 'steelblue', xlab = 'year', ylab = 'GDP')
lines(dt$date[Qind], Ymedian[Qind], col = 'red', lty = 2)
# Our prediction for the next GDP release is three periods ahead of the last release

prediction <- Ymedian[max(which(is.finite(Y[,1])))+3]





