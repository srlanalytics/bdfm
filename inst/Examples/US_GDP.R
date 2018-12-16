
library(bdfm) #estimation routines
library(tsbox) #used to format time series
library(jsonlite) #for calling the FRED API
library(data.table) #used to format data 

your_key = "9c39f8a4fbcc38c22475bc7f26988367"

# ------------------- Helper functions ------------------------------------------------------

# Unscale a set of data
unscale <- function(x, attrib){
  x*(rep(1,nrow(attrib))%x%t(attr(attrib, "scaled:scale"))) +
    rep(1,nrow(attrib))%x%t(attr(attrib, "scaled:center"))
}

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
est <- dfm(Y, factors = 2, lags = 3, nu_r = nu_r, reps = 1000, burn = 500, loud = TRUE)


ts.plot(Y)
max(Y, na.rm = T)

which(Y == max(Y))

#Make the data stationary by taking log differences. Note that we do not take logs of inventory:sales ratios (indexes 3 and 4), but we do difference them as they are not stationary.
data <- as.matrix(fred[,-1])
data[,c(1,2,5)] <- log(data[,c(1,2,5)])
data <- diff(data)
dates <- fred$DATE[-1] #keep corresponding dates (the first is dropped due to differencing)

#Scale and center data
data <- 100*scale(data)

#Estimate a Bayesian DFM
est <- dfm(data, factors = 1, lags = 3)

#Look at the estimated factor
ts.plot(cbind(data, est$factors), col = c(rep("steelblue", 5), "red"), lwd = c(rep(1,5),2))

#Get fitted values
predicted <- predict(est)

#Convert back to original (differenced) units
predicted <- unscale(predicted/100, data)

#How much does each observed series contribute to the factor?
diag(est$R)

#Use priors to update more aggressively on initial jobless claims
nu_r <- c(0,0,0,0,1)
est2  <- dfm(data, factors = 1, lags = 3, nu_r = nu_r)
diag(est2$R)

#Look at traces for a few of the parameters in our first estimation
par(mfrow=c(2,2))
ts.plot(est$Bstore[1,1,])
ts.plot(est$Hstore[1,1,])
ts.plot(est$Qstore[1,1,])
ts.plot(est$Rstore[1,])

#Forecast 3 periods ahead of last observation using xts time series format
data_xts <- xts(data, order.by = dates)
est_fct  <- dfm(data_xts,factors = 2, lags = 3, forecast = 3)
fct <- predict(est_fct)[,"INDPRO"]
print(tail(fct,4)) #note that the value for 2018-11-01 is a nowcast as IP is not observed but job claims are observed

#Store the posterior distribution of predicted industrial production values
est <- dfm(data, factors = 1, lags = 3, store_idx = 1)
#How do predicted values compare to actual values?
median_parameters <- mean((est$values[,1] - data[,1])^2, na.rm = T)
median_draw       <- mean((est$Ymedian[,1] - data[,1])^2, na.rm = T)

#What does the distribution of predicted values look like in period 321 (the second to last period)?
period = 321
hist(est$Ystore[period,], breaks = 30)

#How much did each series contribute to forecast updates in period 320?
Update <- c(est$Kstore[[320]][1,]*est$PEstore[[320]])
names(Update) <- colnames(data)
print(Update)

#Estimation and forecasting by Maximum Likelihood
est_ml <- dfm(data_xts, factors = 2, lags = 3, forecast = 3, method = 'ml')
fct <- predict(est_ml)[,"INDPRO"]
print(tail(fct,4))

#Two step estimation and forecast
est_pc <- dfm(data_xts, factors = 2, lags = 3, forecast = 3, method = 'pc')
fct <- predict(est_pc)[,"INDPRO"]
print(tail(fct,4))





