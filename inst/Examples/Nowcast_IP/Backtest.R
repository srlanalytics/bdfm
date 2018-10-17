library(jsonlite)
library(data.table)
library(BDFM)

setwd("...")
source("./Nowcast_IP/helper_functions.R")

unique_names <- c('INDPRO', #1 Industrial Production Index, monthly, level
                  'USSLIND', #2 US Leading Index --- this will be adjusted forward one period
                  'UNRATE', #3 Unemployment, monthly
                  'FRGSHPUSM649NCIS', #4 CASS freight index, level, not SA --- comes out one day ahead
                  'CPILFESL', #5 CPI
                  'ICSA', #6 Initial claims, SA, weekly
                  'WILL5000IND', #7 Willshire 5000, daily, level
                  'TWEXB', #8 exchange rate index, weekly  
                  'T10Y3M') #9 10Y to 3M treasury spread, daily


Type <- c('Production', #1 Industrial Production
          'Survey', #2 US Leading Index
          'Employment', #3 Employment
          'Trade', #4 CASS freight index
          'Prices', #5 CPI
          'Employment', #6 Initial Claims
          'Prices', #7 Willshire 5000
          'Prices', #8 Exchange Rate
          'Prices') #9 treasury spread

#Taking logs and differences is a matter of judgement --- feel free to mess around

take_log <- c(T, #1 Industrial Production
              F, #2 US Leading Index
              F, #3 unemployment
              F, #4 CASS freight index
              T, #5 CPI
              T, #6 Initial Claims
              T, #7 Willshire 5000
              T, #8 Exchange Rate
              F) #9 treasury spread

take_diff <- c(T, #1 Industrial Production
              F, #2 US Leading Index
              T, #3 unemployment
              T, #4 CASS freight index
              T, #5 CPI
              T, #6 Initial Claims
              T, #7 Willshire 5000
              T, #8 Exchange Rate
              F) #9 treasury spread

late_obs <-  c(T, #1 Industrial Production
               F, #2 US Leading Index
               F, #3 unemployment
               T, #4 CASS freight index
               F, #5 CPI
               F, #6 Initial Claims
               F, #7 Willshire 5000
               F, #8 Exchange Rate
               F) #9 treasury spread

#get a FRED api key at https://research.stlouisfed.org/docs/api/api_key.html
api_key = "XXX"

Get_Data <- function(series_name){
    Data <- get_fred_data(series_id = series_name, api_key = api_key, observation_start = "1980-01-15", observation_end = "2018-10-15", realtime_start = "2018-10-15", realtime_end = "2018-10-15", frequency = "m")
  return(Data)
}   

RawData <- lapply(unique_names, FUN = Get_Data) #get data from the Fred API
Data    <- rbindlist(RawData) #stack observations into a data.table
dt      <- call_data(unique_names, dt = Data, series_id = "series_name") #put data into matrix format
dt      <- setcolorder(dt, c("date", unique_names)) #make sure columns are in the right order
dates   <- dt$date #put dates into a date vector
Y       <- as.matrix(dt[,-1, with = F]) #data matrix that will enter BDFM()
Y[,take_log] <- log(Y[,take_log]) #take logs
dY           <- Y[-1,] #dY will be filled in with differenced values
dates        <- dates[-1]
dY[,'USSLIND']  <- Y[-nrow(Y),'USSLIND'] #push the leading index forward one period
dY[,take_diff]  <- diff(Y[,take_diff]) #take differences

backtest_dates <- seq.Date(from = as.Date("2000-10-15"), to = as.Date("2018-10-15"), by = "m")

# ---- function to run backtests given the testing date ------------------

run_nowcast <- function(test_date){

#Drop data that would not have been observed
tDates <- dates[dates<=as.Date(test_date)]  
tY     <- dY[dates<=as.Date(test_date),]
tY[nrow(tY),] <- NA
tY[nrow(tY)-1,late_obs] <- NA  

#Likelihood based seasonal adjustment (for illustrative purposes). Alterntatively, use the package seasonal

N <- Predetermined.m(dates, predetermined = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))
ind  <- !is.na(tY[,'FRGSHPUSM649NCIS'])
y_SA <- SeasAdj_WE(tY[ind,'FRGSHPUSM649NCIS'], N[ind,], p = 2, Loud = F) # y - data, N - seasonal factors, p - lags in transition equation, Loud = T - print convergence of likelihood
tY[ind,'FRGSHPUSM649NCIS'] <- y_SA$y_SA

# ---- scaling and centering data is optional -------------

tY             <- 100*scale(tY)
scale_var      <- attr(tY,"scaled:scale") #save variance
scale_int      <- attr(tY,"scaled:center") #save intercept

# ---------- Estimate the Model -----------------
Est <- BDFM(tY, factors = 2, lags = 2, lam_B = 50, intercept = F)
y   <- scale_var[1]*Est$predicted[,1]/100 + scale_int[1] #get results and re-scale them
names(y) <- tDates
predictions <- tail(y,2)
print(test_date)
return(predictions)
}


OUT <- lapply(backtest_dates, function(x) try(run_nowcast(x))) #run over backtest dates. Wrap run_nowcast in try() so that if predictions fail at any date lapply keeps running.

#save(OUT, file = 'Backtest.Rdata')

true_vals <- get_fred_data(series_id = 'INDPRO', api_key = "9c39f8a4fbcc38c22475bc7f26988367", observation_start = "1980-01-15", observation_end = Sys.Date(), frequency = "m")
y_true    <- diff(log(true_vals$value)) #take log differences
names(y_true) <- true_vals$date[-1] #-1 due to diff

plot_data  <- matrix(NA,length(backtest_dates),2)
plot_dates <- rep(NA, length(backtest_dates))
for(j in 1:length(backtest_dates)){
  nc            <- OUT[[j]][1]
  true_val      <- y_true[names(y_true)%in%names(nc)]
  plot_data[j,] <- c(true_val,nc)
  plot_dates[j] <- names(nc)
}

plot_dates <- as.Date(plot_dates)
#Look at backtest results
plot(plot_dates,plot_data[,1], type = "l", col = "steelblue", lwd = 2)
lines(plot_dates,plot_data[,2], col = "red", lwd = 2, lty = 2)


