library(jsonlite)
library(data.table)
library(BDFM)

setwd("...")
source("./Nowcast_IP/helper_functions.R")

unique_names <- c('INDPRO', #1 Industrial Production Index, monthly, level
                  'USSLIND', #2 US Leading Index --- this will be adjusted forward one period
                  'UNRATE', #3 Unemployment, monthly
                  'FRGSHPUSM649NCIS', #4 CASS freight index, level, not SA --- comes out one day ahead (sometimes)
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

#get a FRED api key at https://research.stlouisfed.org/docs/api/api_key.html
api_key = "XXX"

Get_Data <- function(series_name){
    Data <- get_fred_data(series_id = series_name, api_key = api_key, observation_start = "1980-01-01", frequency = "m")
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

#Likelihood based seasonal adjustment (for illustrative purposes). Alterntatively, use the package seasonal

N <- Predetermined.m(dates, predetermined = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))
ind  <- !is.na(dY[,'FRGSHPUSM649NCIS'])
y_SA <- SeasAdj_WE(dY[ind,'FRGSHPUSM649NCIS'], N[ind,], p = 2, Loud = F) # y - data, N - seasonal factors, p - lags in transition equation, Loud = T - print convergence of likelihood
dY[ind,'FRGSHPUSM649NCIS'] <- y_SA$y_SA

# ---- scaling and centering data is optional -------------

dY             <- 100*scale(dY)
scale_var      <- attr(dY,"scaled:scale") #save variance
scale_int      <- attr(dY,"scaled:center") #save intercept

# ---------- Estimate the Model -----------------
fc  <- 1
Est <- BDFM(dY, factors = 2, lags = 2, lam_B = 50, forecast = fc, intercept = F)
dates <- c(dates,seq.Date(from = tail(dates,1), length.out = fc+1, by = "m")[-1])

#Look at fitted values
plot(dates,c(dY[,'INDPRO'], rep(NA,fc)), type = "l", col = "steelblue", lwd = 2)
lines(dates,Est$predicted[,1], col = "red", lwd = 2, lty = 2)


