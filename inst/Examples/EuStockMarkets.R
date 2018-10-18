library(tsbox)
library(BDFM)
library(Matrix)

# data(EuStockMarkets) #load data
# Y     <- as.matrix(EuStockMarkets) #data in tabular (matrix) format
# dates <- unique(as.Date(ts_data.frame(EuStockMarkets)$time))
#log difference data

#Enter data as ts
Y <- EuStockMarkets
dY    <- 100*diff(log(Y))

#Estimate ML factor model with:
#m = 1 factor
#p = 5 lags
#Est   <- mlDFM(dY,factors = 1,lags = 5, Loud = T)

#Estimate Bayesian factor model with:
#m = 1 factor
#p = 5 lags
Est   <- bdfm(dY,factors = 1,lags = 5)

#Look at the common factor for the last 20 observations
ts.plot(cbind(tail(Est$factors,50), tail(dY,50)), col = c("red", "steelblue", "steelblue", "steelblue", "steelblue" ), lty = c(1,2,2,2,2))

#Compare the common factor to the first principal component
pc <- prcomp(dY, center = F, scale. = F)
ts.plot(cbind(tail(Est$factors,50), tail(pc$x[,1],50)), col = c("red", "steelblue"), lty = c(1,2))

#Re-estimate the factor model shrinking the variance of shocks in the transition equation towards zero (i.e. smooth factors)
Est2  <- bdfm(dY, factors = 1, lags = 5, nu_q = 10)
ts.plot(cbind(tail(Est2$factors,50), tail(pc$x[,1],50)), col = c("red", "steelblue"), lty = c(1,2))

Est2$B
Est2$H
Est2$R
Est2$q

# Est3  <- bdfm(dY, factors = 1, lags = 5, ID = "name", nu_q = 0)
# ts.plot(cbind(tail(Est3$factors,50), tail(dY[,1],50)), col = c("red", "steelblue"), lty = c(1,2))

#Enter data as data.frame
dates  <- unique(as.Date(ts_data.frame(EuStockMarkets)$time))
Y      <- as.matrix(EuStockMarkets)
dY     <- data.frame(dates[-1], diff(log(Y)))

#Estimate factor model with:
#m = 1 factor
#p = 5 lags
Est   <- BDFM(dY,factors = 1,lags = 5, forecast = 5)

tail(Est$predicted)

