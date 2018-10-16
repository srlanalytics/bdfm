library(tsbox)
library(BDFM)

data(EuStockMarkets) #load data
Y     <- as.matrix(EuStockMarkets) #data in tabular (matrix) format
dates <- unique(as.Date(ts_data.frame(EuStockMarkets)$time))
#log difference data
dY    <- diff(log(Y))
ts.plot(dY)

#Estimate factor model with:
#m = 1 factor
#p = 5 lags
Est   <- BDFM(dY,m = 1,p = 5)

#Look at the common factor for the last 20 observations
ts.plot(cbind(tail(Est$Z[,1],20), tail(dY,20)), col = c("red", "steelblue", "steelblue", "steelblue", "steelblue" ), lty = c(1,2,2,2,2))

#Compare the common factor to the first principal component
pc <- prcomp(dY, center = F, scale. = F)
ts.plot(cbind(tail(Est$Z[,1],20), tail(pc$x[,1],20)), col = c("red", "steelblue"), lty = c(1,2))

#Re-estimate the factor model shrinking the variance of shocks in the transition equation towards zero
Est2  <- BDFM(dY, m = 1, p = 5, nu_q = 1000)
ts.plot(cbind(tail(Est2$Z[,1],20), tail(pc$x[,1],20)), col = c("red", "steelblue"), lty = c(1,2))

