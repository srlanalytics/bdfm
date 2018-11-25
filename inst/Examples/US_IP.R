
library(BDFM)
library(readr)
library(tsbox)

#Helper function

unscale <- function(x, attrib){
  x*(rep(1,nrow(attrib))%x%t(attr(attrib, "scaled:scale"))) + 
    rep(1,nrow(attrib))%x%t(attr(attrib, "scaled:center"))
}

#Import data: data is industrial production, manufacturer's new orders, wholsalers' inventory:sales ratio, manufacturers' inventory:sales ratio, initial jobles claims, all seasonally adjusted at a monthly frequency.
fred <- read_csv(system.file("Examples/freddata.csv", package = "BDFM"),
                  col_types = cols(DATE = col_date(format = "%Y-%m-%d"))
)

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

