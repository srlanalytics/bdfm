
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




