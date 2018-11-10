# Spain_IP <- dt[,1:2]
# colnames(Spain_IP) <- c("date", "IP")
# write.csv(Spain_IP, file = "C:/Users/Seth/Documents/io/toolbox/Book/Data/Spain_IP.csv", row.names = F)

#-------- Source files -----------
library(Rcpp)
library(Matrix)
library(readr)
library(BDFM)

# ------- load data --------------
Spain <- read_csv(system.file("Examples/Spain_IP.csv", package = "BDFM"),
  col_types = cols(date = col_date(format = "%Y-%m-%d"))
)
y <- Spain$IP
dates <- Spain$date
# --------------------------------

N <- Predetermined.m(dates, predetermined = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))

y_SA <- SeasAdj_WE(y, N, lags = 2, Loud = T) # y - data, N - seasonal factors, p - lags in transition equation, Loud = T - print convergence of likelihood

plot(as.Date(dates), y, type = "l", col = "red", xlab = "Year", ylab = "Industrial Production")
lines(as.Date(dates), y_SA$y_SA, col = "steelblue")
legend("topright",
  legend = c("Unadjusted IP", "Seasonally Adjusted IP"),
  lty = c(1, 1), col = c("red", "steelblue")
)
