

#-------- Source files -----------
# library(Rcpp)
# library(Matrix)
library(readr)
library(tsbox)
library(bdfm)
# ------- load data --------------
Spain <- read_csv(system.file("Examples/Spain_IP.csv", package = "bdfm"),
  col_types = cols(date = col_date(format = "%Y-%m-%d"))
)
y <- Spain$IP
dates <- Spain$date
# --------------------------------

#----- Simple Example: Seasonally Adjusting IP for Spain --------------

N <- Predetermined.m(dates, predetermined = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December", "Easter", "trading_days"))


# Note that data is scaled up by 100. This avoids small matrix determinents in the
# calculations which can make results inaccurate.
y_SA <- SeasAdj_WE(100*y, N, lags = 2, Loud = T) # 100*y - scaled up data, N - seasonal factors, p - lags in transition equation, Loud = T - print convergence of likelihood

plot(as.Date(dates), y, type = "l", col = "red", xlab = "Year", ylab = "Industrial Production")
lines(as.Date(dates), y_SA$y_SA/100, col = "steelblue")
legend("topright",
  legend = c("Unadjusted IP", "Seasonally Adjusted IP"),
  lty = c(1, 1), col = c("red", "steelblue")
)

#------------------------------------------------------------------------

#----- Example 2: Getting Seasonal Adjustment Forecasts for 5 months ahead ----

# Note that in this example monthly data is dated at the end of the month.
# For uniform frequency data this is not so important, but when dealing
# with mixed frequency data it is essential.

# To forecast seasonal factors just add dates to the vector of predetermined
# variables. N may be longer then y; the length of the output will correspond
# to the number of rows in N.

dates_long <- c(dates, End_This_Month(seq.Date(from = as.Date(paste(format(tail(dates,1), "%Y-%m"), "01", sep = "-")), length.out = 6, by = "month")[-1]) )

N <- Predetermined.m(dates, predetermined = c("January", "February", "March", "April", "May", "June",
"July", "August", "September", "October", "November", "December", "Easter", "trading_days"))

y_SA <- SeasAdj_WE(100*y, N, lags = 2, Loud = T) # 100*y - scaled up data, N - seasonal factors, p - lags in transition equation, Loud = T - print convergence of likelihood

plot(as.Date(dates), y_SA$sa/100, type = "l", col = "red", xlab = "Year", ylab = "Adjustment Factor")

#------------------------------------------------------------------------

#----- Example 3: Weekday Adjustments for Non-Stationary Daily Data ----
wil5000 <- read_csv(system.file("Examples/Wil5000.csv", package = "bdfm"),
                        col_types = cols(Date = col_date(format = "%Y-%m-%d"))
)
y         <- 10000*log(wil5000$Value)
dates     <- wil5000$Date
dtrend    <- loess(y ~ seq(1,length(y)), span = 0.5, na.action = na.exclude) #make input data stationary
N <- Predetermined.d(dates, predetermined = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday"))
SA <- SeasAdj_WE( y-predict(dtrend), N, lags = 5, Loud = T)
df <- cbind.data.frame(weekdays(tail(dates,5)), tail(SA$sa,5))
colnames(df) <- c("day of week", "effect")
print(df)
# Note that because we scaled the data up by 10000 these effects are very small
#------------------------------------------------------------------------
