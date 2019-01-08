

#-------- Source files -----------
# library(Rcpp)
# library(Matrix)
library(readr)
library(tsbox)
library(bdfm)
# ------- load data --------------
IP <- read_csv(system.file("Examples/IPB50001N.csv", package = "bdfm"),
  col_types = cols(date = col_date(format = "%Y-%m-%d"))
)

IP <- fread("C:/Users/Seth/Documents/GitHub/BDFM/inst/Examples/IPB50001N.csv")
IP[,1] <- as.Date(IP$DATE)

# --------------------------------

#----- Simple Example: Seasonally Adjusting US IP --------------

IP <- data.frame(IP)
y  <- ts_ts(IP)
sa <-seas_we(y, transformation = "log", verbose = T) 

plot(as.Date(IP$DATE), y, type = "l", col = "red", xlab = "Year", ylab = "Industrial Production")
lines(as.Date(IP$DATE), predict(sa), col = "steelblue")
legend("topleft",
  legend = c("Unadjusted IP", "Seasonally Adjusted IP"),
  lty = c(1, 1), col = c("red", "steelblue")
)

#------------------------------------------------------------------------

#----- Example 2: Getting Seasonal Adjustment Forecasts for 5 months ahead ----

sa <-seas_we(y, ahead = 5, transformation = "log", verbose = T)
ts.plot(tail(sa$factor,100))

#when using log transformation the SA factor is multiplicative (i.e. multiply forecasts)
#when using no transformation the SA factor is additive (i.e. add to forecasts)

#------------------------------------------------------------------------

#----- Example 3: Weekday Adjustments for Non-Stationary Daily Data ----

library(readr)
library(tsbox)

wil5000 <- read_csv(system.file("Examples/Wil5000.csv", package = "bdfm"),
                        col_types = cols(Date = col_date(format = "%Y-%m-%d"))
)

wil5000 <- data.frame(wil5000)
y <- ts_ts(wil5000)
test <- ts_data.frame(y)

sa <-seas_we(y, effect = "weekdays", holiday = NULL, transformation = "log", verbose = T)


#------------------------------------------------------------------------
library(reprex)
reprex()

library(readr)
library(bdfm)

wil5000 <- read_csv(system.file("Examples/Wil5000.csv", package = "bdfm"),
                    col_types = cols(Date = col_date(format = "%Y-%m-%d"))
)
sa <- seas_we(wil5000, effect = "weekdays", holiday = NULL, transformation = "log", verbose = T)
#> Loading required namespace: tsbox
#> [time]: 'Date' [value]: 'Value' 
#> [time]: 'Date' [value]: 'Value' 
#> [time]: 'Date' [value]: 'Value'

tsbox::ts_plot(wil5000, predict(sa))
#> [time]: 'Date' [value]: 'Value' 
#> [time]: 'Date' [value]: 'Value' 
#> [time]: 'Date' [value]: 'Value'
