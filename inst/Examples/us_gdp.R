library(bdfm)

# Specifying when to take logs and differences is optional. However, user oversight is strongly recomended!
# To see which variables will be differenced or log differenced use:
# log_diffs <- should_log_diff(econ_us)
# Altermatively, specify logs and diffs manually:
logs <- c(
  "W068RCQ027SBEA",
  "PCEDG",
  "PCEND",
  "JTSJOL",
  "INDPRO",
  "CSUSHPINSA",
  "HSN1F",
  "TSIFRGHT",
  "IPG2211S",
  "DGORDER",
  "AMTMNO",
  "CPILFESL",
  "ICSA"
)

diffs <- setdiff(colnames(econ_us), c("A191RL1Q225SBEA", "USSLIND"))

# Forecasts should ALWAYS be made using store_idx if we are interested in forcasting
# one of the series in the model.
m <- dfm(data = econ_us, logs = logs, diffs = diffs, factors = 3, pre_differenced = "A191RL1Q225SBEA", store_idx = "A191RL1Q225SBEA")

# Are we drawing from a stationary distribution?
ts.plot(m$Qstore[1,1,])
ts.plot(m$Hstore[1,1,])

print(m)

# how did observed variables contribute to the nowcast update in January 2018?
window(m$idx_update, start = c(2018, 12), end = c(2018, 12))

# Fill in missing values in the econ_us data set using estimated values
# Note that this is already done for differenced series as level observations are not replaced by estimated values
Y_fill <- econ_us
Y_fill[!is.finite(Y_fill)] <- m$values[1:NROW(Y_fill), ][!is.finite(Y_fill)]


# Comparison of NA filling by substitution (above) and (adjustement)
# [does not work at the moment, as predict(m) is off.]
x <- ts_pick(econ_us, "JTSJOL")
x_subst <- ts_pick(Y_fill, "JTSJOL")
x_adj <- ts_pick(m$adjusted, "JTSJOL")

ts_c(x, x_subst, x_adj)




