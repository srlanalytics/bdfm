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

# Forecasts should ALWAYS be made using keep_posterior if we are interested in forcasting
# one of the series in the model.
m <- dfm(data = econ_us, logs = logs, diffs = diffs, factors = 3, pre_differenced = "A191RL1Q225SBEA", keep_posterior = "A191RL1Q225SBEA")

# Are we drawing from a stationary distribution?
ts.plot(m$Qstore[1,1,])
ts.plot(m$Hstore[1,1,])

# how did observed variables contribute to the nowcast update in January 2018?
window(m$idx_update, start = c(2018, 12), end = c(2018, 12))

predict(m)
adjusted(m)
factors(m)


