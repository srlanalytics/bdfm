library(bdfm)

# CS: I really would like to have helpers for this. Two functions: should_log(),
# should_diff(), that return a vector with the series that require the
# treatment.


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
m <- dfm(data = econ_us, factors = 3, pre_differenced = "A191RL1Q225SBEA", store_idx = "A191RL1Q225SBEA",
         logs = logs, diffs = diffs)

# Are we drawing from a stationary distribution?
ts.plot(m$Qstore[1,1,])
ts.plot(m$Hstore[1,1,])

summary(m)

# how did observed variables contribute to the nowcast update in January 2018?
window(m$idx_update, start = c(2018, 1), end = c(2018, 1))

