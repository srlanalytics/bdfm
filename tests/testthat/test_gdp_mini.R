library(testthat)
library(bdfm)

context("minimal GDP example")

test_that("dfm works on realistic example", {

  econ_us_mini <- econ_us[,1:4]

  logs <- c(
    "W068RCQ027SBEA",
    "PCEDG"
  )

  diffs <- setdiff(colnames(econ_us_mini), c("A191RL1Q225SBEA", "USSLIND"))

  # Forecasts should ALWAYS be made using store_idx if we are interested in forcasting
  # one of the series in the model.
  m <- dfm(data = econ_us_mini, factors = 1, pre_differenced = "A191RL1Q225SBEA", store_idx = "A191RL1Q225SBEA",
           logs = logs, diffs = diffs)

  summary(m)


})
