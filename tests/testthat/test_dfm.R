library(testthat)
library(bdfm)

context("interface functions")

test_that("dfm works on minimal example", {
  library(tsbox)
  fdeaths0 <- fdeaths
  fdeaths0[length(fdeaths0)] <- NA
  dta <- cbind(fdeaths0, mdeaths)

  library(bdfm)
  m <- dfm(dta, forecast = 2)
  expect_is(m, "dfm")
  a <- predict(m)
  expect_is(a, "ts")
})


# test_that("dfm works with ml, pc method", {
#   # https://github.com/srlanalytics/bdfm/issues/38
#   library(bdfm)
#   m0 <- dfm(cbind(fdeaths, mdeaths), method = "pc")
#   #> Error in dfm(cbind(fdeaths, mdeaths), method = "pc"): Mixed freqeuncy models are only supported for Bayesian estimation
#   m1 <- dfm(fdeaths, method = "pc")
#   #> Error in dfm(fdeaths, method = "pc"): Mixed freqeuncy models are only supported for Bayesian estimation
#   m2 <- dfm(fdeaths, method = "ml")

#   expect_is(predict(m0), "ts")
#   expect_is(predict(m1), "ts")
#   expect_is(predict(m2), "ts")

# })

