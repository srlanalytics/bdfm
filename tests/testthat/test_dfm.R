library(testthat)
library(bdfm)

context("interface functions")

test_that("dfm works on minimal example", {
  library(tsbox)
  fdeaths0 <- fdeaths
  fdeaths0[length(fdeaths0)] <- NA
  dta <- cbind(fdeaths0, mdeaths)

  library(bdfm)
  m0 <- dfm(dta, forecast = 2)
  expect_is(m0, "dfm")
  a0 <- predict(m0)
  expect_is(a0, "ts")

  m1 <- dfm(fdeaths, forecast = 2)
  expect_is(m1, "dfm")
  a1 <- predict(m1)
  expect_is(a1, "ts")


})

test_that("dfm works with ml, pc method", {
  # https://github.com/srlanalytics/bdfm/issues/38
  library(bdfm)
  m0 <- dfm(cbind(fdeaths, mdeaths), method = "pc") #multivariate example
  m1 <- dfm(cbind(fdeaths, mdeaths), method = "ml") #multivariate example
  m2 <- dfm(fdeaths, method = "pc") #univariate example
  m3 <- dfm(fdeaths, method = "ml") #univariate example

  expect_is(predict(m0), "ts")
  expect_is(predict(m1), "ts")
  expect_is(predict(m2), "ts")
  expect_is(predict(m3), "ts")

})




test_that("forecast > 0 works", {
  dfm(cbind(mdeaths, fdeaths), factors = 2, lags = 3, forecast = 3)
  dfm(mdeaths, factors = 1, lags = 3)
})


test_that("tsbox is supported", {
  library(tsbox)
  data_df <- ts_data.frame(cbind(mdeaths, fdeaths))
  m <- dfm(data_df, forecast = 3)
  expect_is(predict(m), "data.frame")
})


test_that("manual logs, diffs specification works", {
  m <- dfm(cbind(fdeaths, mdeaths), logs = "fdeaths", diffs = "fdeaths")
  expect_is(predict(m), "ts")
})


test_that("prior setting for obs_df works", {
  m <- dfm(cbind(fdeaths, mdeaths), obs_df = c("fdeaths" = 1))
  expect_is(predict(m), "ts")
})


test_that("non-stationary series are differenced", {
  m <- dfm(econ_us[,c(4, 5)])
  expect_is(predict(m), "ts")
})


