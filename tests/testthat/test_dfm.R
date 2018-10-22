library(testthat)
library(BDFM)

context("interface functions")

test_that("dfm works on minimal example", {
  library(tsbox)
  fdeaths0 <- fdeaths
  fdeaths0[length(fdeaths0)] <- NA
  dta <- cbind(fdeaths0, mdeaths)

  library(BDFM)
  m <- dfm(dta, forecast = 2)
  expect_is(m, "dfm")
  a <- predict(m)
  expect_is(a, "ts")
})
