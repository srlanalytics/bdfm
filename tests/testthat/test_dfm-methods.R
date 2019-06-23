library(testthat)
library(bdfm)

context("output methods")

test_that("output methods work", {
  # works
  m <- dfm(cbind(mdeaths, fdeaths))
  expect_output(print(m))
  expect_output(summary(m))
  expect_is(factors(m), "ts")
  expect_is(predict(m), "ts")
  expect_is(adjusted(m), "ts")
})

