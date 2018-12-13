library(testthat)
library(bdfm)

context("helper functions")

test_that("AnyNA works as expected", {
  expect_false(AnyNA(c(2, 3, 3)))
  expect_true(AnyNA(c(2, NA, 3)))
})
