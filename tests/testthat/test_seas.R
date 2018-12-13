library(testthat)
library(bdfm)

context("seas_we")

test_that("seas_we works on minimal example", {
  m <- seas_we(mdeaths)
  adj_we <- predict(m)
  expect_is(adj_we, "ts")
})
