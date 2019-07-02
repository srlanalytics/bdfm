library(testthat)
library(bdfm)

context("MLdfm")

test_that("alternative identification works", {
  # works
  m <- dfm(cbind(mdeaths, fdeaths), identification = "pc_long", method = "ml")
  m <- dfm(cbind(mdeaths, fdeaths), identification = "pc_wide", method = "ml")

  # should work?
  m <- dfm(cbind(mdeaths, fdeaths), identification = "mdeaths", method = "ml")
  m <- dfm(cbind(mdeaths, fdeaths), identification = 1, method = "ml")

})


test_that("orthogonal_shocks works", {
  m <- dfm(cbind(mdeaths, fdeaths), orthogonal_shocks = TRUE, method = "ml")
  expect_is(m, "dfm")
})

