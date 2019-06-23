library(testthat)
library(bdfm)

context("PCdfm")

test_that("alternative identification works", {
  # works
  m <- dfm(cbind(mdeaths, fdeaths), identification = "pc_long", method = "pc")
  m <- dfm(cbind(mdeaths, fdeaths), identification = "pc_wide", method = "pc")

  # should work?
  m <- dfm(cbind(mdeaths, fdeaths), identification = "mdeaths", method = "pc")
  m <- dfm(cbind(mdeaths, fdeaths), identification = 1, method = "pc")

})


test_that("orthogonal_shocks works", {
  # works
  m <- dfm(cbind(mdeaths, fdeaths), orthogonal_shocks = TRUE, method = "pc")
  expect_is(m, "dfm")
})

