context("test-output")


x <- matrix(sample(1:4, 64^2, replace = TRUE), nrow = 64)
y <- x
for (i in 1:64) y[i, i] = 1
for (i in 1:63) y[i, i + 1] = 1

test_that("Same input should have 1 as its result", {
  expect_equal(catmssim_2d(x,x), 1.0)
})

test_that("Zero weights should have result of 1.0", {
  expect_equal(catmssim_2d(x,y, weights = c(0.0)), 1.0)
})