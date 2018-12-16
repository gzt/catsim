context("test-output")

set.seed(20181215)
x <- matrix(sample(1:4, 32^2, replace = TRUE), nrow = 32)
y <- x
for (i in 1:32) y[i, i] = 1
for (i in 1:31) y[i, i + 1] = 1

test_that("Same input should have 1 as its result", {
  expect_equal(catmssim_2d(x,x, weights = 1), 1.0)
})

test_that("Zero weights should have result of 1.0", {
  expect_equal(catmssim_2d(x,y, weights = c(0.0)), 1.0)
  expect_equal(binssim(x,x, sqrtgini = TRUE), 1.0)
})

x <- matrix(c(0,1), nrow = 512, ncol = 512)
y <- matrix(c(0,0,1,1), nrow = 512, ncol = 512)
test_that("Large matrices should work", {
  expect_true(AdjRandIndex(x,y)$Cohen < .01)
  expect_true(AdjRandIndex(x,y)$AdjRand < .01)
})