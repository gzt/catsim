context("test-output")

set.seed(20181215)
x <- matrix(sample(1:4, 32^2, replace = TRUE), nrow = 32)
y <- x
for (i in 1:32) y[i, i] = 1
for (i in 1:31) y[i, i + 1] = 1

test_that("Same input should have 1 as its result", {
    expect_equal(catmssim_2d(x,x, weights = 1), 1.0)
    expect_equal(catmssim_2d(x,x, weights = 1, method = "Rand"), 1.0)
    expect_equal(catmssim_2d(x,x, weights = 1, method = "Jaccard"), 1.0)
    expect_warning(catmssim_2d(x,x))
    
})

test_that("Inputs are symmetric", {
    expect_equal(catmssim_2d(x,y, weights = 1), catmssim_2d(y,x, weights = 1))
      expect_equal(catmssim_2d(x,y, weights = c(.5,.5), method = "j"), catmssim_2d(y,x, weights = c(.5,.5), method = "j"))
    expect_equal(catmssim_2d(x,y, weights = 1, method = "rand"), catmssim_2d(y,x, weights = 1, method = "rand"))


    
 set.seed(20181207)
 dim = 16
 x <- array(sample(1:4, dim^3, replace = TRUE), dim = c(dim,dim,dim))
 y <- x
 for (j in 1:dim){
 for (i in 1:dim) y[i, i, j] = 1
 for (i in 1:(dim-1)) y[i, i+1, j] = 1
 }

    expect_equal(catmssim_3d_slice(x,y, weights = c(.75,.25)), catmssim_3d_slice(y,x, weights = c(.75,.25)))
    expect_equal(catmssim_3d_cube(x,y, weights = c(.75,.25), method = "j"), catmssim_3d_cube(y,x, weights = c(.75,.25), method = "j"))
    expect_warning(catmssim_3d_slice(x,y))
     expect_warning(catmssim_3d_cube(x,y))
})


set.seed(20181215)
x <- matrix(sample(1:4, 32^2, replace = TRUE), nrow = 32)
y <- x
for (i in 1:32) y[i, i] = 1
for (i in 1:31) y[i, i + 1] = 1

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
