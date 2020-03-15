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
    expect_equal(binssim(x,x, method = "Jaccard"), 1.0)
    expect_equal(binssim(x,x, method = "rand"), 1.0)
    expect_equal(binssim(x,x), 1.0)
    expect_warning(catmssim_2d(x,x))
    
})

test_that("Gini works",{
    x = c(1,1)
    y = c(1,0)
    expect_equal(gini(x),0)
    expect_equal(gini(y),.5)
    expect_equal(ginicorr(y,2),1)
    expect_equal(ginicorr(x,2),0)
    expect_equal(ginicorr(x,1),1)
    expect_equal(sqrtgini(x),0)
    expect_equal(sqrtgini(y),1-sqrt(.5))
    expect_equal(sqrtginicorr(y,2),1)
    expect_equal(sqrtginicorr(x,2),0)
    expect_equal(sqrtginicorr(x,1),1)
    expect_equal(catsim:::jaccard(x,y),.5)
    expect_equal(ginicorr(1:15,15),1)
    expect_error(ginicorr(y,1))
    expect_error(sqrtginicorr(sample(letters,size = 100, replace = TRUE), 5))
    })



test_that("dimensions 2D work",{
    expect_error(catmssim_2d(x,y[,1:10]))
    expect_error(binssim(x,y[,1:10]))
    expect_error(AdjRandIndex(x,y[,1:10]))
    expect_error(catmssim_2d(x,y[1:10,]))
    expect_warning(catmssim_2d(x[1:2,],y[1:2,], weights = 1, method = "rand"))
    expect_error(catmssim_2d(x[1:3,],y[1:2,], weights = 1, method = "rand"))
    expect_error(catmssim_2d(x[1:2,],y[1,], weights = 1, method = "rand"))
    expect_error(catmssim_2d(y[1,],x[1,], weights = 1, method = "rand"))
    expect_warning(catmssim_2d(y[1,,drop=FALSE],x[1,,drop=FALSE], weights = 1, method = "jaccard"))
    expect_warning(AdjRandIndex(y[1,,drop=FALSE],x[1,,drop=FALSE]))
    expect_warning(catmssim_2d(x[,1:8],y[,1:8], window = 2, weights = c(.5,.25,.25,.25)))
    expect_error(catsim:::jaccard(x,y[,1:10]))
})

test_that("3D is not 2D", {
    expect_error(catmssim_3d_slice(x,y, weights = 1))
    expect_error(catmssim_3d_cube(x,y, weights = 1))
})

test_that("Inputs are symmetric 2D", {
    expect_equal(catmssim_2d(x,y, weights = 1), catmssim_2d(y,x, weights = 1))
      expect_equal(catmssim_2d(x,y, weights = c(.5,.5), method = "j"), catmssim_2d(y,x, weights = c(.5,.5), method = "j"))
    expect_equal(catmssim_2d(x,y, weights = 1, method = "rand"), catmssim_2d(y,x, weights = 1, method = "rand"))
  })

  test_that("Inputs are symmetric 3D",{  
 set.seed(20181207)
 dim = 16
 x <- array(sample(1:4, dim^3, replace = TRUE), dim = c(dim,dim,dim))
 y <- x
 for (j in 1:dim){
     for (i in 1:dim) y[i, i, j] = 1
     for (i in 1:(dim-1)) y[i, i+1, j] = 1
 }
 
 expect_equal(catmssim_3d_slice(x,y, weights = c(.75)), catmssim_3d_slice(y,x, weights = c(.75)))
 expect_equal(catmssim_3d_cube(x,y, weights = c(.75,.25), method = "j"), catmssim_3d_cube(y,x, weights = c(.75,.25), method = "j"))
 expect_warning(catmssim_3d_slice(x,y))
 expect_warning(catmssim_3d_cube(x,y))
 expect_error(catmssim_3d_slice(x,y[,,1:2]))
  expect_error(catmssim_3d_cube(x,y[,,1:2]))

  })


test_that("dimensions 3D work", {
    set.seed(20181207)
    dim = 16
    x <- array(sample(1:4, dim^3, replace = TRUE), dim = c(dim,dim,dim))
    y <- x
    for (j in 1:dim){
        for (i in 1:dim) y[i, i, j] = 1
        for (i in 1:(dim-1)) y[i, i+1, j] = 1
    }
    expect_error(catmssim_3d_slice(x,y[1:10,,]))
    expect_error(catmssim_3d_slice(x,y[,1:10,]))
    expect_error(catmssim_3d_cube(x,y[,1:10,]))
    expect_error(catmssim_3d_cube(x,y[1:10,,]))
    expect_error(catmssim_3d_slice(x,y[1,,]))
    expect_error(catmssim_3d_slice(x,y[,1,]))
    expect_error(catmssim_3d_cube(x,y[,1,]))
    expect_error(catmssim_3d_cube(x,y[1,,]))
    expect_error(catmssim_3d_cube(x[1,,],y[1,,]))
    expect_error(catmssim_3d_slice(x[1,,],y[1,,]))
    expect_error(catmssim_3d_slice(x[,1,1],y[,1,1]))
    expect_error(catmssim_3d_slice(x[,,1],y[,1,1]))
    expect_error(catmssim_3d_cube(x[,1,1],y[,1,1]))
    expect_error(catmssim_3d_cube(x[,,1],y[,1,1]))

    expect_error(catmssim_3d_slice(x[,1:3,],y[,1:3,]))
    expect_warning(catmssim_3d_cube(x[,1:3,],y[,1:3,]))

    expect_warning(catmssim_3d_cube(x[,1:3,],y[,1:3,], window = 3))

   # expect_warning(catmssim_3d_slice(x[,1:3,],y[,1:3,], window = 3, weight = 1))
    expect_warning(catmssim_3d_cube(x[,1:3,],y[,1:3,], window = 3, weight = 1))

#    expect_error(catmssim_3d_slice(x[,1:6,],y[,1:6,], window = 3, weight = 1))
#    expect_error(catmssim_3d_cube(x[,1:6,],y[,1:6,], window = 3, weight = 1))

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
