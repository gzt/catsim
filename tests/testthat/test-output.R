###   test-output.R
###   CatSIM : Test output for CatSIM
###   Copyright (C) 2020  GZ Thompson <gzthompson@gmail.com>
###
###   This program is free software; you can redistribute it and/or modify
###   it under the terms of the GNU General Public License as published by
###   the Free Software Foundation; either version 3 of the License, or
###   (at your option) any later version.
###
###   This program is distributed in the hope that it will be useful,
###   but WITHOUT ANY WARRANTY; without even the implied warranty of
###   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###   GNU General Public License for more details.
###
###     You should have received a copy of the GNU General Public License
###   along with this program; if not, a copy is available at
###   https://www.R-project.org/Licenses/
context("test-output")

set.seed(20181215)
x <- matrix(sample(1:4, 32^2, replace = TRUE), nrow = 32)
y <- x
for (i in 1:32) y[i, i] <- 1
for (i in 1:31) y[i, i + 1] <- 1

test_that("Same input should have 1 as its result", {
  expect_equal(catmssim_2d(x, x, weights = 1), 1.0)
  expect_equal(catmssim_2d(x, x, weights = 1, method = "Rand"), 1.0)
  expect_equal(catmssim_2d(x, x, weights = 1, method = "hamming"), 1.0)
  expect_equal(catmssim_2d(x, x,
    weights = 1,
    method = "Accuracy"
  ), 1.0)
  expect_equal(catsim(x, x,
    weights = 1,
    method = "jaccard"
  ), 1.0)
  expect_equal(catsim(x, x, weights = 1, method = "dice"), 1.0)
  expect_equal(binssim(x, x, method = "Jaccard"), 1.0)
  expect_equal(binssim(x, x, method = "adjrand"), 1.0)
  expect_equal(binssim(x, x), 1.0)
  expect_warning(catmssim_2d(x, x))
})

test_that("Bad dimensions fail", {
  expect_equal(catmssim_2d(x, x, weights = 1L), 1.0)
  expect_error(catmssim_2d(x, c(x), weights = 1, random = NULL))
  expect_error(catmssim_2d(c(x), (x), weights = 1))
  expect_error(catsim(x, c(x), weights = 1))
  expect_error(catsim(c(x), (x), weights = 1))
  expect_error(catmssim_3d_slice(x, (x), weights = 1))
  expect_error(catmssim_3d_cube(x, (x), weights = 1))
  expect_error(catsim(x, x, method = "bob", weights = 1))
})

test_that("Weights and levels work 2D", {
  expect_error(catsim(x, y, weights = 1:3, levels = 5))
  expect_equal(
    {
      set.seed(20200408)
      catsim(x, y, random = "random", weights = c(.5, .5))
    },
    {
      set.seed(20200408)
      catsim(x, y, random = "random", levels = 2)
    }
  )
  expect_equal(
    catsim(x, y, random = "pseudo", weights = c(.5, .5, .5), levels = 2),
    catsim(x, y, random = "pseudo", weights = c(.5, .5))
  )
  expect_equal(
    catsim(x, y, window = 2, random = NULL),
    catsim(x, y, window = 2, weights = rep(1, 5) / 5, random = NULL)
  )
  expect_equal(
    catsim(x, y, window = 2, random = NULL),
    catsim(x, y, window = 2, levels = 5, random = NULL)
  )
})

test_that("Gini works", {
  x <- c(1, 1)
  y <- c(1, 0)
  expect_equal(gini(x), 0)
  expect_equal(gini(y), .5)
  expect_equal(ginicorr(y, 2), 1)
  expect_equal(ginicorr(x, 2), 0)
  expect_equal(ginicorr(x, 1), 1)
  expect_equal(sqrtgini(x), 0)
  expect_equal(sqrtgini(y), 1 - sqrt(.5))
  expect_equal(sqrtginicorr(y, 2), 1)
  expect_equal(sqrtginicorr(x, 2), 0)
  expect_equal(sqrtginicorr(x, 1), 1)
  expect_equal(catsim:::jaccard(x, y), .5)
  expect_equal(ginicorr(1:15, 15), 1)
  expect_error(ginicorr(y, 1))
  expect_error(sqrtginicorr(sample(letters, size = 100, replace = TRUE), 5))
})

set.seed(20181215)
x <- matrix(sample(1:4, 32^2, replace = TRUE), nrow = 32)
y <- x
for (i in 1:32) y[i, i] <- 1
for (i in 1:31) y[i, i + 1] <- 1

test_that("dimensions 2D work", {
  expect_error(catmssim_2d(x, y[, 1:10]))
  expect_error(binssim(x, y[, 1:10]))
  expect_error(adj_rand(x, y[, 1:10]))
  expect_error(rand_index(x, y[, 1:10]))
  expect_error(normalized_mi(x, y[, 1:10]))
  expect_error(adjusted_mi(x, y[, 1:10]))
  expect_error(cohen_kappa(x, y[, 1:10]))
  expect_error(catmssim_2d(x, y[1:10, ]))
  expect_warning(catmssim_2d(x[1:2, ], y[1:2, ], weights = 1, method = "rand"))
  expect_error(catmssim_2d(x[1:3, ], y[1:2, ], weights = 1, method = "jaccard"))
  expect_error(catsim(x[1:2, ], y[1, ], weights = 1, method = "adjrand"))
  expect_error(catmssim_2d(y[1, ], x[1, ], weights = 1, method = "accuracy"))
  expect_warning(catmssim_2d(y[1, , drop = FALSE],
    x[1, , drop = FALSE],
    weights = 1, method = "jaccard"
  ))
  expect_warning(catmssim_2d(x[, 1:8], y[, 1:8],
    window = 2, weights = c(.5, .25, .25, .25)
  ))
  expect_error(catsim:::jaccard(x, y[, 1:10]))
})

test_that("3D is not 2D", {
  expect_error(catmssim_3d_slice(x, y, weights = 1))
  expect_error(catmssim_3d_cube(x, y, weights = 1))
})

test_that("Inputs are symmetric 2D", {
  expect_equal(catmssim_2d(x, y, weights = 1, random = NULL), catmssim_2d(y, x, weights = 1, random = NULL))
  expect_equal(
    catmssim_2d(x, y, weights = 1, method = "NMI"),
    catmssim_2d(y, x, weights = 1, method = "NMI")
  )
  expect_equal(
    catmssim_2d(x, y, weights = 1, method = "AMI"),
    catmssim_2d(y, x, weights = 1, method = "AMI")
  )
  expect_equal(
    catmssim_2d(x, y, weights = c(.5, .5), method = "j", random = "pseudo"),
    catmssim_2d(y, x, weights = c(.5, .5), method = "j", random = "pseudo")
  )
  expect_equal(
    catmssim_2d(x, y, weights = 1, method = "dice"),
    catmssim_2d(y, x, weights = 1, method = "dice")
  )
  expect_equal(normalized_mi(x, y), normalized_mi(y, x))
  expect_equal(adjusted_mi(x, y), adjusted_mi(y, x))
})

test_that("catmssim_2d equivalent to catsim in 2D", {
  expect_equal(catmssim_2d(x, y, weights = 1), catsim(x, y, weights = 1))
  expect_equal(
    catmssim_2d(x, y, weights = c(.23, 77), random = NULL),
    catsim(x, y, weights = c(.23, 77), random = NULL)
  )
  expect_equal(
    catmssim_2d(x, y, window = 5, weights = 1),
    catsim(x, y, window = 5, weights = 1)
  )
  expect_equal(
    catmssim_2d(x, y, window = c(5, 6), weights = 1),
    catsim(x, y, window = c(5, 6), weights = 1)
  )
  expect_equal(
    catmssim_2d(x, y, window = c(5), weights = 1),
    catsim(x, y, window = c(5, 5), weights = 1)
  )
})

test_that("Inputs are symmetric 3D", {
  set.seed(20181207)
  dim <- 16
  x <- array(sample(1:4, dim^3, replace = TRUE), dim = c(dim, dim, dim))
  y <- x
  for (j in 1:dim) {
    for (i in 1:dim) y[i, i, j] <- 1
    for (i in 1:(dim - 1)) y[i, i + 1, j] <- 1
  }
  expect_equal(catmssim_3d_slice(x, x,
    weights = c(.75, .25),
    window = 2, random = "pseudo"
  ), 1.0)
  expect_equal(catmssim_3d_cube(x, x,
    weights = c(.75, .25),
    window = 2, random = NULL
  ), 1.0)
  expect_equal(
    catmssim_3d_slice(x, y, weights = c(.75, .25), window = 2, random = NULL),
    catmssim_3d_slice(y, x, weights = c(.75, .25), window = 2, random = NULL)
  )

  expect_equal(
    catmssim_3d_slice(x, y, weights = c(.75, .25), window = 2, random = NULL),
    catsim(y, x, weights = c(.75, .25), window = 2, cube = FALSE, random = NULL)
  )

  expect_equal(
    catmssim_3d_cube(x, y, weights = 1),
    catsim(x, y, weights = 1)
  )
  expect_equal(
      catmssim_3d_cube(x, y, weights = c(.75, .25), method = "dice", random = NULL),
      catmssim_3d_cube(y, x, weights = c(.75, .25), method = "dice", random = NULL)
  )
  expect_warning(catmssim_3d_slice(x, y, random = NULL))
  expect_warning(catmssim_3d_cube(x, y, random = NULL))
  expect_error(catmssim_3d_slice(x, y[, , 1:2], random = NULL))
  expect_error(catmssim_3d_cube(x, y[, , 1:2], random = NULL))
})


test_that("dimensions 3D work", {
  set.seed(20181207)
  dim <- 16
  x <- array(sample(1:4, dim^3, replace = TRUE), dim = c(dim, dim, dim))
  y <- x
  for (j in 1:dim) {
    for (i in 1:dim) y[i, i, j] <- 1
    for (i in 1:(dim - 1)) y[i, i + 1, j] <- 1
  }
  expect_equal(catsim(x, x,
    weights = c(.5, .5),
    window = c(5, 5, 5), method = "accuracy", random = NULL
  ), 1.0)
  expect_equal(catsim(x, x,
    weights = c(.5, .5),
    window = c(5, 5, 5), method = "NMI", random = NULL
  ), 1.0)
  expect_error(catsim(x, x,
    weights = c(.5, .5),
    window = c(5, 5), method = "NMI", random = NULL
  ))
  expect_equal(catsim(x, x,
    weights = c(.5, .5),
    cube = FALSE, window = c(4, 4), random = NULL
  ), 1.0)
  expect_error(catsim(x, x[1:15, , ], weights = c(.5, .5)))
  expect_error(catmssim_3d_slice(x, y[1:10, , ]))
  expect_error(catmssim_3d_slice(x, y[, 1:10, ]))
  expect_error(catmssim_3d_cube(x, y[, 1:10, ]))
  expect_error(catmssim_3d_cube(x, y[1:10, , ]))
  expect_error(catmssim_3d_slice(x, y[1, , ]))
  expect_error(catmssim_3d_slice(x, y[, 1, ]))
  expect_error(catmssim_3d_cube(x, y[, 1, ]))
  expect_error(catmssim_3d_cube(x, y[1, , ]))
  expect_error(catmssim_3d_cube(x[1, , ], y[1, , ]))
  expect_error(catmssim_3d_slice(x[1, , ], y[1, , ]))
  expect_error(catmssim_3d_slice(x[, 1, 1], y[, 1, 1]))
  expect_error(catmssim_3d_slice(x[, , 1], y[, 1, 1]))
  expect_error(catmssim_3d_cube(x[, 1, 1], y[, 1, 1]))
  expect_error(catmssim_3d_cube(x[, , 1], y[, 1, 1]))
  expect_error(catmssim_3d_slice(x[, 1:3, ], y[, 1:3, ]))
  expect_warning(catmssim_3d_cube(x[, 1:3, ], y[, 1:3, ]))
  expect_error(catmssim_2d(x, y[1, , ]))
  expect_error(catsim(x, y[1, , ]))
  expect_warning(catmssim_3d_cube(x[, 1:3, ], y[, 1:3, ], window = 3))

  expect_warning(catmssim_3d_cube(x[, 1:3, ],
    y[, 1:3, ],
    window = 3, weight = 1
  ))
})

set.seed(20181215)
x <- matrix(sample(1:4, 32^2, replace = TRUE), nrow = 32)
y <- x
for (i in 1:32) y[i, i] <- 1
for (i in 1:31) y[i, i + 1] <- 1

test_that("Zero weights should have result of 1.0", {
  expect_equal(catmssim_2d(x, y, weights = c(0, 0), random = "random"), 1.0)
  expect_equal(binssim(x, x, sqrtgini = TRUE), 1.0)
  expect_equal(catsim(x, y, weights = c(0.0)), 1.0)
})

test_that("NA RM works", {
  x <- c(rep(0:5, 5), NA)
  y <- c(rep(0:5, 4), rep(0, 6), 3)
  xtrim <- x[1:30]
  ytrim <- y[1:30]
  expect_warning(rand_index(x, y))
  expect_warning(adj_rand(x, y))
  expect_warning(cohen_kappa(x, y))
  expect_warning(normalized_mi(x, y))
  expect_warning(adjusted_mi(x, y))
  expect_equal(rand_index(x, y, na.rm = TRUE), rand_index(xtrim, ytrim))
  expect_equal(adj_rand(x, y, na.rm = TRUE), adj_rand(xtrim, ytrim))
  expect_equal(cohen_kappa(x, y, na.rm = TRUE), cohen_kappa(xtrim, ytrim))
  expect_equal(normalized_mi(x, y, na.rm = TRUE), normalized_mi(xtrim, ytrim))
  expect_equal(adjusted_mi(x, y, na.rm = TRUE), adjusted_mi(xtrim, ytrim))
})


test_that("jaccard works with NA", {
  expect_equal(catsim(hoffmanphantom[, , 1], hoffmanphantom[, , 1], random = "pseudo", levels = 2), 1)
})
