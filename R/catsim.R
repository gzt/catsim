###   catsim.R
###   CatSIM : CatSIM main functions
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

#' @useDynLib catsim, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


#' Means Function (internal)
#'
#' @param x,y binary or categorical image or vector
#' @param c1 small constant
#'
#' @return measure of similarity of means
#' @keywords internal
#'
#' @noRd
meansfunc <- function(x, y, c1 = 0.01) {
  c_meansfunc(x, y, c1)
}

#' @title Diversity Indices
#' @name gini
#'
#' @description `gini()` is a measure of diversity that goes by a
#' number of different names, such as the probability of interspecific encounter
#' or the Gibbs-Martin index. It is \eqn{1 - sum(p_i^2)}, where \eqn{p_i} is the
#' probability of observing class i.
#'
#' @param x binary or categorical image or vector
#'
#' @return The index (between 0 and 1), with 0 indicating no variation and 1
#'    being maximal. The Gini index is bounded above by \eqn{1-1/k} for a group
#'    with `k` categories. The modified index is bounded above by
#'    \eqn{1-1/\sqrt{k}}.  The corrected indices fix this by dividing by the
#'    maximum.
#' @export

#' @examples
#' x <- rep(c(1:4), 5)
#' gini(x)
gini <- function(x) {
  x <- as.numeric(x)
  c_gini(x)
}


#' @name Corrected Gini-Simpson index
#'
#' @description The corrected Gini-Simpson index, `ginicorr` takes the
#' index and  corrects it so that the maximum possible is 1. If there are
#' `k` categories, the maximum possible of the uncorrected index is
#' \eqn{1-1/k}. It corrects the index by dividing by the maximum.
#' `k` must be specified.
#'
#' @param k number of categories
#' @export
#' @rdname gini
#' @examples
#'
#' x <- rep(c(1:4), 5)
#' ginicorr(x, 4)
ginicorr <- function(x, k) {
  if (k < length(unique(x))) stop("more unique values in x than k")
  if (k > 1) {
    c_gini(x) / (1 - 1 / k)
  } else {
    1
  }
}

#' @name Modified Gini-Simpson index
#'
#' @description The modified Gini-Simpson index is similar to the unmodified,
#' except it uses the square root of the summed squared
#' probabilities, that is, \eqn{1 - \sqrt{ sum(p_i^2)}}, where \eqn{p_i} is the
#' probability of observing class i.
#'
#' @export
#' @rdname gini
#' @examples
#'
#' x <- rep(c(1:4), 5)
#' sqrtgini(x)
sqrtgini <- function(x) {
  1 - sqrt(1 - c_gini(x))
}

#' @name Modified Corrected Gini index
#'
#' @description The modified corrected Gini index then
#' corrects the modified index for the number of categories, `k`.
#'
#' @export
#' @rdname gini
#' @examples
#'
#' x <- rep(c(1:4), 5)
#' sqrtginicorr(x, 4)
sqrtginicorr <- function(x, k) {
  if (k < length(unique(x))) stop("more unique values in x than k")

  if (k > 1) {
    sqrtgini(x) / (1 - 1 / sqrt(k))
  } else {
    1
  }
}

#' Variance function (internal)
#'
#' @param x,y binary or categorical image or vector
#' @param c2 small constant
#' @noRd
#' @return variance function
#' @keywords internal
#' @examples
#'
#' x <- rep(1:4, 4)
#' y <- c(rep(1:4, 3), rep(4, 4))
#' cfunc(x, y, k = 4)
cfunc <- function(x, y, c2 = 0.01, k, sqrtgini = TRUE) {
  c_cfunc(x, y, c2, k, sqrtgini)
}

#' Jaccard Index
#'
#' @param x binary image or vector
#' @param y binary image or vector
#'
#' @return Jaccard index - note this privileges 1 over 0
#' @keywords internal
#' @noRd
#' @examples
#'
#' x <- rep(c(0, 1), 7)
#' y <- rep(c(1, 1), 7)
#' jaccard(x, y)
jaccard <- function(x, y) {
  if (length(x) != length(y)) stop("x and y have differing lengths.")
  sumxy <- sum(x | y, na.rm = TRUE)
  if (is.na(sumxy) || !sumxy) {
    return(NA)
  }

  jaccard <- sum(x & y, na.rm = TRUE) / sumxy
  jaccard
}

#' Dice Index
#'
#' @param x,y binary image or vector
#'
#' @return Dice index
#' @keywords internal
#' @noRd
#' @examples
#' x <- rep(c(0, 1), 7)
#' y <- rep(c(1, 1), 7)
#' dice(x, y)
dice <- function(x, y) {
  jacc <- jaccard(x, y)
  2 * jacc / (1 + jacc)
}

#' Accuracy/Hamming Index
#'
#' @param x,y binary image or vector
#'
#' @return Accuracy (Hamming index). Proportion of pixels (voxels) that agree.
#' @keywords internal
#' @noRd
#' @examples
#' x <- rep(c(0, 1), 7)
#' y <- rep(c(1, 1), 7)
#' hamming(x, y)
hamming <- function(x, y) {
  if (length(x) != length(y)) stop("x and y have differing lengths.")
  n <- sum(!is.na(x) | !is.na(y))
  a <- sum(x == y, na.rm = TRUE)
  a / n
}


#' Covariance function (internal)
#'
#' @param x,y binary or categorical image or vector
#' @noRd
#'
#' @return Covariance function
#' @keywords internal
#' @examples
#' x <- rep(1:4, 4)
#' y <- c(rep(1:4, 3), rep(4, 4))
#' sfunc(x, y)
sfunc <- function(x, y, methodflag = c_cohen) {
  x <- as.vector(x)
  y <- as.vector(y)
  methodflag(x, y)
}


#' Categorical Structural Similarity Index Metric (whole image)
#'
#' This computes the categorical or binary structural similarity index metric
#' on a whole-image scale. The difference between this and the default 2-D
#' method is that this considers the whole image at once and one scale rather
#' than computing the index over a sliding window and downsampling to consider
#' it at other scales.
#'
#' @param x,y binary or categorical image
#' @param alpha normalizing parameter, by default 1
#' @param beta normalizing parameter, by default 1
#' @param gamma normalizing parameter, by default 1
#' @param c1 small normalization constant for the `c` function,
#' by default 0.01
#' @param c2 small normalization constant for the `s` function,
#' by default 0.01
#' @param method whether to use Cohen's kappa (`Cohen`),
#'   Jaccard Index (`Jaccard`), Dice index (`Dice`),
#'   accuracy (`accuracy`),  Rand index (`Rand`),
#'   Adjusted Rand Index (`AdjRand` or `ARI`), or normalized mutual
#'   information (`NMI` or `MI`) as the similarity index.
#'   Note Jaccard and Dice should only be used on binary data.
#' @param ... Constants can be passed to the components of the index.
#'
#' @return Structural similarity index.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' x <- matrix(sample(1:4, 10000, replace = TRUE), nrow = 100)
#' y <- x
#' for (i in 1:100) y[i, i] <- 1
#' for (i in 1:99) y[i, i + 1] <- 1
#' binssim(x, y)
binssim <- function(x, y, alpha = 1, beta = 1, gamma = 1,
                    c1 = 0.01, c2 = 0.01, method = "Cohen", ...) {
  if (length(x) != length(y)) stop("x and y must be the same size.")
  naxy <- (!is.na(x) & !is.na(y))
  k <- length(unique(c(x[naxy], y[naxy])))
  dotlist <- dots_parser(...)
  methodflag <- method_parser(method)
  (meansfunc(x[naxy], y[naxy], c1)^alpha) *
    (cfunc(
      x = x[naxy], y = y[naxy], c2 = c2, k = k,
      sqrtgini = dotlist[["sqrtgini"]]
    )^beta) *
    (sfunc(x[naxy], y[naxy], methodflag)^gamma)
}

#' Categorical SSIM Components
#'
#' @param x,y binary or categorical image
#' @param method whether to use Cohen's kappa, Jaccard, or Adjusted Rand Index
#'   as the similarity index. Note Jaccard should only be used on binary data.
#' @param c1 constant for the means function
#' @param c2 constant for the chrominance function
#' @param sqrtgini whether to use sqrtgini or gini, by default `TRUE`
#' @noRd
#' @return the three components of the Categorical SSIM.
#' @keywords internal
#'
ssimcomponents <- function(x, y, k, method = c_cohen,
                           c1 = 0.01, c2 = 0.01, sqrtgini = TRUE) {
  naxy <- (!is.na(x) & !is.na(y))
  if (length(x[naxy]) == 0) {
    return(c(NA, NA, NA))
  }

  c(
    meansfunc(x[naxy], y[naxy], c1),
    cfunc(x = x[naxy], y = y[naxy], c2 = c2, k = k, sqrtgini),
    sfunc(x[naxy], y[naxy], methodflag = method)
  )
}


#' Picks the first mode
#'
#' @param x binary or categorical image
#' @param rand if `random`, samples using [sample()]. if `pseudo`,
#'             cycles through options depending on an augmented
#'             global variable (ie a simple PRNG)
#'             If \code{NULL},
#'             will choose the first mode.
#' @param modepick changes each time you call the function and need to break
#'             a tie if rand = FALSE using a simple PRNG.
#'
#' @return The first mode in the default order returned by R and the state of
#'         the PRNG is applicable
#' @noRd
#' @keywords internal
#'
pickmode <- function(x, rand = NULL, modepick = 1) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  if (is.null(rand)) rand <- "none"
  y <- seq_along(tab)[tab == max(tab)]
  if (length(y) > 1L) {
    if (rand == "pseudo") {
      modepick <- ((75 * modepick) %% 65537) + 1
      c(
        ux[(modepick %% length(y)) + 1],
        modepick
      )
    } else if (rand == "random") {
      c(ux[sample(y, 1L)], modepick)
    } else {
      c(ux[y], modepick)
    }
  } else {
    c(ux[y], modepick)
  }
}

#' Downsampling by a factor of 2 for a 2D categorical image
#'
#' Cuts the image into a grid of \eqn{2 \times 2}{2x2} squares and selects
#' the mode of each (discarding any odd boundary). In case there is more than
#' one mode, it selects the first in lexicographic order.
#' @param x an \eqn{n \times m}{n x m} binary or categorical image
#' @param random whether to have deterministic PRNG (\code{pseudo})
#'               or to use [sample()] (\code{random}). If \code{NULL},
#'               will choose the first mode. For complete reproducibility,
#'               use \code{pseudo} or \code{NULL}.
#'
#' @return an \eqn{n/2 \times m/2}{n/2 x m/2} binary or categorical image
#'
#' @keywords internal
#'
#' @noRd
downsample_2d <- function(x, random = "random") {
  dims <- dim(x)
  modepick <- 1
  newdims <- floor(dims / 2)
  newx <- matrix(nrow = newdims[1], ncol = newdims[2])
  for (i in seq(newdims[1])) {
    for (j in seq(newdims[2])) {
      xstart <- 2 * i - 1
      ystart <- 2 * j - 1
      tmpvec <- pickmode(
        c(x[xstart:(xstart + 1), ystart:(ystart + 1)]),
        random, modepick
      )
      newx[i, j] <- tmpvec[1]
      modepick <- tmpvec[2]
    }
  }
  newx
}

#' Downsampling by a factor of 2 for a 3D categorical image
#'
#' This function presumes that only the x and y axis are subsampled and the
#' z-axis is preserved. Cuts the image into a grid of \eqn{2 \times 2}{2 x 2}
#' squares and selects the mode of each (discarding any odd boundary).
#' It treats each level of z as an independent slice. In case there is more
#' than one mode, it selects the first in lexicographic order.
#' @param x an \eqn{n \times m \times q}{n x m x q} binary or categorical image
#' @param random whether to have deterministic PRNG (\code{pseudo})
#'               or to use [sample()] (\code{random}). If \code{NULL},
#'               will choose the first mode. For complete reproducibility,
#'               use \code{pseudo} or \code{NULL}.
#' @return a an \eqn{n/2 \times m/2 \times q}{n/2 x m/2 x q}
#'     binary or categorical image
#'
#' @keywords internal
#'
#' @noRd
downsample_3d_slice <- function(x, random = "random") {
  dims <- dim(x)
  newdims <- floor(c(dims[1:2] / 2, dims[3]))
  newx <- array(dim = newdims)
  for (i in seq(dims[3])) {
    newx[, , i] <- downsample_2d(x[, , i], random)
  }
  newx
}


#' Downsampling by a factor of 2 for a 3D categorical image
#'
#' This function presumes that only the x and y axis are subsampled and the
#' z-axis is preserved. Cuts the image into a grid of
#' \eqn{2 \times 2 \times 2}{2 x 2 x 2} cubes and selects
#' the mode of each (discarding any odd boundary). It treats each direction as
#' equal. In case there is more than one mode, it selects
#' the first in lexicographic order.
#' @param x an \eqn{n \times m \times q}{n x m x q} binary or categorical image
#' @param random whether to have deterministic PRNG (\code{pseudo})
#'               or to use [sample()] (\code{random}). If \code{NULL},
#'               will choose the first mode. For complete reproducibility,
#'               use \code{pseudo} or \code{NULL}.
#' @return  an \eqn{n/2 \times m/2 \times q/2}{n/2 x m/2 x q/2}
#' binary or categorical image
#'
#' @keywords internal
#'
#' @noRd
downsample_3d_cube <- function(x, random = "random") {
  dims <- dim(x)
  modepick <- 1
  newdims <- floor(dims / 2)
  newx <- array(dim = newdims)
  for (i in seq(newdims[1])) {
    for (j in seq(newdims[2])) {
      for (k in seq(newdims[3])) {
        xstart <- 2 * i - 1
        ystart <- 2 * j - 1
        zstart <- 2 * k - 1
        tmpvec <- pickmode(c(x[
          xstart:(xstart + 1),
          ystart:(ystart + 1),
          zstart:(zstart + 1)
        ]), random, modepick)
        newx[i, j, k] <- tmpvec[1]
        modepick <- tmpvec[2]
      }
    }
  }
  newx
}


#' Categorical Structural Similarity Index Measure (2D)
#'
#' The categorical structural similarity index measure for 2D categorical or
#' binary images for a single scale. This computes it using moving
#' \eqn{11 \times 11}{11 x 11} windows and is suitable for modestly-sized
#' images which are not large enough to warrant
#' looking at multiple scales.
#'
#' @inheritParams catmssim_2d
#'
#' @return the three components of the index, all less than 1.
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' set.seed(20181207)
#' x <- matrix(sample(1:4, 400, replace = TRUE), nrow = 20)
#' y <- x
#' for (i in 1:20) y[i, i] <- 1
#' for (i in 1:19) y[i, i + 1] <- 1
#' catssim_2d(x, y)
catssim_2d <- function(x, y, window = c(11, 11), method = "Cohen", ...) {
  k <- length(unique(c(x, y)))
  dims <- dim(x)
  nrow <- dims[1]
  ncol <- dims[2]
  methodflag <- method_parser(method)
  dotlist <- dots_parser(...)

  resultmatrix <- array(0, c(nrow - window[1] + 1, ncol - window[2] + 1, 3))

  for (i in seq(nrow - (window[1] - 1))) {
    for (j in seq(ncol - (window[2] - 1))) {
      subx <- x[i:(i + (window[1] - 1)), j:(j + (window[2] - 1))]
      suby <- y[i:(i + (window[1] - 1)), j:(j + (window[2] - 1))]

      resultmatrix[i, j, ] <- ssimcomponents(
        x = subx, y = suby,
        k = k, method = methodflag, c1 = dotlist[["c1"]],
        c2 = dotlist[["c2"]], sqrtgini = dotlist[["sqrtgini"]]
      )
    }
  }
  resultmatrix[resultmatrix < 0.0] <- 0.0

  colMeans(resultmatrix, na.rm = TRUE, dims = 2)
}


#' Multiscale Categorical Structural Similarity Index Measure (2D)
#'
#' The categorical structural similarity index measure for 2D categorical or
#' binary images for multiple scales. The default is to compute over 5 scales.
#'
#' @param x,y a binary or categorical image
#' @param levels how many levels of downsampling to use. By default, 5. If
#'        `weights` is specified and this is left blank, the argument
#'        will be inferred from the number of weights specified.
#' @param weights a vector of weights for the different scales. By default,
#'        equal to `rep(1,levels)/levels`. If specified, there must
#'        at least as many  weights as there are levels and the first
#'        `levels` weights will be used.
#' @param window by default 11 for 2D and 5 for 3D images,
#'     but can be specified as a
#'     vector if the window sizes differ by dimension.
#'     The vector must have the same number of
#'     dimensions as the inputted `x` and `y`.
#' @param method whether to use Cohen's kappa (`Cohen`),
#'     Jaccard Index (`Jaccard`), Dice index (`Dice`),
#'     accuracy (`accuracy`),  Rand index (`Rand`),
#'     Adjusted Rand Index (`AdjRand` or `ARI`), normalized mutual
#'     information (`NMI` or `MI`) or the adjusted mutual
#'     information, `AMI` and `ami`, as
#'     the similarity index. Note Jaccard and Dice should only be used on
#'     binary data.
#'
#' @param ... additional constants can be passed to internal functions.
#' @param random whether to have deterministic PRNG (\code{pseudo})
#'               or to use [sample()] (\code{random}). If \code{NULL},
#'               will choose the first mode. For complete reproducibility,
#'               use \code{pseudo} or \code{NULL}.
#'
#' @return a value less than 1 indicating the similarity between the images.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' x <- matrix(sample(0:3, 128^2, replace = TRUE), nrow = 128)
#' y <- x
#' for (i in 1:128) y[i, i] <- 0
#' for (i in 1:127) y[i, i + 1] <- 0
#' catmssim_2d(x, y, method = "Cohen", levels = 2) # the default
#' # now using a different similarity score (Jaccard Index)
#' catmssim_2d(x, y, method = "NMI")
catmssim_2d <- function(x, y, levels = NULL, weights = NULL, window = 11,
                        method = "Cohen", ..., random = "random") {
  if (is.null(dim(x))) stop("x is 1-dimensional")
  if (is.null(dim(y))) stop("y is 1-dimensional")
  if (length(dim(x)) != length(dim(y))) {
    stop("x and y have nonconformable dimensions.")
  }
  if (any(dim(x) != dim(y))) stop("x and y have nonconformable dimensions.")
  if (length(window) == 1) window <- c(window, window)
  weights <- level_parser(weights = weights, levels = levels)
  dotlist <- dots_parser(...)
  levels <- length(weights)
  mindim <- min(dim(x))
  minwindow <- min(window[1:2])
  if (mindim < minwindow) {
    warning("Minimum dimension should be greater than window size.
             Using only one level.")
    return(binssim(x = x, y = y, method = method, ...))
  }

  if (mindim < (2^(levels - 1)) * minwindow) {
    levels <- min(c(floor(log2(dim(x) / window[1:2]) + 1), levels))
    warning("Truncating levels because of minimum dimension.")
  }
  weights <- weights[1:levels]
  results <- matrix(0, nrow = levels, ncol = 3)
  results[1, ] <- catssim_2d(
    x = x, y = y,
    window = window, method = method, ...
  )

  if (levels > 1) {
    for (i in 2:levels) {
      x <- downsample_2d(x, random = random)
      y <- downsample_2d(y, random = random)
      results[i, ] <- catssim_2d(
        x = x, y = y,
        window = window, method = method, ...
      )
    }
  }
  results[is.na(results)] <- 1 # fix Jaccard NA results
  # use "luminosity" only from top level, use C and S from all levels
  csresults <- prod(results[, 2:3]^(weights))

  (results[levels, 1]^weights[levels]) * csresults
}


#' CatSIM (single level) for 3D slices
#'
#' Performs the 2D CatSSIM for each slice of the 3D image.
#'
#' @param x,y a binary or categorical image
#' @param window by default 11
#' @param method whether to use Cohen's kappa (`Cohen`),
#'    Jaccard Index (`Jaccard`), Dice index (`Dice`),
#'       accuracy (`accuracy`),  Rand index (`Rand`),
#'     Adjusted Rand Index (`AdjRand` or `ARI`), or normalized mutual
#'   information (`NMI` or `MI`) as
#'     the similarity index. Note Jaccard and Dice should only be used on
#' binary data.
#' @param ...
#'
#' @return SSIM componenets for the cube.
#' @keywords internal
#'
#' @noRd
#'
catssim_3d_slice <- function(x, y, window = c(11, 11), method = "Cohen", ...) {
  dims <- dim(x)
  sliceresults <- matrix(0, nrow = dims[3], ncol = 3)
  for (i in seq(dims[3])) {
    sliceresults[i, ] <- catssim_2d(
      x = x[, , i], y = y[, , i],
      window = window, method = method, ...
    )
  }
  colMeans(sliceresults, na.rm = TRUE)
}

catssim_3d_cube <- function(x, y, window = c(5, 5, 5), method = "Cohen", ...) {
  k <- length(unique(c(x, y)))
  dims <- dim(x)
  methodflag <- method_parser(method)
  dotlist <- dots_parser(...)

  cuberesults <- array(0, c(dims - window[1:3] + 1, 3))
  for (i in seq(dims[1] - (window[1] - 1))) {
    for (j in seq(dims[2] - (window[2] - 1))) {
      for (k in seq(dims[3] - (window[3] - 1))) {
        subx <- x[
          i:(i + (window[1] - 1)),
          j:(j + (window[2] - 1)),
          k:(k + (window[3] - 1))
        ]
        suby <- y[
          i:(i + (window[1] - 1)),
          j:(j + (window[2] - 1)),
          k:(k + (window[3] - 1))
        ]
        cuberesults[i, j, k, ] <- ssimcomponents(
          x = subx, y = suby,
          k = k, method = methodflag, c1 = dotlist[["c1"]],
          c2 = dotlist[["c2"]], sqrtgini = dotlist[["sqrtgini"]]
        )
      }
    }
  }
  colMeans(cuberesults, na.rm = TRUE, dims = 3)
}

#' Multiscale Categorical Structural Similarity Index Measure by Slice (3D)
#'
#' The categorical structural similarity index measure for 3D categorical or
#' binary images for multiple scales. The default is to compute over 5 scales.
#' This computes a 2D measure for each x-y slice of the z-axis
#' and then averages over the z-axis.
#'
#' @inheritParams catmssim_2d
#'
#' @return a value less than 1 indicating the similarity between the images.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' dim <- 8
#' x <- array(sample(0:4, dim^5, replace = TRUE), dim = c(dim^2, dim^2, dim))
#' y <- x
#' for (j in 1:(dim)) {
#'   for (i in 1:(dim^2)) y[i, i, j] <- 0
#'   for (i in 1:(dim^2 - 1)) y[i, i + 1, j] <- 0
#' }
#' catmssim_3d_slice(x, y, weights = c(.75, .25)) # by default method = "Cohen"
#' # compare to some simple metric:
#' mean(x == y)
catmssim_3d_slice <- function(x, y, levels = NULL, weights = NULL,
                              window = 11, method = "Cohen", ..., random = "random") {
  if (is.null(dim(x))) {
    stop("x is 1-dimensional")
  }
  if (is.null(dim(y))) {
    stop("y is 1-dimensional")
  }
  if (length(dim(x)) != length(dim(y))) {
    stop("x and y have nonconformable dimensions.")
  }
  if (any(dim(x) != dim(y))) {
    stop("x and y have nonconformable dimensions.")
  }
  weights <- level_parser(weights = weights, levels = levels)
  dotlist <- dots_parser(...)
  levels <- length(weights)
  dims <- dim(x)
  if (length(window) == 1) window <- c(window, window)
  if (length(dims) < 3) {
    stop("x and y are not 3-dimensional.")
  }
  mindim <- min(dim(x)[1:2])
  minwindow <- min(window[1:2])
  if (any(dims[1:2] < window[1:2])) {
    stop("Minimum dimension must be greater than window size.")
  }
  if (mindim < (2^(levels - 1)) * minwindow) {
    levels <- min(c(floor(log2(dims[1:2] / window[1:2]) + 1), levels))
    warning("Truncating levels because of minimum dimension.")
  }

  weights <- weights[1:levels]
  results <- matrix(0, nrow = levels, ncol = 3)

  results[1, ] <- catssim_3d_slice(
    x = x, y = y,
    window = window, method = method, ...
  )

  if (levels > 1) {
    for (j in 2:levels) {
      x <- downsample_3d_slice(x, random = random)
      y <- downsample_3d_slice(y, random = random)
      results[j, ] <- catssim_3d_slice(
        x = x, y = y,
        window = window, method = method, ...
      )
    }
  }
  results[is.na(results)] <- 1 # fixing Jaccard NAs
  # use "luminosity" only from top level, use C and S from all levels
  csresults <- prod(results[, 2:3]^(weights))

  (results[levels, 1]^weights[levels]) * csresults
}

#' Multiscale Categorical Structural Similarity Index Measure for a Cube (3D)
#'
#' The categorical structural similarity index measure for 3D
#' categorical or binary images for multiple scales.
#' The default is to compute over 5 scales.
#' This computes a 3D measure based on \eqn{5 \times 5 \times 5}{5x5x5}
#' windows by default with 5 levels of downsampling.
#'
#' @inheritParams catmssim_2d
#'
#' @return a value less than 1 indicating the similarity between the images.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' dim <- 16
#' x <- array(sample(0:4, dim^3, replace = TRUE), dim = c(dim, dim, dim))
#' y <- x
#' for (j in 1:dim) {
#'   for (i in 1:dim) y[i, i, j] <- 0
#'   for (i in 1:(dim - 1)) y[i, i + 1, j] <- 0
#' }
#' catmssim_3d_cube(x, y, weights = c(.75, .25))
#' # Now using a different similarity score
#' catmssim_3d_cube(x, y, weights = c(.75, .25), method = "Accuracy")
catmssim_3d_cube <- function(x, y, levels = NULL, weights = NULL, window = 5,
                             method = "Cohen", ..., random = "random") {
  if (is.null(dim(x))) {
    stop("x is 1-dimensional")
  }
  if (is.null(dim(y))) {
    stop("y is 1-dimensional")
  }
  if (length(dim(x)) != length(dim(y))) {
    stop("x and y have nonconformable dimensions.")
  }
  if (any(dim(x) != dim(y))) {
    stop("x and y have nonconformable dimensions.")
  }
  weights <- level_parser(weights = weights, levels = levels)
  dotlist <- dots_parser(...)
  levels <- length(weights)
  if (length(window) == 1) window <- rep(window, 3)
  if (length(window) != 3) {
    stop("Window not of length 1 or 3")
  }
  dims <- dim(x)
  if (length(dims) < 3) {
    stop("x and y are not 3-dimensional.")
  }
  mindim <- min(dim(x))
  if (mindim < 1.5 * min(window)) {
    warning("Minimum dimension must be greater than 1.5 * minimum window.")
    if (min(window) > mindim) window <- rep(mindim, 3)
    return(catssim_3d_cube(x = x, y = y, window = window, method = method, ...))
  }
  if (mindim < (2^(levels - 1)) * min(window)) {
    levels <- min(c(floor(log2(dims / window[1:3]) + 1), levels))
    warning("Truncating levels because of minimum dimension.")
  }

  weights <- weights[1:levels]
  results <- matrix(0, nrow = levels, ncol = 3)

  results[1, ] <- catssim_3d_cube(
    x = x, y = y,
    window = window, method = method, ...
  )

  if (levels > 1) {
    for (j in 2:levels) {
      x <- downsample_3d_cube(x, random = random)
      y <- downsample_3d_cube(y, random = random)
      results[j, ] <- catssim_3d_cube(
        x = x, y = y,
        window = window, method = method, ...
      )
    }
  }
  results[is.na(results)] <- 1 # fixing Jaccard NAs
  # use "luminosity" only from top level, use C and S from all levels
  csresults <- prod(results[, 2:3]^(weights))

  (results[levels, 1]^weights[levels]) * csresults
}

#' Multiscale Categorical Structural Similarity Index Measure
#'
#' The categorical structural similarity index measure for 2D or 3D categorical or
#' binary images for multiple scales. The default is to compute over 5 scales.
#' This determines whether this is a 2D or 3D image and applies the appropriate
#' windowing, weighting, and scaling. Additional arguments can be passed.
#' This is a wrapper function for the 2D and 3D functions whose functionality
#' can be accessed through the ... arguments. This function is a wrapper for the
#' [catmssim_2d()], [catmssim_3d_slice()], and
#' [catmssim_3d_cube()]
#' functions.
#'
#' @param x,y a binary or categorical image
#' @param ... additional arguments, such as window, can be passed
#'        as well as arguments for internal functions.
#' @param cube for the 3D method, whether to use the true 3D method
#'        (cube or [catmssim_3d_cube()])
#'        or compute the metric using 2D slices which are then averaged
#'        ([catmssim_3d_slice()]). By default, `TRUE`,
#'        which evaluates as a cube. `FALSE` will treat it as
#'        2D slices.
#' @param levels how many levels of downsampling to use. By default, 5. If
#'        `weights` is specified and this is left blank, the argument
#'        will be inferred from the number of weights specified.
#' @param weights a vector of weights for the different scales. By default,
#'        equal to `rep(1,levels)/levels`. If specified, there must
#'        at least as many  weights as there are levels and the first
#'        `levels` weights will be used.
#' @param method whether to use Cohen's kappa (`Cohen`),
#'     Jaccard Index (`Jaccard`), Dice index (`Dice`),
#'     accuracy (`accuracy`),  Rand index (`Rand`),
#'     Adjusted Rand Index (`AdjRand` or `ARI`), normalized mutual
#'     information (`NMI` or `MI`), or adjusted mutual information
#'     (`AMI`) as the similarity index.
#'     Note Jaccard and Dice should only be used on binary data.
#' @param window by default 11 for 2D and 5 for 3D images, but can be
#'     specified as a vector if the window sizes differ by dimension.
#'     The vector must have the same number of
#'     dimensions as the inputted `x` and `y`.
#' @return a value less than 1 indicating the similarity between the images.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' dim <- 16
#' x <- array(sample(0:4, dim^3, replace = TRUE), dim = c(dim, dim, dim))
#' y <- x
#' for (j in 1:dim) {
#'   for (i in 1:dim) y[i, i, j] <- 0
#'   for (i in 1:(dim - 1)) y[i, i + 1, j] <- 0
#' }
#' catsim(x, y, weights = c(.75, .25))
#' # Now using a different similarity score
#' catsim(x, y, levels = 2, method = "accuracy")
#' # with the slice method:
#' catsim(x, y, weights = c(.75, .25), cube = FALSE, window = 8)
catsim <- function(x, y, ..., cube = TRUE, levels = NULL, weights = NULL,
                   method = "Cohen", window = NULL) {
  ##  old weights: c(0.0448, 0.2856, 0.3001, 0.2363, 0.1333)
  if (is.null(dim(x))) {
    stop("x is 1-dimensional")
  }
  if (is.null(dim(y))) {
    stop("y is 1-dimensional")
  }
  if (length(dim(x)) != length(dim(y))) {
    stop("x and y have nonconformable dimensions.")
  }
  if (any(dim(x) != dim(y))) {
    stop("x and y have nonconformable dimensions.")
  }
  dims <- dim(x)
  if (is.null(window)) {
    if (cube && length(dims) == 3) {
      window <- 5
    } else {
      window <- 11
    }
  }
  weights <- level_parser(weights = weights, levels = levels)
  if (length(dims) == 2) {
    catmssim_2d(x, y,
      weights = weights, method = method,
      levels = levels, window = window, ...
    )
  } else if (cube) {
    catmssim_3d_cube(x, y,
      weights = weights, levels = levels,
      method = method, window = window, ...
    )
  } else {
    catmssim_3d_slice(x, y,
      weights = weights, levels = levels,
      method = method, window = window, ...
    )
  }
}
