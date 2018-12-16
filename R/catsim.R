#' @useDynLib catsim, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


#' Means Function (internal)
#'
#' @param x binary or categorical image or vector
#' @param y binary or categorical image or vector
#' @param c1 small constant
#'
#' @return measure of similarity of means
#' @keywords internal
#'
#' @noRd
meansfunc <- function(x,y, c1 = 0.01) {
  C_meansfunc(x, y, c1)
}

#' Gini-Simpson index
#'
#' A measure of diversity that goes by a number of different names, such as
#' the probability of interspecific encounter or the Gibbs-Martin index.
#' It is 1 - sum(p_i^2), where p_i is the probability of observing class i.
#'
#' @param x binary or categorical image or vector
#'
#' @return Gini-Simpson index (between 0 and 1)
#' @export
#'
#'
#' @examples
#'
#' x <- rep(c(1:4),5)
#' gini(x)
#'
#'
gini <- function(x){
  x <- as.numeric(x)
  C_gini(x)
}


#' Corrected Gini-Simpson index
#'
#' @param x binary or categorical image or vector
#' @param k number of categories
#' @return Gini-Simpson index (corrected based on the number of categories
#' so that max possible is 1)
#' @export
#'
#' @examples
#'
#' x <- rep(c(1:4),5)
#' ginicorr(x, 4)
#'
ginicorr <- function(x, k) {
  if (k > 1) {
    C_gini(x)/(1 - 1/k)
  } else 1
}

#' Modified Gini-Simpson index
#'
#' The Gini-Simpson index, except with the square root of the summed squared
#' probabilities.
#'
#' @param x binary or categorical image or vector
#' @return Modified Gini-Simpson index (square root of squared frequencies)
#' @export
#'
#' @examples
#'
#' x <- rep(c(1:4),5)
#' sqrtgini(x)
#'
sqrtgini <-  function(x){
  1 - sqrt(1 - C_gini(x))
}

#' Modified Corrected Gini index
#'
#' The Gini-Simpson index, except with the square root of the summed squared
#' probabilities.
#'
#' @param x binary or categorical image or vector
#' @param k number of categories
#' @return Modified corrected Gini index (square root of squared frequencies) -
#' max possible is 1.
#' @export
#'
#' @examples
#'
#' x <- rep(c(1:4),5)
#' sqrtginicorr(x, 4)
#'
sqrtginicorr <- function(x, k){
  if (k > 1) {
    sqrtgini(x)/(1 - 1/sqrt(k))
  } else 1
}

#' Variance function (internal)
#'
#' @param x binary or categorical image or vector
#' @param y binary or categorical image or vector
#' @param c2 small constant
#' @noRd
#' @return variance function
#' @keywords internal
#' @example
#' \dontrun{
#' x <- rep(1:4,4)
#' y <- c(rep(1:4,3),rep(4,4))
#' cfunc(x,y, k = 4)
#' }

cfunc <- function(x, y, c2 = 0.01, k, sqrtgini = FALSE){
  C_cfunc(x, y, c2, k, sqrtgini)
}

#' Covariance function (internal)
#'
#' @param x binary or categorical image or vector
#' @param y binary or categorical image or vector
#' @noRd
#'
#' @return Covariance function (Cohen's Kappa)
#' @keywords internal
#' \dontrun{
#' x <- rep(1:4,4)
#' y <- c(rep(1:4,3),rep(4,4))
#' sfunc(x,y)
#' }
sfunc <- function(x, y){
    x <- as.vector(x)
    y <- as.vector(y)
    C_Cohen(x,y)
}


#' Categorical Structural Similarity Index Metric (whole image)
#'
#' This computes the categorical or binary structural similarity index metric
#' on a whole-image scale.
#'
#' @param x binary or categorical image
#' @param y binary or categorical image
#' @param alpha normalizing parameter, by default 1
#' @param beta normalizing parameter, by default 1
#' @param gamma normalizing parameter, by default 1
#' @param c1 small normalization constant for the c function, by default 0.01
#' @param c2 small normalization constant for the s function, by default 0.01
#' @param ... Constants can be passed to the components of the index.
#'
#' @return Structural similarity index.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' x <- matrix(sample(1:4, 10000, replace = TRUE), nrow=100)
#' y <- x
#' for (i in 1:100) y[i, i] = 1
#' for (i in 1:99) y[i, i+1] = 1
#' binssim(x,y)
binssim <- function(x, y, alpha = 1, beta = 1, gamma = 1, c1 = 0.01, c2 = 0.01, ...){
  if (length(x) != length(y)) stop("x and y must be the same size.")
  k = length(unique(c(x,y)))
  (meansfunc(x, y, c1)^alpha)*(cfunc(x, y, c2, k = k, ...)^beta)*(sfunc(x, y)^gamma)
}

#' Categorical SSIM Components
#'
#' @param x binary or categorical image
#' @param y binary or categorical image

#' @param ... constants can be passed to the internal functions
#' @noRd
#' @return the three components of the Categorical SSIM.
#' @keywords internal
#'
ssimcomponents <- function(x, y, k, c1 = 0.01, c2 = 0.01, ...){
  #k = length(levels)
  #levels <- levels(factor(c(x,y)))
  c((meansfunc(x, y, c1)),(cfunc(x, y, c2, k = k, ...)),(sfunc(x, y)))
}


#' Picks the first mode
#'
#' @param x binary or categorical image
#'
#' @return The first mode in the default order returned by R
#' @noRd
#' @keywords internal
#'
pickmode <- function(x){
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  #sample(ux[tab == max(tab)],1)
  # the random selection made things look too bad! so just taking the first
  # mode programmatically.
  ux[which.max(tab)]
}


#' Downsampling by a factor of 2 for a 2D categorical image
#'
#' Cuts the image into a grid of 2x2 squares and selects the mode of each
#' (discarding any odd boundary). In case there is more than one mode, it selects
#' the first in lexicographic order.
#' @param x an n x m binary or categorical image
#'
#' @return a an n/2 x m/2 binary or categorical image
#'
#' @keywords internal
#'
#' @noRd
downsample_2d <- function(x){
  dims = dim(x)
  newdims = floor(dims/2)
  if (any(newdims < 1)) {
    warning("Cannot downsample, at least one dimension is of length 1")
    return(x)
  }
  newx <- matrix(nrow = newdims[1], ncol = newdims[2])
  for (i in 1:newdims[1]) {
    for (j in 1:newdims[2]) {
      xstart <- 2*i - 1
      ystart <- 2*j - 1
      newx[i,j] = pickmode(c(x[xstart:(xstart + 1), ystart:(ystart + 1)]))
    }
  }
  newx
}

#' Downsampling by a factor of 2 for a 3D categorical image
#'
#' This function presumes that only the x and y axis are subsampled and the
#' z-axis is preserved. Cuts the image into a grid of 2x2 squares and selects
#' the mode of each (discarding any odd boundary). It treats each level of z as
#' an independent slice. In case there is more than one mode, it selects
#' the first in lexicographic order.
#' @param x an n x m x q binary or categorical image
#'
#' @return a an n/2 x m/2 x q binary or categorical image
#'
#' @keywords internal
#'
#' @noRd
downsample_3d_slice <- function(x){
  dims = dim(x)
  newdims = floor(c(dims[1:2]/2,dims[3]))
  if (any(newdims < 1)) {
    warning("Cannot downsample, at least one dimension is of length 1.")
    return(x)
  }
  newx <- array(dim = newdims)
  for (i in 1:dims[3]) {
    newx[,,i] = downsample_2d(x[,,i])
  }
  newx
}


#' Downsampling by a factor of 2 for a 3D categorical image
#'
#' This function presumes that only the x and y axis are subsampled and the
#' z-axis is preserved. Cuts the image into a grid of 2x2x2 cubes and selects
#' the mode of each (discarding any odd boundary). It treats each direction as
#' equal. In case there is more than one mode, it selects
#' the first in lexicographic order.
#' @param x an n x m x q binary or categorical image
#'
#' @return a an n/2 x m/2 x q/2 binary or categorical image
#'
#' @keywords internal
#'
#' @noRd
downsample_3d_cube <- function(x){
  dims = dim(x)
  newdims = floor(dims/2)
  if (any(newdims < 1)) {
    warning("Cannot downsample, at least one dimension is of length 1.")
    return(x)
  }
  newx <- array(dim = newdims)
  for (i in 1:newdims[1]) {
    for (j in 1:newdims[2]) {
      for (k in 1:newdims[3]) {
        xstart <- 2*i - 1
        ystart <- 2*j - 1
        zstart <- 2*k - 1
        newx[i, j, k] = pickmode(c(x[xstart:(xstart + 1), ystart:(ystart + 1), zstart:(zstart + 1)]))
      }
    }
  }
  newx
}


#' Categorical Structural Similarity Index Measure (2D)
#'
#' The categorical structural similary index measure for 2D categorical or binary
#' images for a single scale. This computes it using moving 8x8 windows and is
#' suitable for modestly-sized images which are not large enough to warrant
#' looking at multiple scales.
#'
#' @param x a binary or categorical image
#' @param y a binary or categorical image
#' @param window window size, by default 8
#' @param ... additional constants can be passed to internal functions.
#'
#' @return the three components of the index, all less than 1.
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' set.seed(20181207)
#' x <- matrix(sample(1:4, 400, replace = TRUE), nrow=20)
#' y <- x
#' for (i in 1:20) y[i, i] = 1
#' for (i in 1:19) y[i, i+1] = 1
#' catssim_2d(x,y)
#' }
catssim_2d <- function(x,y, window = 8,...){
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions.')
  #levels <- levels(factor(c(x,y)))
  k = length(unique(c(x,y)))
  dims = dim(x)
  nrow = dims[1]
  ncol = dims[2]
  if (any(dims < (window - 1))) return(ssimcomponents((x), (y), k, ...))
  resultmatrix = c(0,0,0)

  for (i in 1:(nrow - (window - 1))) {
    for (j in 1:(ncol - (window - 1))) {
      #place = j + (ncol - (window - 1))*(i - 1)
      subx = x[i:(i + (window - 1)), j:(j + (window - 1))]
      suby = y[i:(i + (window - 1)), j:(j + (window - 1))]

      resultmatrix = resultmatrix + ssimcomponents(subx, suby, k, ...)
    }
  }
  resultmatrix[resultmatrix < 0.0] <- 0.0

  resultmatrix / ((nrow - (window - 1)) * (ncol - (window - 1)))

}


#' Multiscale Categorical Structural Similarity Index Measure (2D)
#'
#' The categorical structural similary index measure for 2D categorical or binary
#' images for multiple scales. The default is to compute over 5 scales.
#'
#' @param x a binary or categorical image
#' @param y a binary or categorical image
#' @param weights a vector of weights for the different scales. By default,
#'     five different scales are used.
#' @param window window size, by default 8.
#' @param ... additional constants can be passed to internal functions.
#'
#' @return a value less than 1 indicating the similarity between the images.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' x <- matrix(sample(1:4, 128^2, replace = TRUE), nrow=128)
#' y <- x
#' for (i in 1:128) y[i, i] = 1
#' for (i in 1:127) y[i, i+1] = 1
#' catmssim_2d(x,y)
#'
catmssim_2d <- function(x, y, weights = c(0.0448, 0.2856, 0.3001, 0.2363, 0.1333), window = 8, ...){
  # the weights are from the original MS-SSIM program
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions.')
  levels = length(weights)
  mindim <- min(dim(x))
  if (mindim < window) {
    warning("Minimum dimension should be greater than window size. Using only one level.")
    return(binssim(x,y,...))
  }

  if (mindim < (2^(levels - 1))*window) {
    levels = min(c(floor(log2(mindim/window) + 1),levels))
    warning("Truncating levels because of minimum dimension.")
  }
  weights = weights[1:levels]
  results = matrix(0, nrow = levels, ncol = 3)
  results[1,] = catssim_2d(x, y, window,...)

  if ( levels > 1) {
    for (i in 2:levels) {
      x = downsample_2d(x)
      y = downsample_2d(y)
      results[i,] = catssim_2d(x, y, window, ...)
    }
  }

  # use "luminosity" only from top level, use C and S from all levels
  csresults <- prod(results[,2:3]^(weights))
  #print(results)

  (results[levels,1]^weights[levels])*csresults

}


#' CatSSIM for 3D slices
#'
#' Performs the 2D CatSSIM for each slice of the 3D image.
#'
#' @param x a binary or categorical image
#' @param y a binary or categorical image
#' @param window by default 8
#' @param ...
#'
#' @return SSIM componenets for the cube.
#' @keywords internal
#'
#' @noRd
#'
catssim_3d_slice <- function(x, y, window = 8,...){
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions.')
  dims = dim(x)
  sliceresults = matrix(0, nrow = dims[3], ncol = 3)
  for (i in 1:dims[3]) sliceresults[i,] = catssim_2d(x[,,i],y[,,i], window, ...)
  colMeans(sliceresults)
}

catssim_3d_cube <- function(x, y, window = 4, ...){
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions.')
  #levels <- levels(factor(c(x,y)))
  k = length(unique(c(x,y)))
  dims = dim(x)
  if (any(dims < window)) return(ssimcomponents((x), (y)), k, ...)

  cuberesults = c(0,0,0)
  for (i in 1:(dims[1] - (window-1))) {
    for (j in 1:(dims[2] - (window-1))) {
      for (k in 1:(dims[3] - (window-1))) {
        subx = x[i:(i + (window-1)), j:(j + (window-1)), k:(k + (window-1))]
        suby = y[i:(i + (window-1)), j:(j + (window-1)), k:(k + (window-1))]

        cuberesults = cuberesults + ssimcomponents(subx, suby, k, ...)
      }
    }
  }

  (cuberesults / prod(dims - (window-1)) )
}

#' Multiscale Categorical Structural Similarity Index Measure by Slice (3D)
#'
#' The categorical structural similary index measure for 3D categorical or binary
#' images for multiple scales. The default is to compute over 5 scales.
#' This computes a 2D measure for each x-y slice of the z-axis
#' and then averages over the z-axis.
#'
#' @param x a binary or categorical image
#' @param y a binary or categorical image
#' @param weights a vector of weights for the different scales. By default,
#'     five different scales are used.
#' @param window window size, by default 8.
#' @param ... additional constants can be passed to internal functions.
#'
#' @return a value less than 1 indicating the similarity between the images.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' dim = 16
#' x <- array(sample(1:4, dim^3, replace = TRUE), dim = c(dim,dim,dim))
#' y <- x
#' for (j in 1:dim){
#' for (i in 1:dim) y[i, i, j] = 1
#' for (i in 1:(dim-1)) y[i, i+1, j] = 1
#' }
#' catmssim_3d_slice(x,y, weights = c(.75,.25))
catmssim_3d_slice <- function(x, y, weights = c(0.0448, 0.2856, 0.3001, 0.2363, 0.1333), window = 8, ...){
  # the weights are from the original MS-SSIM program
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions.')
  levels = length(weights)
  dims = dim(x)
  if (length(dims) < 3) stop('x and y are not 3-dimensional.')
  mindim <- min(dim(x)[1:2])
  if (mindim < window) stop("Minimum dimension must be greater than window size.")
  if (mindim < (2^(levels - 1))*window) {
    levels = min(c(floor(log2(mindim/window) + 1),levels))
    warning("Truncating levels because of minimum dimension.")
  }

  weights = weights[1:levels]
  results = matrix(0, nrow = levels, ncol = 3)

  results[1,] = catssim_3d_slice(x, y, window,...)

  if ( levels > 1) {
    for (j in 2:levels) {
      x = downsample_3d_slice(x)
      y = downsample_3d_slice(y)
      results[j,] = catssim_3d_slice(x, y, window,...)
    }
  }

  # use "luminosity" only from top level, use C and S from all levels
  csresults <- prod(results[,2:3]^(weights))
  #print(results)

  (results[levels,1]^weights[levels])*csresults

}

#' Multiscale Categorical Structural Similarity Index Measure for a Cube (3D)
#'
#' The categorical structural similary index measure for 3D categorical or binary
#' images for multiple scales. The default is to compute over 5 scales.
#' This computes a 3D measure based on 4x4x4 windows by default with 5 levels of downsampling.
#'
#' @param x a binary or categorical image
#' @param y a binary or categorical image
#' @param weights a vector of weights for the different scales. By default,
#'     five different scales are used.
#' @param window size of window, by default 4
#' @param ... additional constants can be passed to internal functions.
#'
#' @return a value less than 1 indicating the similarity between the images.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' dim = 16
#' x <- array(sample(1:4, dim^3, replace = TRUE), dim = c(dim,dim,dim))
#' y <- x
#' for (j in 1:dim){
#' for (i in 1:dim) y[i, i, j] = 1
#' for (i in 1:(dim-1)) y[i, i+1, j] = 1
#' }
#' catmssim_3d_cube(x,y, weights = c(.75,.25))
catmssim_3d_cube <- function(x, y, weights = c(0.0448, 0.2856, 0.3001, 0.2363, 0.1333), window = 4, ...){
  # the weights are from the original MS-SSIM program
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions.')
  levels = length(weights)
  dims = dim(x)
  if (length(dims) < 3) stop('x and y are not 3-dimensional.')
  mindim <- min(dim(x))
  if (mindim < 2*window) {
    warning("Minimum dimension must be greater than 2 * window.")
    catssim_3d_cube(x,y,...)
  }
  if (mindim < (2^(levels - 1))*window) {
    levels = min(c(floor(log2(mindim/window) + 1),levels))
    warning("Truncating levels because of minimum dimension.")
  }

  weights = weights[1:levels]
  results = matrix(0, nrow = levels, ncol = 3)

  results[1,] = catssim_3d_cube(x,y,window,...)

  if ( levels > 1) {
    for (j in 2:levels) {
      x = downsample_3d_cube(x)
      y = downsample_3d_cube(y)
      results[j,] = catssim_3d_cube(x,y,window,...)
    }
  }

  # use "luminosity" only from top level, use C and S from all levels
  csresults <- prod(results[,2:3]^(weights))
  #print(results)

  (results[levels,1]^weights[levels])*csresults

}

#' Adjusted Rand Index
#'
#' Computes the adjusted Rand index for two
#' inputs. These inputs should be binary or categorical and of the same length.
#' It also computes the PSNR, which is generalized here as simply
#' -10 log10(MSE). The adjusted Rand index and Cohen's Kappa are used as a measure of
#' the similarity of the structure of the two images. A small constant is added to the numerator
#' and denominator of the Rand index to ensure stability, as it is possible to have a zero
#' denominator. The PSNR can be infinite if the error rate is 0.
#'
#'
#' @param x a numeric or factor vector or image
#' @param y a numeric or factor vector or image
#'
#' @return The Rand index, the Adjusted Rand Index, the PSNR, and Cohen's Kappa.
#'
#' @references Lawrence Hubert and Phipps Arabie (1985).
#' "Comparing partitions". Journal of Classification. 2 (1): 193–218. \doi{10.1007/BF01908075}
#'
#'  W. M. Rand (1971). "Objective criteria for the evaluation of clustering methods".
#'  Journal of the American Statistical Association. American Statistical Association. 66 (336): 846–850.
#'  \doi{10.2307/2284239}
#'
#'  Cohen, Jacob (1960). "A coefficient of agreement for nominal scales". Educational and Psychological Measurement. 20 (1): 37–46. \doi{10.1177/001316446002000104}
#'
#' @export
#'
#' @examples
#'
#' x <- rep(1:5, 5)
#' y <- c(rep(1:5, 4),rep(1,5))
#' AdjRandIndex(x, y)
AdjRandIndex <- function(x,y){
  if (length(x) != length(y)) stop("x and y have differing lengths.")
  n <- length(x)
  a <- sum(x == y)
  Rand <- a/n
  # this disagrees with some implementations I've seen but I don't know how this is wrong.
  Cohen <- C_Cohen(x,y)
  x <- as.numeric(x)
  y <- as.numeric(y)
  AdjRand <- C_AdjRand(x,y)
  list(Rand = Rand,
       AdjRand = AdjRand,
       PSNR = -10 * log10(Rand),
       Cohen = Cohen)
}
