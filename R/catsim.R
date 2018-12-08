
#' Means Function
#'
#' @param x binary or categorical image or vector
#' @param y binary or categorical image or vector
#' @param c1 small constant
#'
#' @return measure of similarity of means
#' @keywords internal
#'
#' @noRd
meansfunc <- function(x,y, c1 = 0.01){
  levels <- levels(factor(c(x,y)))
  x = factor(x, levels = levels)
  y = factor(y, levels = levels)
  tablexy <- table(x,y)
  tablex <- rowSums(tablexy)
  tabley <- colSums(tablexy)

  (2*sum(tablex*tabley) + c1)/(sum(tablex^2) + sum(tabley^2) + c1)
}

#' Gini index
#'
#' @param x binary or categorical image or vector
#'
#' @return Gini index
#' @keywords internal
#'
#' @noRd
#' @examples
#' \dontrun{
#' x <- rep(c(1:4),5)
#' gini(x)
#' }
#'
gini <- function(x){
 1 - (sum(table(x)^2)/length(x)^2)
}


#' Corrected Gini index
#'
#' @param x binary or categorical image or vector
#' @return Gini index (corrected so that max possible is 1)
#' @keywords internal
#' @noRd
#' @examples
#' \dontrun{
#' x <- rep(c(1:4),5)
#' ginicorr(x, 4)
#' }
ginicorr <- function(x, k){
  # k <- length(table(x))
  if (k > 1){
  gini(x)/(1 - 1/k)
  } else gini(x)
}

#' Modified Gini index
#'
#' @param x binary or categorical image or vector
#' @return Modified Gini index (square root of squared frequencies)
#' @keywords internal
#' @noRd
#' @examples
#' \dontrun{
#' x <- rep(c(1:4),5)
#' sqrtgini(x)
#' }
sqrtgini <-  function(x){
  1 - sqrt(sum(table(x)^2)/length(x)^2)
}

#' Modified Corrected Gini index
#'
#' @param x binary or categorical image or vector
#' @return Modified corrected Gini index (square root of squared frequencies) - max possible is 1.
#' @keywords internal
#' @noRd
#' @examples
#' \dontrun{
#' x <- rep(c(1:4),5)
#' sqrtginicorr(x, 4)
#' }
sqrtginicorr <- function(x, k){
  # k <- length(table(x))
  if (k > 1) {
  sqrtgini(x)/(1 - 1/sqrt(k))
  } else sqrtgini(x)
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

cfunc <- function(x, y, c2 = 0.01, k){
  varx <- ginicorr(x, k)
  vary <- ginicorr(y, k)

  (2*sqrt(varx * vary)+c2)/(varx + vary + c2)

}

#' Covariance function (internal)
#'
#' @param x binary or categorical image or vector
#' @param y binary or categorical image or vector
#' @noRd
#'
#' @return Covariance function (Adjusted Rand Index)
#' @keywords internal
#' \dontrun{
#' x <- rep(1:4,4)
#' y <- c(rep(1:4,3),rep(4,4))
#' sfunc(x,y)
#' }
sfunc <- function(x, y){
  # should take care of this beforehand!
  numx <- as.numeric(x)
  numy <- as.numeric(y)
  EMCluster::RRand(numx, numy)$adjRand
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
#' @param ... Constants can be passed to the components of the index.
#'
#' @return Structural similarity index.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' x <- matrix(sample(1:4, 400, replace = TRUE), nrow=20)
#' y <- x
#' for (i in 1:20) y[i, i] = 1
#' for (i in 1:19) y[i, i+1] = 1
#' binssim(x,y)
binssim <- function(x, y, alpha = 1, beta = 1, gamma = 1, ...){
  if (length(x) != length(y)) stop("x and y must be the same size.")
  k = length(unique(c(x,y)))
  (meansfunc(x, y,...)^alpha)*(cfunc(x, y, k = k, ...)^beta)*(sfunc(x, y)^gamma)
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
ssimcomponents <- function(x, y, ...){
  k = length(unique(c(x,y)))
  c((meansfunc(x, y,...)),(cfunc(x, y, k = k, ...)),(sfunc(x, y)))
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
  if (any(newdims < 1)){
    warning("Cannot downsample, at least one dimension is of length 1")
    return(x)
  }
  newx <- matrix(nrow = newdims[1], ncol = newdims[2])
  for (i in 1:newdims[1]){
    for (j in 1:newdims[2]){
      xstart <- 2*i -1
      ystart <- 2*j - 1
      newx[i,j] = pickmode(c(x[xstart:(xstart+1), ystart:(ystart+1)]))
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
  if (any(newdims < 1)){
    warning("Cannot downsample, at least one dimension is of length 1")
    return(x)
  }
  newx <- array(dim = newdims)
  for (i in 1:dims[3]){
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
  if (any(newdims < 1)){
    warning("Cannot downsample, at least one dimension is of length 1")
    return(x)
  }
  newx <- array(dim = newdims)
  for (i in 1:newdims[1]){
    for (j in 1:newdims[2]){
      for (k in 1:newdims[3]){
      xstart <- 2*i -1
      ystart <- 2*j - 1
      zstart <- 2*k - 1
      newx[i,j,k] = pickmode(c(x[xstart:(xstart+1), ystart:(ystart+1), zstart:(zstart+1)]))
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
#' @param ... additional constants can be passed to internal functions.
#'
#' @return a value less than 1 indicating the similarity between the images.
#' @keywords internal
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
catssim_2d <- function(x,y,...){
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions')
  levels <- levels(factor(c(x,y)))
  dims = dim(x)
  nrow = dims[1]
  ncol = dims[2]
  if (any(dims < 8)) return(ssimcomponents(factor(x, levels),factor(y,levels)))

  resultmatrix = array(0, dim = c((nrow-7), (ncol-7),3))
  for (i in 1:(nrow-7)){
    for (j in 1:(ncol-7)){
      place = j+(ncol-7)*(i-1)
      subx = x[i:(i+7),j:(j+7)]
      suby = y[i:(i+7),j:(j+7)]

      subx = factor(subx, levels = levels)
      suby = factor(suby, levels = levels)
      resultmatrix[i, j, ] = ssimcomponents(subx, suby, ...)
    }
  }
  colMeans(resultmatrix, dims = 2)
  #resultmatrix
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
catmssim_2d <- function(x, y, weights = c(0.0448, 0.2856, 0.3001, 0.2363, 0.1333), ...){
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions')
  levels = length(weights)
  mindim <- min(dim(x))
  if(mindim < 8) stop("Minimum dimension must be greater than 8.")
  if (mindim < 128) levels = min(c(floor(log2(mindim) - 2),levels))
  weights = weights[1:levels]
  results = matrix(0, nrow = levels, ncol = 3)
  results[1,] = catssim_2d(x,y)

  if ( levels > 1){
    for (i in 2:levels){
        x = downsample_2d(x)
        y = downsample_2d(y)
        results[i,] = catssim_2d(x,y)
      }
    }

  # use "luminosity" only from top level, use C and S from all levels
  csresults <- prod(results[,2:3]^(weights))
  #print(results)

  (results[levels,1]^weights[levels])*csresults

}


catssim_3d_slice <- function(x,y,...){
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions')
  dims = dim(x)
  sliceresults = matrix(0, nrow = dims[3], ncol = 3)
  for (i in 1:dims[3]) sliceresults[i,] = catssim_2d(x[,,i],y[,,i],...)
  colMeans(sliceresults)
}

catssim_3d_cube <- function(x,y,...){
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions')
  levels <- levels(factor(c(x,y)))
  dims = dim(x)
  if (any(dims < 4)) return(ssimcomponents(factor(x, levels), factor(y,levels)))

  # this gets messy
  cuberesults = array(0, dim = c((dims - 3), 3))
  for (i in 1:(dims[1] - 3)) {
    for (j in 1:(dims[2] - 3)) {
      for (k in 1:(dims[3] - 3)) {
        subx = x[i:(i + 3), j:(j + 3), k:(k + 3)]
        suby = y[i:(i + 3), j:(j + 3), k:(k + 3)]

        subx = factor(subx, levels = levels)
        suby = factor(suby, levels = levels)
        cuberesults[i, j, k, ] = ssimcomponents(subx, suby, ...)
      }
    }
  }
  #for (i in 1:dims[3]) cuberesults[i,] = catssim_2d(x[,,i],y[,,i],...)
  colMeans(cuberesults, dims = 3)
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
#' @param ... additional constants can be passed to internal functions.
#'
#' @return a value less than 1 indicating the similarity between the images.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' x <- array(sample(1:4, 16^3, replace = TRUE), dim = c(32,32,32))
#' y <- x
#' for (j in 1:32){
#' for (i in 1:32) y[i, i, j] = 1
#' for (i in 1:31) y[i, i+1, j] = 1
#' }
#' catmssim_3d_slice(x,y, weights = c(.75,.25))
catmssim_3d_slice <- function(x, y, weights = c(0.0448, 0.2856, 0.3001, 0.2363, 0.1333), ...){
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions')
  levels = length(weights)
  dims = dim(x)
  if (length(dims) < 3) stop('x and y are not 3-dimensional')
  mindim <- min(dim(x)[1:2])
  if (mindim < 8) stop("Minimum dimension must be greater than 8.")
  if (mindim < 128) levels = min(c(floor(log2(mindim) - 2),levels))

  weights = weights[1:levels]
  results = matrix(0, nrow = levels, ncol = 3)

  results[1,] = catssim_3d_slice(x,y,...)

  if ( levels > 1){
    for (j in 2:levels){
      x = downsample_3d_slice(x)
      y = downsample_3d_slice(y)
      results[j,] = catssim_3d_slice(x,y,...)
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
#' This computes a 3D measure
#'
#' @param x a binary or categorical image
#' @param y a binary or categorical image
#' @param weights a vector of weights for the different scales. By default,
#'     five different scales are used.
#' @param ... additional constants can be passed to internal functions.
#'
#' @return a value less than 1 indicating the similarity between the images.
#' @export
#'
#' @examples
#' x <- array(sample(1:4, 16^3, replace = TRUE), dim = c(32,32,32))
#' y <- x
#' for (j in 1:32){
#' for (i in 1:32) y[i, i, j] = 1
#' for (i in 1:31) y[i, i+1, j] = 1
#' }
#' catmssim_3d_cube(x,y, weights = c(.75,.25))
catmssim_3d_cube <- function(x, y, weights = c(0.0448, 0.2856, 0.3001, 0.2363, 0.1333), ...){
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions')
  levels = length(weights)
  dims = dim(x)
  if (length(dims) < 3) stop('x and y are not 3-dimensional')
  mindim <- min(dim(x))
  if (mindim < 8) stop("Minimum dimension must be greater than 8.")
  if (mindim < 128) levels = min(c(floor(log2(mindim) - 2),levels))

  weights = weights[1:levels]
  results = matrix(0, nrow = levels, ncol = 3)

  results[1,] = catssim_3d_cube(x,y,...)

  if ( levels > 1){
    for (j in 2:levels){
      x = downsample_3d_cube(x)
      y = downsample_3d_cube(y)
      results[j,] = catssim_3d_cube(x,y,...)
    }
  }

  # use "luminosity" only from top level, use C and S from all levels
  csresults <- prod(results[,2:3]^(weights))
  #print(results)

  (results[levels,1]^weights[levels])*csresults

}

