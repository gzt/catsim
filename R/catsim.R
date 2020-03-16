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

#' @title Diversity Indexes
#' @name gini
#'
#' @description \code{gini()} is a measure of diversity that goes by a number of different names, such as
#' the probability of interspecific encounter or the Gibbs-Martin index.
#' It is \eqn{1 - sum(p_i^2)}, where \eqn{p_i} is the probability of observing class i.
#'
#' @param x binary or categorical image or vector
#'
#' @return The index (between 0 and 1), with 0 indicating no variation and 1 being maximal.
#'    The Gini index is bounded above by \eqn{1-1/k} for a group with \code{k} categories.
#'    The modified index is bounded above by \eqn{1-1/\sqrt{k}}.
#'    The corrected indexes fix this by dividing by the maximum.
#' @export

#' @examples
#' x <- rep(c(1:4),5)
#' gini(x)
#'
#'
gini <- function(x){
  x <- as.numeric(x)
  C_gini(x)
}


#' @name Corrected Gini-Simpson index
#'
#' @description The corrected Gini-Simpson index, \code{ginicorr} takes the index and  corrects it
#' so that the maximum possible is 1. If there are \code{k} categories,
#' the maximum possible of the uncorrected index is \eqn{1-1/k}.
#' It corrects the index by dividing by the maximum. \code{k} must be specified.
#'
#' @param k number of categories
#' @export
#' @rdname gini
#' @examples
#'
#' x <- rep(c(1:4),5)
#' ginicorr(x, 4)
#'
ginicorr <- function(x, k) {
    if(k < length(unique(x))) stop("more unique values in x than k")
  if (k > 1) {
    C_gini(x)/(1 - 1/k)
  } else 1
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
#' x <- rep(c(1:4),5)
#' sqrtgini(x)
#'
sqrtgini <-  function(x){
  1 - sqrt(1 - C_gini(x))
}

#' @name Modified Corrected Gini index
#'
#' @description The modified corrected Gini index then
#' corrects the modified index for the number of categories, \code{k}.
#'
#' @export
#' @rdname gini
#' @examples
#'
#' x <- rep(c(1:4),5)
#' sqrtginicorr(x, 4)
#'
sqrtginicorr <- function(x, k){
    if(k < length(unique(x))) stop("more unique values in x than k")

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
#' @examples
#' \dontrun{
#' x <- rep(1:4,4)
#' y <- c(rep(1:4,3),rep(4,4))
#' cfunc(x,y, k = 4)
#'
#' }
cfunc <- function(x, y, c2 = 0.01, k, sqrtgini = TRUE){
  C_cfunc(x, y, c2, k, sqrtgini)
}



#' Jaccard Index
#'
#' @param x binary image or vector
#' @param y binary image or vector
#'
#' @return Jaccard index - note this privileges 1 over 0
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' x <- rep(c(0,1), 7)
#' y <- rep(c(1,1), 7)
#' jaccard(x,y)
#'}
jaccard <- function(x, y){
    if (length(x) != length(y)) stop("x and y have differing lengths.")
    sumxy <- sum(x|y,na.rm=TRUE)
    if(is.na(sumxy) || !sumxy) return(NA)
    
  Jaccard = sum(x&y,na.rm=TRUE)/sumxy
  Jaccard
}

#' Dice Index
#'
#' @param x binary image or vector
#' @param y binary image or vector
#'
#' @return Dice index
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' x <- rep(c(0,1), 7)
#' y <- rep(c(1,1), 7)
#' dice(x,y)
#'}
dice <- function(x, y){
    jacc=jaccard(x,y)
    2*jacc/(1+jacc)
}

#' Accuracy (Hamming Index)
#'
#' @param x binary image or vector
#' @param y binary image or vector
#'
#' @return Accuracy (Hamming index)
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' x <- rep(c(0,1), 7)
#' y <- rep(c(1,1), 7)
#' hamming(x,y)
#'}
hamming <- function(x, y){
    if (length(x) != length(y)) stop("x and y have differing lengths.")
  n <- sum(!is.na(x)|!is.na(y))
  a <- sum(x == y,na.rm=TRUE)
 a/n
}


#' Covariance function (internal)
#'
#' @param x binary or categorical image or vector
#' @param y binary or categorical image or vector
#' @noRd
#'
#' @return Covariance function
#' @keywords internal
#' \dontrun{
#' x <- rep(1:4,4)
#' y <- c(rep(1:4,3),rep(4,4))
#' sfunc(x,y)
#' }
sfunc <- function(x, y, methodflag = "Cohen"){
    x <- as.vector(x)
    y <- as.vector(y)
    if (methodflag == "Cohen"){
    return(C_Cohen(x,y))
    } else if (methodflag == "Jaccard"){
      return(jaccard(x,y))
    } else if (methodflag == "Dice"){
        return(dice(x,y))
    } else if (methodflag== "hamming"){
        return(hamming(x,y))
               }  else C_AdjRand(x,y)
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
#' @param c1 small normalization constant for the \code{c} function, by default 0.01
#' @param c2 small normalization constant for the \code{s} function, by default 0.01
#' @param method whether to use Cohen's kappa, Jaccard Index, Dice Index, accuracy (Hamming index), or Adjusted Rand Index as
#'     the similarity index. Note Jaccard and Dice should only be used on binary data.
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
binssim <- function(x, y, alpha = 1, beta = 1, gamma = 1, c1 = 0.01, c2 = 0.01, method = "Cohen", ...){
  if (length(x) != length(y)) stop("x and y must be the same size.")
  k = length(unique(c(x,y)))

  methodflag = "Cohen"
  if (method %in% c("AdjRand", "Rand", "rand", "adjrand","R","r","a","A")) methodflag = "AdjRand"
  if (method %in% c("Jaccard", "jaccard", "j", "J")) methodflag = "Jaccard"
  if (method %in% c("Dice","dice","D","d")) methodflag="Dice"
  if (method %in% c("Accuracy","accuracy","Hamming","hamming","H","h")) methodflag="hamming"

  (meansfunc(x, y, c1)^alpha)*(cfunc(x=x, y=y, c2=c2, k = k, ...)^beta)*(sfunc(x, y, methodflag)^gamma)
}

#' Categorical SSIM Components
#'
#' @param x binary or categorical image
#' @param y binary or categorical image
#' @param method whether to use Cohen's kappa, Jaccard, or Adjusted Rand Index as
#'     the similarity index. Note Jaccard should only be used on binary data.
#' @param c1 constant for the means function
#' @param c2 constant for the chrominance function
#' @param ... constants can be passed to the internal functions
#' @noRd
#' @return the three components of the Categorical SSIM.
#' @keywords internal
#'
ssimcomponents <- function(x, y, k, method = "Cohen", c1 = 0.01, c2 = 0.01, ...){
  #k = length(levels)
  #levels <- levels(factor(c(x,y)))
  methodflag = "Cohen"
  if (method %in% c("AdjRand", "Rand", "rand", "adjrand","R","r","a","A")) methodflag = "AdjRand"
  if (method %in% c("Jaccard", "jaccard", "j", "J")) methodflag = "Jaccard"
  if (method %in% c("Dice","dice","D","d")) methodflag="Dice"
  if (method %in% c("Accuracy","accuracy","Hamming","hamming","H","h")) methodflag="hamming"

  c((meansfunc(x, y, c1)),(cfunc(x=x, y=y, c2=c2, k = k, ...)),(sfunc(x, y, methodflag)))
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
#' Cuts the image into a grid of \eqn{2 \times 2}{2x2} squares and selects the mode of each
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
    if(is.null(dims)) stop("x is 1-dimensional")
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
#' @param x an \eqn{n \times m \times q}{n x m x q} binary or categorical image
#'
#' @return a an \eqn{n/2 \times m/2 \times q}{n/2 x m/2 x q} binary or categorical image
#'
#' @keywords internal
#'
#' @noRd
downsample_3d_slice <- function(x){
    dims = dim(x)
    if (is.null(dims)) stop("x is 1-dimensional")
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
#' @param x an \eqn{n \times m \times q}{n x m x q} binary or categorical image
#'
#' @return  an \eqn{n/2 \times m/2 \times q/2}{n/2 x m/2 x q/2} binary or categorical image
#'
#' @keywords internal
#'
#' @noRd
downsample_3d_cube <- function(x){
    dims = dim(x)
    if (is.null(dims)) stop("x is 1-dimensional")
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
#' @param window window size, by default 11
#' @param method whether to use Cohen's kappa, Jaccard Index, accuracy (Hamming index),  or Adjusted Rand Index as
#'     the similarity index. Note Jaccard should only be used on binary data.
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
catssim_2d <- function(x,y, window = 11, method = "Cohen", ...){
    if (is.null(dim(x))) stop("x is 1-dimensional")
    if (is.null(dim(y))) stop("y is 1-dimensional")
    if (length(dim(x)) != length(dim(y)))  stop('x and y have nonconformable dimensions.')
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions.')
  #levels <- levels(factor(c(x,y)))
  k = length(unique(c(x,y)))
  dims = dim(x)
  nrow = dims[1]
  ncol = dims[2]
  if (any(dims < (window - 1))) return(ssimcomponents(x=(x), y=(y), k=k, method=method, ...))
  resultmatrix = array(0, c(nrow-window+1,ncol-window+1,3))# c(0,0,0)

  for (i in 1:(nrow - (window - 1))) {
    for (j in 1:(ncol - (window - 1))) {
      #place = j + (ncol - (window - 1))*(i - 1)
      subx = x[i:(i + (window - 1)), j:(j + (window - 1))]
      suby = y[i:(i + (window - 1)), j:(j + (window - 1))]

      resultmatrix[i,j,] = ssimcomponents(x=subx, y=suby, k=k, method=method, ...)
    }
  }
  resultmatrix[resultmatrix < 0.0] <- 0.0

  #resultmatrix / ((nrow - (window - 1)) * (ncol - (window - 1)))
  colMeans(resultmatrix,na.rm=TRUE,dims=2)
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
#' @param window window size, by default 11.
#' @param method whether to use Cohen's kappa, Jaccard Index, Dice index,  accuracy (Hamming index),  or Adjusted Rand Index as
#'     the similarity index. Note Jaccard and Dice should only be used on binary data.
#' @param ... additional constants can be passed to internal functions.
#'
#' @return a value less than 1 indicating the similarity between the images.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' x <- matrix(sample(0:3, 128^2, replace = TRUE), nrow=128)
#' y <- x
#' for (i in 1:128) y[i, i] = 0
#' for (i in 1:127) y[i, i+1] = 0
#' catmssim_2d(x,y)
#' # now using a different similarity score (Jaccard Index)
#' catmssim_2d(x,y, method = "Jaccard")
#' # now using a different similarity score (Adjusted Rand Index)
#'  catmssim_2d(x,y, method = "Rand")
#'
catmssim_2d <- function(x, y, weights = c(0.0448, 0.2856, 0.3001, 0.2363, 0.1333), window = 11,
                        method = "Cohen", ...){
    # the weights are from the original MS-SSIM program
    if (is.null(dim(x))) stop("x is 1-dimensional")
    if (is.null(dim(y))) stop("y is 1-dimensional")
    if (length(dim(x)) != length(dim(y)))  stop('x and y have nonconformable dimensions.')
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions.')
  levels = length(weights)
  mindim <- min(dim(x))
  if (mindim < window) {
    warning("Minimum dimension should be greater than window size. Using only one level.")
    return(binssim(x=x,y=y, method=method,...))
  }

  if (mindim < (2^(levels - 1))*window) {
    levels = min(c(floor(log2(mindim/window) + 1),levels))
    warning("Truncating levels because of minimum dimension.")
  }
  weights = weights[1:levels]
  results = matrix(0, nrow = levels, ncol = 3)
  results[1,] = catssim_2d(x=x, y=y, window=window, method=method,...)

  if ( levels > 1) {
    for (i in 2:levels) {
      x = downsample_2d(x)
      y = downsample_2d(y)
      results[i,] = catssim_2d(x=x, y=y, window=window, method=method,...)
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
#' @param window by default 11
#' @param method whether to use Cohen's kappa, Jaccard Index, Dice Index, accuracy (Hamming index),  or Adjusted Rand Index as
#'     the similarity index. Note Jaccard or Dice should only be used on binary data.
#' @param ...
#'
#' @return SSIM componenets for the cube.
#' @keywords internal
#'
#' @noRd
#'
catssim_3d_slice <- function(x, y, window = 11, method = "Cohen", ...){
    if (is.null(dim(x))) stop("x is 1-dimensional")
    if (is.null(dim(y))) stop("y is 1-dimensional")
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions.')
  dims = dim(x)
  sliceresults = matrix(0, nrow = dims[3], ncol = 3)
  for (i in 1:dims[3]) sliceresults[i,] = catssim_2d(x=x[,,i],y=y[,,i], window=window, method=method,...)
  colMeans(sliceresults)
}

catssim_3d_cube <- function(x, y, window = 5, method = "Cohen", ...){
    if (is.null(dim(x))) stop("x is 1-dimensional")
    if (is.null(dim(y))) stop("y is 1-dimensional")
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions.')
  #levels <- levels(factor(c(x,y)))
  k = length(unique(c(x,y)))
  dims = dim(x)
  if (any(dims < window)) return(ssimcomponents(x=(x), y=(y),k=k, method=method, ...))

  cuberesults = array(0,c(dims-window+1,3))# c(0,0,0)
  for (i in 1:(dims[1] - (window-1))) {
    for (j in 1:(dims[2] - (window-1))) {
      for (k in 1:(dims[3] - (window-1))) {
        subx = x[i:(i + (window-1)), j:(j + (window-1)), k:(k + (window-1))]
        suby = y[i:(i + (window-1)), j:(j + (window-1)), k:(k + (window-1))]

        cuberesults[i,j,k,] = ssimcomponents(x=subx, y=suby, k=k, method=method, ...)
      }
    }
  }
    colMeans(cuberesults,na.rm=TRUE,dims=3)
#  (cuberesults / prod(dims - (window-1)) )
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
#' @param window window size, by default 11.
#' @param method whether to use Cohen's kappa, Jaccard Index, Dice Index,  accuracy (Hamming index),  or Adjusted Rand Index as
#'     the similarity index. Note Jaccard or Dice should only be used on binary data.
#' @param ... additional constants can be passed to internal functions.
#'
#' @return a value less than 1 indicating the similarity between the images.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' dim = 16
#' x <- array(sample(0:4, dim^3, replace = TRUE), dim = c(dim,dim,dim))
#' y <- x
#' for (j in 1:dim){
#' for (i in 1:dim) y[i, i, j] = 0
#' for (i in 1:(dim-1)) y[i, i+1, j] = 0
#' }
#' catmssim_3d_slice(x,y, weights = c(.75,.25))
#' # Now using a different similarity score
#' catmssim_3d_slice(x,y, weights = c(.75,.25), method = "Jaccard")
#' # And using the last possible similarity score
#' catmssim_3d_slice(x,y, weights = c(.75,.25), method = "Rand")
catmssim_3d_slice <- function(x, y, weights = c(0.0448, 0.2856, 0.3001, 0.2363, 0.1333),
                              window = 11, method = "Cohen", ...){
  # the weights are from the original MS-SSIM program
    if (is.null(dim(x))) stop("x is 1-dimensional")
    if (is.null(dim(y))) stop("y is 1-dimensional")
    if (length(dim(x)) != length(dim(y)))  stop('x and y have nonconformable dimensions.')
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

  results[1,] = catssim_3d_slice(x=x, y=y, window=window, method=method,...)

  if ( levels > 1) {
    for (j in 2:levels) {
      x = downsample_3d_slice(x)
      y = downsample_3d_slice(y)
      results[j,] = catssim_3d_slice(x=x, y=y, window=window, method=method,...)
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
#' This computes a 3D measure based on \eqn{4 \times 4 \times 4}{4x4x4}
#' windows by default with 5 levels of downsampling.
#'
#' @param x a binary or categorical image
#' @param y a binary or categorical image
#' @param weights a vector of weights for the different scales. By default,
#'     five different scales are used.
#' @param window size of window, by default 5
#' @param method whether to use Cohen's kappa, Jaccard Index, Dice Index,  accuracy (Hamming index), or Adjusted Rand Index as
#'     the similarity index. Note Jaccard or Dice should only be used on binary data.
#' @param ... additional constants can be passed to internal functions.
#'
#' @return a value less than 1 indicating the similarity between the images.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' dim = 16
#' x <- array(sample(0:4, dim^3, replace = TRUE), dim = c(dim,dim,dim))
#' y <- x
#' for (j in 1:dim){
#' for (i in 1:dim) y[i, i, j] = 0
#' for (i in 1:(dim-1)) y[i, i+1, j] = 0
#' }
#' catmssim_3d_cube(x,y, weights = c(.75,.25))
#' # Now using a different similarity score
#' catmssim_3d_cube(x,y, weights = c(.75,.25), method = "Jaccard")
#' # And using the last possible similarity score
#' catmssim_3d_cube(x,y, weights = c(.75,.25), method = "Rand")


catmssim_3d_cube <- function(x, y, weights = c(0.0448, 0.2856, 0.3001, 0.2363, 0.1333), window = 5,
                             method = "Cohen", ...){
   # the weights are from the original MS-SSIM program
    if (is.null(dim(x))) stop("x is 1-dimensional")
    if (is.null(dim(y))) stop("y is 1-dimensional")
    if (length(dim(x)) != length(dim(y)))  stop('x and y have nonconformable dimensions.')
  if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions.')
  levels = length(weights)
  dims = dim(x)
  if (length(dims) < 3) stop('x and y are not 3-dimensional.')
  mindim <- min(dim(x))
  if (mindim < 2*window) {
      warning("Minimum dimension must be greater than 2 * window.")
      if(window > mindim) windox = mindim
    return(catssim_3d_cube(x=x,y=y,window=window,method=method,...))
  }
  if (mindim < (2^(levels - 1))*window) {
    levels = min(c(floor(log2(mindim/window) + 1),levels))
    warning("Truncating levels because of minimum dimension.")
  }

  weights = weights[1:levels]
  results = matrix(0, nrow = levels, ncol = 3)

  results[1,] = catssim_3d_cube(x=x,y=y,window=window, method=method,...)

  if ( levels > 1) {
    for (j in 2:levels) {
      x = downsample_3d_cube(x)
      y = downsample_3d_cube(y)
      results[j,] = catssim_3d_cube(x=x,y=y,window=window, method=method,...)
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
#' \eqn{-10 log_{10}(MSE)}. The adjusted Rand index, Jaccard Index, and Cohen's Kappa are used as a measure of
#' the similarity of the structure of the two images. A small constant is added to the numerator
#' and denominator of the Adjusted Rand index to ensure stability, as it is possible to have a zero
#' denominator. The PSNR can be infinite if the error rate is 0. The Jaccard index and accuracy can
#' account for NA values, but the Adjusted Rand and Cohen's Kappa cannot.
#'
#'
#' @param x a numeric or factor vector or image
#' @param y a numeric or factor vector or image
#'
#' @return The accuracy, Jaccard index, the Adjusted Rand Index, the PSNR, and Cohen's Kappa. Note:
#'     The Jaccard index will not make sense if this is not binary.
#'
#' @references Lawrence Hubert and Phipps Arabie (1985).
#' "Comparing partitions". Journal of Classification. 2 (1): 193–218. \doi{10.1007/BF01908075}
#'
#'  W. M. Rand (1971). "Objective criteria for the evaluation of clustering methods".
#'  Journal of the American Statistical Association. American Statistical Association. 66 (336): 846–850.
#'  \doi{10.2307/2284239}
#'
#'  Cohen, Jacob (1960). "A coefficient of agreement for nominal scales".
#'   Educational and Psychological Measurement. 20 (1): 37–46. \doi{10.1177/001316446002000104}
#'
#' @export
#'
#' @examples
#'
#' x <- rep(0:5, 5)
#' y <- c(rep(0:5, 4),rep(0,6))
#' AdjRandIndex(x, y)
AdjRandIndex <- function(x,y){
    if (length(x) != length(y)) stop("x and y have differing lengths.")
    if(!all(!is.na(x),!is.na(y))) warning("NAs present in x or y, Adjusted Rand and Cohen don't account for NA values.")
    if (length(table(c(x,y))) > 2)
        warning("Jaccard index may not make sense if more than two classes are present.")
  n <- sum(!is.na(x)|!is.na(y))
  a <- sum(x == y,na.rm=TRUE)
  Accuracy <- a/n
  BinJaccard <- jaccard(x,y)
  Cohen <- C_Cohen(x,y)
  x <- as.numeric(x)
  y <- as.numeric(y)
  AdjRand <- C_AdjRand(x,y)
    list(Accuracy = Accuracy,
         Jaccard = BinJaccard,
       AdjRand = AdjRand,
       PSNR = -10 * log10(1-Accuracy),
       Cohen = Cohen)
}


#' Multiscale Categorical Structural Similarity Index Measure
#'
#' The categorical structural similary index measure for 2D or 3D categorical or binary
#' images for multiple scales. The default is to compute over 5 scales.
#' This determines whether this is a 2D or 3D image and applies the appropriate
#' windowing, weighting, and scaling. Additional arguments can be passed.
#' This is a wrapper function for the 2D and 3D functions whose functionality
#' can be accessed through the ... arguments.
#'
#' @param x a binary or categorical image
#' @param y a binary or categorical image
#' @param ... additional arguments, such as method and window, can be passed
#'        as well as arguments for internal functions.
#' @param cube for the 3D method, whether to use the true 3D method (cube)
#'        or compute the metric using 2D slices which are then averaged.
#'        By default, TRUE, which evaluates as a cube. FALSE will treat it as
#'        2D slices.
#' @param weights  a vector of weights for the different scales. By default,
#'     five different scales are used.
#' @return a value less than 1 indicating the similarity between the images.
#' @export
#'
#' @examples
#' set.seed(20181207)
#' dim = 16
#' x <- array(sample(0:4, dim^3, replace = TRUE), dim = c(dim,dim,dim))
#' y <- x
#' for (j in 1:dim){
#' for (i in 1:dim) y[i, i, j] = 0
#' for (i in 1:(dim-1)) y[i, i+1, j] = 0
#' }
#' catsim(x,y, weights = c(.75,.25))
#' # Now using a different similarity score
#' catsim(x,y, weights = c(.75,.25), method = "Jaccard")
#' # And using the last possible similarity score
#' catsim(x,y, weights = c(.75,.25), method = "Rand")
#' # with the slice method:
#' catsim(x,y, weights = c(.75,.25), cube = FALSE)
catsim <- function(x,y,...,cube = TRUE, weights =  c(0.0448, 0.2856, 0.3001, 0.2363, 0.1333)){
    if (is.null(dim(x))) stop("x is 1-dimensional")
    if (is.null(dim(y))) stop("y is 1-dimensional")
    if (length(dim(x)) != length(dim(y)))  stop('x and y have nonconformable dimensions.')
    if (any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions.')
    dims = dim(x)
    if(length(dims)==2) catmssim_2d(x,y,weights=weights,...)
    else if (cube) catmssim_3d_cube(x,y,weights=weights,...)
    else catmssim_3d_slice(x,y,weights=weights,...)
}
