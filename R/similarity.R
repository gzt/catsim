###   similarity.R
###   CatSIM : Similarity Measures for CatSIM
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

#' @title Similarity Indices
#' @name rand
#'
#' @description The Rand index, [rand_index], computes the agreement
#' between two different clusterings or partitions of the same  set of objects.
#' The inputs to the function should be binary or categorical and of the same
#' length.
#'
#' @param x,y  a numeric or factor vector or array
#' @param na.rm whether to remove `NA` values. By default, `FALSE`.
#' If `TRUE`, will perform pair-wise deletion.
#'
#' @return the similarity index, which is between 0 and 1 for most of the
#' options. The adjusted Rand and Cohen's kappa can be negative, but are
#' bounded above by 1.
#'
#' @references
#'  W. M. Rand (1971). "Objective criteria for the evaluation of clustering
#'  methods". Journal of the American Statistical Association.
#' American Statistical Association. 66 (336): 846–850.
#'  \doi{10.2307/2284239}
#'
#' Lawrence Hubert and Phipps Arabie (1985).
#' "Comparing partitions". Journal of Classification. 2 (1): 193–218.
#' \doi{10.1007/BF01908075}
#'
#'  Cohen, Jacob (1960). "A coefficient of agreement for nominal scales".
#'   Educational and Psychological Measurement. 20 (1): 37–46.
#' \doi{10.1177/001316446002000104}
#'
#'  Jaccard, Paul (1912). "The distribution of the flora in the alpine zone,”
#' New Phytologist, vol. 11, no. 2, pp. 37–50.
#'   \doi{10.1111/j.1469-8137.1912.tb05611.x}
#'
#' Nguyen Xuan Vinh, Julien Epps, and James Bailey (2010).
#' Information Theoretic Measures for Clusterings Comparison:
#' Variants, Properties, Normalization and Correction for Chance.
#' J. Mach. Learn. Res. 11 (December 2010), 2837–2854.
#' \url{http://www.jmlr.org/papers/v11/vinh10a}
#' @export
#'
#' @examples
#' x <- rep(0:5, 5)
#' y <- c(rep(0:5, 4), rep(0, 6))
#' # Simple Matching, or Accuracy
#' mean(x == y)
#' # Hamming distance
#' sum(x != y)
#' rand_index(x, y)
#' adj_rand(x, y)
#' cohen_kappa(x, y)
#' normalized_mi(x, y)
#' adjusted_mi(x, y)
rand_index <- function(x, y, na.rm = FALSE) {
  if (length(x) != length(y)) stop("x and y have differing lengths.")
  if (na.rm) {
    naxy <- (!is.na(x) & !is.na(y))
    x <- x[naxy]
    y <- y[naxy]
  }
  if (!all(!is.na(x), !is.na(y))) {
    warning("NAs present in x or y,
              the Rand index doesn't account for NA values.")
  }
  x <- as.numeric(x)
  y <- as.numeric(y)
  c_rand(x, y)
}



#' @name Adjusted Rand Index
#'
#' @description The adjusted Rand index,  `adj_rand`,
#' computes a corrected version
#' of the Rand index, adjusting for the probability
#' of chance agreement of clusterings. A small constant is added to the
#' numerator and denominator of the adjusted Rand index to ensure stability
#' when there is a small or 0 denominator, as it is possible to have a zero
#' denominator.
#' @rdname rand
#' @export
adj_rand <- function(x, y, na.rm = FALSE) {
  if (length(x) != length(y)) stop("x and y have differing lengths.")
  if (na.rm) {
    naxy <- (!is.na(x) & !is.na(y))
    x <- x[naxy]
    y <- y[naxy]
  }
  if (!all(!is.na(x), !is.na(y))) {
    warning("NAs present in x or y,
               Adjusted Rand doesn't account for NA values.")
  }
  x <- as.numeric(x)
  y <- as.numeric(y)
  c_adj_rand(x, y)
}

#' @name Cohen's kappa
#'
#' @description Cohen's kappa, `cohen_kappa`,
#' is an inter-rater agreement metric for two raters which
#' corrects for the probability of chance agreement. Note
#' there is a difference here
#' between this measure and the Rand indices and mutual information:
#' those consider the similarities of the groupings of points,
#' while this considers how often the
#' raters agreed on individual points.
#' @rdname rand
#' @export
cohen_kappa <- function(x, y, na.rm = FALSE) {
  if (length(x) != length(y)) stop("x and y have differing lengths.")
  if (na.rm) {
    naxy <- (!is.na(x) & !is.na(y))
    x <- x[naxy]
    y <- y[naxy]
  }
  if (!all(!is.na(x), !is.na(y))) {
    warning("NAs present in x or y, Cohen's Kappa doesn't
               account for NA values.")
  }
  x <- as.numeric(x)
  y <- as.numeric(y)
  c_cohen(x, y)
}

#' @name Normalized Mutual Information
#'
#' @description Like the Rand index, the mutual information
#' computes the agreement between two different clusterings or
#' partitions of the same set of objects. If \eqn{H(X)} is the
#' entropy of some probability distribution \eqn{X}, then
#' the mutual information of two distributions is
#' \eqn{I(X;Y) = -H(X,Y) +H(X) + H(Y)}.
#' The normalized mutual information, `normalized_mi`, is defined here as:
#' \eqn{2I(X;Y)/(H(X)+H(Y)),} but is set to be 0 if both H(X) and H(Y) are 0.
#' @rdname rand
#' @export
normalized_mi <- function(x, y, na.rm = FALSE) {
  if (length(x) != length(y)) stop("x and y have differing lengths.")
  if (na.rm) {
    naxy <- (!is.na(x) & !is.na(y))
    x <- x[naxy]
    y <- y[naxy]
  }
  if (!all(!is.na(x), !is.na(y))) {
    warning("NAs present in x or y, normalized mutual
               information doesn't account for NA values.")
  }
  x <- as.numeric(x)
  y <- as.numeric(y)
  c_nmi(x, y)
}

#' @name Adjusted Mutual Information
#'
#' @description The adjusted mutual information, `adjusted_mi`,
#' is a correction of the mutual information to account
#' for the probability of chance agreement in a manner similar to the
#' adjusted Rand index
#' or Cohen's kappa.
#' @rdname rand
#' @export
adjusted_mi <- function(x, y, na.rm = FALSE) {
  if (length(x) != length(y)) stop("x and y have differing lengths.")
  if (na.rm) {
    naxy <- (!is.na(x) & !is.na(y))
    x <- x[naxy]
    y <- y[naxy]
  }
  if (!all(!is.na(x), !is.na(y))) {
    warning("NAs present in x or y, adjusted mutual
               information doesn't account for NA values.")
  }
  x <- as.numeric(x)
  y <- as.numeric(y)
  c_ami(x, y)
}
