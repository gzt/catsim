###   data.R
###   CatSIM : Data Loading for CatSIM
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

#' A hand-constructed image from Besag (1986)
#'
#' An \eqn{100 \times 88} matrix representing a two-color
#' hand-drawn scene designed specifically to contain some
#' awkward features for an image reconstruction method
#' evaluated in the paper.
#'
#' @format an \eqn{100 \times 88} matrix with entries `1` and `2`
#' denoting the color of the corresponding pixels. The example code will
#' produce the image as it is in the original paper. To use as a `0-1` binary
#' dataset, either use `besag - 1` or `besag %% 2`.
#' @docType data
#' @usage data(besag)
#' @references
#' J. Besag, “On the statistical analysis of dirty pictures,”
#' Journal of the Royal Statistical Society: Series B (Methodological),
#' vol. 48, no. 3, pp. 259–279, 1986. \doi{10.1111/j.2517-6161.1986.tb01412.x}
#' @keywords datasets
#' @examples
#' image(besag[, 88:1])
"besag"

#' An example fMRI phantom
#'
#' A \eqn{128 \times 128}{128 x 128} activation map for a slice of an fMRI
#' phantom and an anatomical reference.
#'
#' @format a \eqn{128 \times 128 \times 2}{128 x 128 x 2} array with the
#' first slice an activation map for an MRI phantom and the second an
#' anatomical overlay. \code{NA} values are outside the surface. The
#' activation map (\code{hoffmanphantom[,,1]}) is 1 if activated,
#' 0 otherwise. The second layer  (\code{hoffmanphantom[,,2]})
#' indicates the anatomical structure. Approximately 3.8 percent of the
#' pixels are activated in this slice.
#' @docType data
#' @keywords datasets
#' @usage data(hoffmanphantom)
#' @references
#' E. Hoffman, P. Cutler, W. Digby, and J. Mazziotta, “3-D phantom
#' to simulate cerebral blood flow and metabolic images for PET,”
#' Nuclear Science, IEEE Transactions on, vol. 37, pp. 616 – 620, 05
#' 1990.
#'
#'  I. A. Almodóvar-Rivera and R. Maitra, “FAST adaptive smoothed
#' thresholding for improved activation detection in low-signal
#' fMRI,” IEEE Transactions on Medical Imaging, vol. 38, no. 12, pp.
#' 2821–2828, 2019.
#' @examples
#' image(hoffmanphantom[, , 2], col = rev(gray(0:15 / 16))[1:4], axes = FALSE)
#' image(hoffmanphantom[, , 1],
#'   add = TRUE, zlim = c(0.01, 1),
#'   col = c("yellow", "maroon")
#' )
"hoffmanphantom"
