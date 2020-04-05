#' A hand-constructed image from Besag (1986)
#'
#' An \eqn{100 \times 88} matrix representing a two-color
#' hand-drawn scene designed specifically to contain some
#' awkward features for an image reconstruction method
#' evaluated in the paper.
#'
#' @format an \eqn{100 \times 88} matrix with entries `1` and `2`
#' denoting the color of the corresponding pixels. The example code will
#' produce the image as it is in the original paper.
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
