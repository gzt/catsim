
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Coverage
status](https://codecov.io/gh/gzt/catsim/branch/master/graph/badge.svg)](https://app.codecov.io/github/gzt/catsim?branch=master)
[![R build
status](https://github.com/gzt/catsim/workflows/R-CMD-check/badge.svg)](https://github.com/gzt/catsim/actions)
<!-- badges: end -->

# `catsim`: a Categorical Image Similarity Index

The goal of `catsim` is to provide a similarity measure for binary or
categorical images in either 2D or 3D similar to the [MS-SSIM
index](https://en.wikipedia.org/wiki/Structural_similarity) for color
images. Suppose you have a ground truth segmentation of some image that
has been segmented into regions - perhaps a brain scan with different
types of tissues or a map with different types of terrain - and a
segmentation produced by some classification method. Comparing the two
pixel-by pixel (or voxel-by-voxel) might work well, but a method that
captures structural similarities might work better for your purposes.
MS-SSIM is an image comparison metric that tries to match the assessment
of the human visual system by considering structural similarities across
multiple scales. CatSIM applies a similar logic in the case of 2-D and
3-D binary and multicategory images, such as might be found in image
segmentation or classification problems.

## Installation

You can install the released version of catsim from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("catsim")
#### or the dev version with:
#devtools::install_github("gzt/catsim")
```

## Usage

If you have two images, `x` and `y`, the simplest method of comparing
them is:

``` r
library(catsim)
set.seed(20200505)
x <- besag
y <- x
y[10:20,10:20] <- 1
catsim(x, y, levels = 3)
#> [1] 0
```

By default, this performs 5 levels of downsampling and uses Cohen’s
kappa as the local similarity metric on `11 x 11` windows for a
2-dimensional image and `5 x 5 x 5` windows for a 3-D image. Those can
be adjusted using the `levels`, `method`, and `window` arguments.

Please note that the `catsim` project is released with a [Contributor
Code of Conduct](https://gzt.github.io/catsim/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.
