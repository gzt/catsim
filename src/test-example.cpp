/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.

#include <Rcpp.h>
#include <testthat.h>
using namespace Rcpp;
#include "internalfunction.h"

context("Simple function checks") {
  test_that("Diversity measures correct") {
    Rcpp::NumericVector x(2);
    Rcpp::NumericVector y(2);
    for (int i = 0; i < 2; i++) x[i] = 0;
    y[0] = 0;
    y[1] = 1;
    expect_true(std::abs(c_gini(x)) < 1e-5);
    expect_true(std::abs(c_gini(y) - .5) < 1e-5);
    expect_true(std::abs(c_ginicorr(x, 2)) < 1e-5);
    expect_true(std::abs(c_ginicorr(y, 2) - 1.0) < 1e-5);
    expect_true(std::abs(c_cfunc(x, x, .001, 2, TRUE) - 1.0) < 1e-5);
    expect_true(std::abs(c_cfunc(x, y, 0.0, 2, TRUE)) < 1e-5);
    expect_true(std::abs(c_cfunc(x, x, .001, 2, FALSE) - 1.0) < 1e-5);
    expect_true(std::abs(c_cfunc(x, y, 0.0, 2, FALSE)) < 1e-5);
  }
  test_that("Difference measures correct") {
    Rcpp::NumericVector x(2);
    Rcpp::NumericVector y(2);
    for (int i = 0; i < 2; i++) x[i] = 0;
    y[0] = 0;
    y[1] = 1;
    expect_true(std::abs(c_cohen(x, x) - 1.0) < 1e-6);
    expect_true(std::abs(c_cohen(y, y) - 1.0) < 1e-6);
    expect_true(std::abs(c_adj_rand(x, x) - 1.0) < 1e-3);
    expect_true(std::abs(c_adj_rand(y, y) - 1.0) < 1e-3);
    expect_true(std::abs(c_cohen(x, y)) < 1e-6);
    expect_true(std::abs(c_cohen(x, y) - c_cohen(y, x)) < 1e-6);
    expect_true(std::abs(c_adj_rand(x, y)) < 1e-2);
    expect_true(std::abs(c_adj_rand(x, y) - c_adj_rand(y, x)) < 1e-6);
    expect_true(std::abs(c_adj_rand(x, x) - c_adj_rand(y, y)) < 1e-6);
    expect_true(std::abs(c_rand(x, y)) < 1e-6);
    expect_true(std::abs(c_rand(x, y) - c_rand(y, x)) < 1e-6);
    expect_true(std::abs(c_rand(x, x) - c_rand(y, y)) < 1e-6);
    expect_true(std::abs(c_nmi(x, y)) < 1e-5);
    expect_true(std::abs(c_nmi(y, y) - 1.0) < 1e-5);
    expect_true(std::abs(c_ami(x, y)) < 1e-5);
    expect_true(std::abs(c_ami(y, y) - 1.0) < 1e-5);
  }
}
