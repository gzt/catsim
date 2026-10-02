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
  
  }
}
