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

#include <testthat.h>
#include <Rcpp.h>
using namespace Rcpp;
#include "internalfunction.h"



context("Simple function checks"){
  test_that("Diversity measures correct"){
    Rcpp::NumericVector x(2);
    Rcpp::NumericVector y(2);
    for(int i = 0; i < 2; i++) x[i] = 0;
    y[0] = 0;
    y[1] = 1;
    expect_true(std::abs(C_gini(x))< 1e-5);
    expect_true(std::abs(C_gini(y)-.5) < 1e-5);
    expect_true(std::abs(C_ginicorr(x, 2))< 1e-5);
    expect_true(std::abs(C_ginicorr(y, 2)-1.0)< 1e-5);
    expect_true(std::abs(C_cfunc(x,x,.001,2,TRUE)-1.0) < 1e-5);
    expect_true(std::abs(C_cfunc(x,y,0.0,2,TRUE)) < 1e-5);


  }
  test_that("Difference measures correct"){
    Rcpp::NumericVector x(2);
    Rcpp::NumericVector y(2);
    for(int i = 0; i < 2; i++) x[i] = 0;
    y[0] = 0;
    y[1] = 1;
    expect_true(std::abs(C_Cohen(x,x)-1.0) < 1e-6);
    expect_true(std::abs(C_Cohen(y,y)-1.0) < 1e-6);
    expect_true(std::abs(C_AdjRand(x,x)-1.0) < 1e-3);
    expect_true(std::abs(C_AdjRand(y,y)-1.0) < 1e-3);
    expect_true(std::abs(C_Cohen(x,y)) < 1e-6);
    expect_true(std::abs(C_Cohen(x,y)-C_Cohen(y,x)) < 1e-6);
    expect_true(std::abs(C_AdjRand(x,y)) < 1e-2);
    expect_true(std::abs(C_AdjRand(x,y) - C_AdjRand(y,x)) < 1e-6);
    expect_true(std::abs(C_AdjRand(x,x) - C_AdjRand(y,y)) < 1e-6);

  }
}
