###   utils.R
###   CatSIM : Utility Functions for CatSIM
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

##' Method Parser
##'
##' Parses input of method to a standardized name
##' @param method The method used as a similarity metric.
##'    Certain abbreviations work.
##'    `Cohen`, `cohen`, `C`, `c`, `Kappa` and
##'    `kappa` yield Cohen's kappa.
##'    `AdjRand`, `adjrand`, `Adj`, `adj`, `a`,
##'    `A`, `ARI`, and `ari` yield the adjusted Rand index.
##'    `Rand`, `rand`, `r`, and `R` yield the Rand index.
##'    `Jaccard`, `jaccard`, `j`, and `J` yield the
##'    Jaccard index.
##'    `Dice`, `dice`, `D`, and `d` yield the Dice index.
##'    `Accuracy`, `accuracy`, `Hamming`, `hamming`,
##'    `H`, and `h` yield the accuracy.
##'    `NMI`, `MI`, `mutual`, `information`, `nmi`,
##'    `mi` yield the normalized mutual information.
##'    `AMI` and `ami` yield the adjusted mutual information.
##' @return the name of the similarity metric.
##' @keywords internal
##' @noRd
method_parser <- function(method) {
  methodflag <- NULL
  if (method %in% c("Cohen", "cohen", "C", "c", "kappa", "Kappa")) {
    methodflag <- c_cohen
  } # "Cohen"
  if (method %in% c(
    "AdjRand", "adjrand", "Adj", "adj",
    "a", "A", "ARI", "ari"
  )) {
    methodflag <- c_adj_rand
  } # "AdjRand"
  if (method %in% c("Rand", "rand", "r", "R")) {
    methodflag <- c_rand
  } # "Rand"
  if (method %in% c("Jaccard", "jaccard", "j", "J")) {
    methodflag <- jaccard
  } # "Jaccard"
  if (method %in% c("Dice", "dice", "D", "d")) {
    methodflag <- dice
  } # "Dice"
  if (method %in% c("Accuracy", "accuracy", "acc", "Hamming", "hamming", "H", "h")) {
    methodflag <- hamming
  } # "hamming"
  if (method %in% c("NMI", "MI", "mutual", "information", "nmi", "mi")) {
    methodflag <- c_nmi
  } # normalized mutual information
  if (method %in% c("AMI", "ami")) {
    methodflag <- c_ami
  } # adjusted mutual information
  if (is.null(methodflag)) stop("Error: invalid method")
  methodflag
}

##' Level Parser
##'
##' Parses the weights and levels.
##' @param weights the vector of weights
##' @param levels the levels
##' @keywords internal
##' @noRd
##' @return the parsed weights
level_parser <- function(weights, levels) {
  if (!is.null(weights) && !is.null(levels)) {
    if (levels > length(weights)) {
      stop("Inconsistent weight and levels specified.")
    }
    weights <- weights[1:levels]
  }
  if (is.null(weights) && is.null(levels)) {
    levels <- 5
  }
  if (is.null(weights)) {
    weights <- rep(1, levels) / levels
  }
  weights
}

##' Parser for ... arguments.
##'
##' Parses the ... arguments to the defaults for internal functions
##' if missing.
##' @title Dots Parser
##' @param ... arguments passed
##' @return a list with defaults filled in for missing arguments
##' @keywords internal
##' @noRd
dots_parser <- function(...) {
  dotlist <- list(...)
  ## list of args: with defaults
  ## alpha: 1; beta: 1; gamma: 1
  ## c1: 0.01
  ## c2: 0.01
  ## sqrtgini: TRUE
  ## random or rand: FALSE
  if (is.null(dotlist[["c1"]])) dotlist[["c1"]] <- 0.01
  if (is.null(dotlist[["c2"]])) dotlist[["c2"]] <- 0.01
  if (is.null(dotlist[["sqrtgini"]])) dotlist[["sqrtgini"]] <- TRUE
  if (is.null(dotlist[["random"]])) dotlist[["random"]] <- FALSE
  dotlist
}
