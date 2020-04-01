# catsim 0.1.0.9040

* Add `besag` image as a data example.
* Changed names of similarity functions to use "_" instead of "camel case".
* Features still being added.

# catsim 0.1.0.9039

* Added `levels` argument to `catsim()`. Weights or levels can be specifed. 
  `levels` alone will define `weights = rep(1, levels)/levels`. If `weights` 
  and `levels` are both specified, `weights = weights[1:levels]`. 
* Added documentation.

# catsim 0.1.0.9038

* Added normalized mutual information and adjusted mutual information
  as similarity metrics and into AdjRandIndex function, which really 
  should be renamed.
* Split out the various similarity measures and exported them. Deprecated
  `AdjRandIndex()` function.

# catsim 0.1.0.9037

* Added a `NEWS.md` file to track changes to the package.
* Added more tests, documentation.
* The package is mostly complete in its features. 
