# catsim 0.2.1
* Fix some cosmetic issues in README.Rmd
* Fix a spurious warning in compiled code about an unused variable. 

# catsim 0.2.0
* Last cosmetic cleanups and bumping version number because
this really should have been done for the change in 0.1.5,
as removing deprecated functions needs a bigger version bump.

# catsim 0.1.6
* Clean up final documentation issues in preparation for submission to CRAN.

# catsim 0.1.5
* Remove deprecated `AdjRandIndex` files.

# catsim 0.1.4
* Further refinement of `pickmode(random)` handling
* Further documentation and vignette updates.
* Move to Markdown documentation.

# catsim 0.1.2
* Improve handling of `...` arguments.
* Easier access to `pickmode(random)` argument.
* Included `hoffmanphantom` dataset.

# catsim 0.1.1
* Change deterministic PRNG for downsampling to not assign to parent
  environment to preserve state, made the option to use `sample`
  accessible to advanced users through `...` arguments.
* Fixed some documentation. 

# catsim 0.1.0
* Mostly documentation fixes.
* Version used for paper submission.

# catsim 0.1.0.9040

* Add `besag` image as a data example.
* Changed names of similarity functions to use "_" instead of "camel case".
* Features still being added.

# catsim 0.1.0.9039

* Added `levels` argument to `catsim()`. Weights or levels can be specified. 
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
