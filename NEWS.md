# smile 1.2.0

* New distances between spatial units, in the plane and on the sphere, with
  fast C++ code: `dist_hausdorff()` (Hausdorff distance) and `dist_cr()`
  (center-radii distance). `get_mbc_params()` returns the Chebyshev centers
  and radii. Unit tests compare the results with `sf`.

# smile 1.1.0

* Simplified logic for the log-lik functions

* Add CITATION file

* Keeping only covariance functions valid on $R^n$ (as $n \to \infty$). This
  keeps the code easier to mantain and extend. In addition, allows us to fit
  models on different domains (on the square-root scale).

* Function for empirical validation of the Matern Covariance function

# smile 1.1.0

* Changed how likelihood calculations work (for increased performance)

* Transitioned from `RcppArmadillo` to `RcppEigen`

* Function for empirical validation of the Powered Exponential Covariance
  function

# smile 1.0.5

* Minor bug fixes and documentation enhancement.

* Removed NYC datasets (they are too large).

* Spatial areal interpolation vignette rewritten.

# smile 1.0.4.2

* New dataset included
* Parallel version of `fit_spm2`

# smile 1.0.4.1

* Speeding up `sf_to_spm`.
* Added icon
* Improved documentation
* Added examples

# smile 1.0.4

* Generalized Wendland covariance function
* More flexible prediction function

# smile 1.0.3.1

* bug fixes for prediction

# smile 1.0.3

* fixed Wendland (1) covariance function
* support for sparse covariance matrices through the `Matrix` package

# smile 1.0.2

* added Wendland (1) covariance function
* added "tapered" Matern covariance function

# smile 1.0.1

* Added a `NEWS.md` file to track changes to the package.
* Fixed Spherical covariance function.
* Added Cubic spline covariance function.
