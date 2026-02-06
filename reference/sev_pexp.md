# Calculate Smallest Eigenvalue for Power Exponential Correlation Matrices

This function computes the smallest eigenvalue of a correlation matrix
derived from the power exponential correlation function. It evaluates
this across a grid of values for the power parameter (`nu`) and the
practical range parameter (`rho`), based on a provided distance matrix.

## Usage

``` r
sev_pexp(range_nu, range_rho, grid_len = 50, dmat)
```

## Arguments

- range_nu:

  A numeric vector of length 2, specifying the minimum and maximum
  values for the power parameter `nu`. `nu` typically ranges between 0
  and 2 (e.g., `nu = 1` for exponential, `nu = 2` for Gaussian).

- range_rho:

  A numeric vector of length 2, specifying the minimum and maximum
  values for the practical range parameter `rho`. `rho` must be
  positive.

- grid_len:

  An integer specifying the number of points to create for both `nu` and
  `rho` sequences. The total number of grid combinations will be
  `grid_len^2`. Default is 50.

- dmat:

  A numeric matrix representing the distance matrix between locations.
  The distances should be non-negative.

## Value

A tibble with three columns:

- rho:

  The practical range parameter value.

- nu:

  The power parameter value.

- lambda:

  The smallest eigenvalue of the power exponential correlation matrix
  corresponding to the `rho` and `nu` pair.

## Details

The practical range `rho` is defined here as the distance at which the
correlation is 0.1. The internal scale parameter `phi` is calculated as
`phi = rho / (log(10)^(1/nu))`. The power exponential correlation
function is assumed to be of the form C(h) = exp(-(h/phi)^nu), where h
is distance. The function `smile:::pexp_cov` is used internally to
compute the covariance/correlation matrix with a sill of 1.

The function first creates a grid of `nu` and `rho` parameters. For each
pair of (`rho`, `nu`) in the grid: 1. It calculates the scale parameter
`phi` for the power exponential correlation function, where
`phi = rho / (log(10)^(1/nu))`. This definition implies that the
correlation is 0.1 at the distance `rho`. 2. It computes the power
exponential correlation matrix using
`smile:::pexp_cov(dists = dmat, sill = 1, range = phi, smooth = nu)`.
Note the use of an internal function from the `smile` package. 3. It
calculates the eigenvalues of this correlation matrix. 4. The minimum
eigenvalue is extracted. The final output is a tibble containing all
parameter combinations and their corresponding minimum eigenvalues.
