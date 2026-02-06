# Find phi parameter for the Exponential spatial auto-correlation function

Function designed to find the phi parameter such that the correlation
between points within a given distance `d` is at most a given value.

## Usage

``` r
find_phi(
  d,
  nu,
  kappa,
  mu2,
  family = "matern",
  range = c(1e-04, 1000),
  cut = 0.05
)
```

## Arguments

- d:

  maximum distance for spatial dependence equal to `cut`.

- nu:

  smoothness parameter associated with the Matern cov. function.

- kappa:

  one of the smoothness parameters associated with the Generalized
  Wendland covariance function

- mu2:

  one of the smoothness parameters associated with the Generalized
  Wendland covariance function

- family:

  covariance function family, the options are
  `c("matern", "gw", "cs", "spher", "pexp", "gaussian")`.

- range:

  Minimum and maximum distance to be considered. The default is
  `range = c(1e-04, 1000)`.

- cut:

  desired spatial correlation at a distance `d`, the default is
  `cut = .05`.

## Value

a `numeric` value indicating the range parameter such that the spatial
correlation between two points at distance `d` is `cut`.
