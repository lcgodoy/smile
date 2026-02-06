# Cubic spline covariance function (scalar)

Computing the Spherical covariance function for a scalar distance.

## Usage

``` r
single_cs(d, sigsq, phi)
```

## Arguments

- d:

  a scalar representing the distance on which it is desired to evaluate
  the covariance function.

- sigsq:

  the \\\sigma^2\\ parameter from the Spherical covariance. function.

- phi:

  the \\\phi\\ parameter from the Spherical covariance function,
  controls the range of the spatial dependence.

## Value

a scalar representing the (gaussian) covariance between two observations
`d` apart of each other.

## See also

[`single_exp`](https://lcgodoy.me/smile/reference/single-matern.md),
[`single_matern`](https://lcgodoy.me/smile/reference/single-matern.md),
[`single_matern3`](https://lcgodoy.me/smile/reference/single-matern.md),
[`single_matern5`](https://lcgodoy.me/smile/reference/single-matern.md),
[`mat_cov`](https://lcgodoy.me/smile/reference/mat_cov.md)
