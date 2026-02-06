# Powered Exponential covariance function (scalar)

Computing the Powered Exponential covariance function for a scalar
distance.

## Usage

``` r
single_pexp(d, sigsq, phi, nu)
```

## Arguments

- d:

  a scalar representing the distance on which it is desired to evaluate
  the covariance function.

- sigsq:

  the \\\sigma^2\\ parameter from the Exponential covariance function.

- phi:

  the \\\phi\\ parameter from the Exponential covariance function,
  controls the range of the spatial dependence.

- nu:

  the \\\nu \in (0, 2\]\\ parameter representing the "power"

## Value

a scalar representing the (exponential) covariance between two
observations `d` apart of each other.

## See also

[`single_exp`](https://lcgodoy.me/smile/reference/single-matern.md),
[`single_matern`](https://lcgodoy.me/smile/reference/single-matern.md),
[`single_matern3`](https://lcgodoy.me/smile/reference/single-matern.md),
[`single_matern5`](https://lcgodoy.me/smile/reference/single-matern.md),
[`mat_cov`](https://lcgodoy.me/smile/reference/mat_cov.md)
