# Matern covariance function (scalar - generic)

Computing the Matern covariance function for a scalar distance, adapted
from `geoR`.

## Usage

``` r
single_exp(d, sigsq, phi)

single_matern3(d, sigsq, phi)

single_matern5(d, sigsq, phi)

single_matern(d, sigsq, phi, nu)

single_gauss(d, sigsq, phi)
```

## Arguments

- d:

  a scalar representing the distance on which it is desired to evaluate
  the covariance function.

- sigsq:

  the \\\sigma^2\\ parameter from the Matern covariance function.

- phi:

  the \\\phi\\ parameter from the Matern covariance function, controls
  the range of the spatial dependence.

- nu:

  the \\\nu\\ parameter from the Matern covariance function, controls
  the differentiability of the process.

## Value

a scalar representing the (matern) covariance between two observations
`d` apart of each other.

## Details

`single_matern3` and `single_matern5` are optimized for when \\\nu\\ is
1.5 or 2.5, respectively. Similarly, `single_exp` and `single_gauss`
represent the cases where \\\nu = 0.5\\ or \\\nu \to \infty\\. In other
words, they are the exponential and Gaussian covariance functions.

## See also

`single_matern3`, `single_matern5` `single_exp`,
[`mat_cov`](https://lcgodoy.me/smile/reference/mat_cov.md)
