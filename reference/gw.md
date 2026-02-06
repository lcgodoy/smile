# Matern Generalized Wendland (GW) covariance function (scalar - generic)

adapted from Bevilacqua et al. 2019.

## Usage

``` r
single_gw0(d, sigsq, phi, mu)

single_gw1(d, sigsq, phi, mu)

single_gw2(d, sigsq, phi, mu)

single_gw3(d, sigsq, phi, mu)

single_gw(d, sigsq, phi, kappa, mu)
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

- mu:

  a parameter that controls the smoothness of the covariance function.
  Note that, \\\mu \geq 1\\.

- kappa:

  \\\kappa \in \\0, \ldots, 3 \\\\.

## Value

a scalar representing the GW covariance between two observations `d`
apart of each other.
