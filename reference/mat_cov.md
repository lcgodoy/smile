# Matern covariance function for a given distance matrix.

Computing the Matern covariance function for a matrix of distances.

Computing the Matern covariance function between polygons.

## Usage

``` r
mat_cov(dists, sigsq, phi, nu)

comp_mat_cov(cross_dists, n, n2, sigsq, phi, nu)

pexp_cov(dists, sigsq, phi, nu)

comp_pexp_cov(cross_dists, n, n2, sigsq, phi, nu)

gauss_cov(dists, sigsq, phi)

comp_gauss_cov(cross_dists, n, n2, sigsq, phi)

spher_cov(dists, sigsq, phi)

comp_spher_cov(cross_dists, n, n2, sigsq, phi)

cs_cov(dists, sigsq, phi)

comp_cs_cov(cross_dists, n, n2, sigsq, phi)

gw_cov(dists, sigsq, phi, kappa, mu)

comp_gw_cov(cross_dists, n, n2, sigsq, phi, kappa, mu)
```

## Arguments

- dists:

  a numeric matrix representing the distance between spatial entities.

- sigsq:

  the \\\sigma^2\\ parameter from the Matern covariance function.

- phi:

  the \\\phi\\ parameter from the Matern covariance function, controls
  the range of the spatial dependence.

- nu:

  the \\\nu\\ parameter from the Matern covariance function, controls
  the differentiability of the process.

- cross_dists:

  a `list` such that each position contains the cross distances between
  points within different polygons.

## Value

The matern covariance function (for a stationary and isotropic process)
associated with the provided distances (`dists`) and the given set of
parameters.

The matern covariance matrix associated with a set of polygons.

## See also

[`single_exp`](https://lcgodoy.me/smile/reference/single-matern.md),
[`single_matern`](https://lcgodoy.me/smile/reference/single-matern.md)

[`single_exp`](https://lcgodoy.me/smile/reference/single-matern.md),
[`single_matern`](https://lcgodoy.me/smile/reference/single-matern.md),
`mat_cov`
