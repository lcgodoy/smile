# Evaluate log-lik

Evaluate the log-likelihood for a given set of parameters - No nugget +
profile likelihood

## Usage

``` r
singl_log_lik_nn(
  theta,
  .dt,
  dists,
  npix,
  model,
  nu = NULL,
  kappa = 1,
  mu2 = 1.5,
  apply_exp = FALSE
)
```

## Arguments

- theta:

  a scalar for the \\\phi\\ parameter.

- .dt:

  a `numeric` vector containing the variable \\Y\\.

- dists:

  a `list` of size three. The first containing the distance matrices
  associated with the regions where \\Y\\ was measured, the second for
  the distance matrices associated with \\X\\, and the last containing
  the cross-distance matrices.

- npix:

  a `integer vector` containing the number of pixels within each
  polygon. (Ordered by the id variables for the polygons).

- model:

  a `character` indicating which covariance function to use. Possible
  values are
  `c("matern", "pexp", "gaussian", "spherical", "cs", "gw", "tapmat")`.

- nu:

  \\\nu\\ parameter. Not necessary if `mode` is `"gaussian"` or
  `"spherical"`

- kappa:

  \\\kappa \in \\0, \ldots, 3 \\\\ parameter for the GW cov function.

- mu2:

  the smoothness parameter \\\mu\\ for the GW function.

- apply_exp:

  a `logical` indicating whether the exponential transformation should
  be applied to variance parameters. This facilitates the optimization
  process.

## Value

a scalar representing `-log.lik`.

## Details

Internal use.
