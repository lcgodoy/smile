# Evaluate log-lik

Evaluate the log-likelihood for a given set of parameters

## Usage

``` r
singl_log_lik(
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

  a `numeric` vector of size 4 (\\\mu, \sigma^2, \alpha, \phi\\)
  containing the parameters associated with the model.

- .dt:

  a `numeric` vector containing the variable \\Y\\.

- dists:

  a `list` of size distance matrices at the point level.

- npix:

  a `integer vector` containing the number of pixels within each
  polygon. (Ordered by the id variables for the polygons).

- model:

  a `character` indicating which covariance function to use. Possible
  values are `c("matern", "pexp", "gaussian", "spherical", "cs", "gw")`.

- nu:

  \\\nu\\ parameter. Not necessary if `model` is `"gaussian"` or
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
