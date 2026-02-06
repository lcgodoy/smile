# Fitting an underlying continuous process to areal data

This function fits a spatial process model to areal data by estimating
the parameters of an underlying continuous process. It leverages the
[stats::optim](https://rdrr.io/r/stats/optim.html) function for Maximum
Likelihood estimation, providing flexibility in optimization algorithms
and control parameters.

## Usage

``` r
fit_spm(x, ...)

# S3 method for class 'spm'
fit_spm(
  x,
  model,
  theta_st,
  nu = NULL,
  kappa = 1,
  mu2 = 1.5,
  apply_exp = FALSE,
  opt_method = "Nelder-Mead",
  control_opt = list(),
  comp_hess = TRUE,
  ...
)

fit_spm2(
  x,
  model,
  nu,
  kappa = 1,
  mu2 = 1.5,
  comp_hess = TRUE,
  phi_min,
  phi_max,
  nphi = 10,
  cores = getOption("mc.cores", 1L)
)
```

## Arguments

- x:

  An object of type `spm`. The dimension of `theta_st` depends on the
  number of variables analyzed and whether the input is an `spm` object.

- ...:

  Additional parameters passed to
  [stats::optim](https://rdrr.io/r/stats/optim.html).

- model:

  A `character` scalar indicating the covariance function family.
  Options are `c("matern", "pexp", "gaussian", "spherical", "gw")`.

- theta_st:

  A named `numeric` vector containing initial parameter values.

- nu:

  A `numeric` value specifying the \\\nu\\ parameter for the Matern or
  Powered Exponential covariance functions. If `model` is "matern" and
  `nu` is not provided, it defaults to 0.5. If `model` is "pexp" and
  `nu` is not provided, it defaults to 1. In both cases, this results in
  the exponential covariance function.

- kappa:

  \\\kappa \in \\0, \ldots, 3 \\\\ parameter for the GW covariance
  function.

- mu2:

  The smoothness parameter \\\mu\\ for the GW function.

- apply_exp:

  A `logical` scalar indicating whether to exponentiate non-negative
  parameters.

- opt_method:

  A `character` scalar specifying the optimization algorithm (see
  [stats::optim](https://rdrr.io/r/stats/optim.html)).

- control_opt:

  A named `list` of control arguments for the optimization algorithm
  (see [stats::optim](https://rdrr.io/r/stats/optim.html)).

- comp_hess:

  A `logical` indicating whether to compute the Hessian matrix.

- phi_min:

  A `numeric` scalar representing the minimum \\phi\\ value for grid
  search.

- phi_max:

  A `numeric` scalar representing the maximum \\phi\\ value for grid
  search.

- nphi:

  A `numeric` scalar indicating the number of \\phi\\ values for grid
  search.

- cores:

  An `integer` scalar indicating the number of cores to use. Defaults to
  `getOption("mc.cores")`. No effect on Windows.

## Value

An `spm_fit` object containing the estimated model parameters.

## Details

The function utilizes optimization algorithms from
[stats::optim](https://rdrr.io/r/stats/optim.html) to determine Maximum
Likelihood estimators (MLEs) and their standard errors for a model
adapted for areal data. Users can customize the optimization process by
providing control parameters through the `control_opt` argument (a named
list, see [stats::optim](https://rdrr.io/r/stats/optim.html) for
details), specifying lower and upper bounds for parameters, and
selecting the desired optimization algorithm via `opt_method` (must be
available in [stats::optim](https://rdrr.io/r/stats/optim.html)).

Additionally, the function supports various covariance functions,
including Matern, Exponential, Powered Exponential, Gaussian, and
Spherical. The `apply_exp` argument, a logical scalar, allows for
exponentiation of non-negative parameters, enabling optimization over
the entire real line.

The underlying model assumes a point-level process: \$\$Y(\mathbf{s}) =
\mu + S(\mathbf{s})\$\$ where \\S ~ GP(0, \sigma^2 C(\lVert \mathbf{s} -
\mathbf{s}\_2 \rVert; \theta))\\. The observed areal data is then
assumed to behave as the average of the process over each area: \$\$Y(B)
= \lvert B \rvert^{-1} \int\_{B} Y(\mathbf{s}) \\ \textrm{d}
\mathbf{s}\$\$.

## Examples

``` r
data(liv_lsoa) ## loading the LSOA data

msoa_spm <- sf_to_spm(sf_obj = liv_msoa, n_pts = 500,
                     type = "regular", by_polygon = FALSE,
                     poly_ids = "msoa11cd",
                     var_ids = "leb_est")
## fitting model
theta_st_msoa <- c("phi" = 1) # initial value for the range parameter

fit_msoa <-
  fit_spm(x = msoa_spm,
          theta_st = theta_st_msoa,
          model = "matern",
          nu = .5,
          apply_exp  = TRUE,
          opt_method = "L-BFGS-B",
          control    = list(maxit = 500))

AIC(fit_msoa)
#> [1] 297.5391

summary_spm_fit(fit_msoa, sig = .05)
#> 
#>  optimization algorithm converged: yes 
#>  
#>     par   estimate        se               ci
#> 1    mu 75.8665481 0.8014059 [74.296; 77.437]
#> 2 sigsq 23.1509888 4.1974039 [14.924; 31.378]
#> 3   phi  0.8211898 0.2152306   [0.399; 1.243]
```
