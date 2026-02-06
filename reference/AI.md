# Areal Interpolation

This function estimates variables observed at a "source" region into a
"target" region. "Source" and "target" regions represent two different
ways to divide a city, for example. For more details, see
<https://lcgodoy.me/smile/articles/sai.html>.

## Usage

``` r
ai(source, target, vars)

ai_var(source, target, vars, vars_var, sc_vars = FALSE, var_method = "CS")
```

## Arguments

- source:

  a `sf` object - source spatial data.

- target:

  a `sf` object - target spatial data.

- vars:

  a `character` representing the variables (observed at the source) to
  be estimated at the target data.

- vars_var:

  a scalar of type `character` representing the name of the variable in
  the source dataset that stores the variances of the variable to be
  estimated at the target data.

- sc_vars:

  boolean indicating whether `vars` should be scaled by its observed
  variance (if available).

- var_method:

  a `character` representing the method to approximate the variance of
  the AI estimates. Possible values are "CS" (Cauchy-Schwartz) or "MI"
  (Moran's I).

## Value

the target (of type `sf`) with estimates of the variables observed at
the source data.
