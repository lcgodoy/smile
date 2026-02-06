# Building weight matrix **W** for Areal Interpolation

internal use. \\W\_{ij} = \| A_i \\ cap \\ B_j \|\\.

## Usage

``` r
w_col(source_unit, target)

build_w(source, target)

est_w(W, source_dt, target)

var_w(W, var_vec, target, method = "CS", rho_mi)
```

## Arguments

- source_unit:

  a single `geometry` from the source dataset.

- target:

  a `sf` object - target spatial data.

- source:

  a `sf` object - source spatial data.

- W:

  the weight matrix.

- source_dt:

  a `data.frame` object representing the source dataset but excluding
  the `geometry`, i.e. the spatial information, column.

- var_vec:

  a `numeric` vector with variances observed at the source data.

- method:

  a `character` representing the method to approximate the variance of
  the AI estimates. Possible values are "CS" (Cauchy-Schwartz) or "MI"
  (Moran's I).

- rho_mi:

  `numeric` calculated Moran's I.

## Value

A \\n \times m\\ `numeric` matrix. Where \\n\\ is the number of
observations in the target and \\m\\ is the sample size in the source
dataset.
