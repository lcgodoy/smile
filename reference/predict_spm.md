# Prediction over the same or a different set of regions (or points).

Realizes predictions that can be useful when researchers are interested
in predict a variable observed in one political division of a city (or
state) on another division of the same region.

## Usage

``` r
predict_spm(x, ...)

# S3 method for class 'spm_fit'
predict_spm(x, .aggregate = TRUE, ...)

# S3 method for class 'sf'
predict_spm(x, spm_obj, n_pts, type, outer_poly = NULL, id_var, ...)
```

## Arguments

- x:

  a `sf` object such that its geometries are either points or polygons.

- ...:

  additional parameters

- .aggregate:

  `logical`. Should the predictions be aggregated? In case the input is
  only a "fit" object, the aggregation is made over the polygons on
  which the original data was observed. In case the input `x` is
  composed by `sf POLYGONS`, the aggregation is made over this new
  partition of the study region.

- spm_obj:

  an object of either class `spm_fit` or `mspm_fit`

- n_pts:

  a `numeric` scalar standing for number of points to form a grid over
  the whole region to make the predictions

- type:

  `character` type of grid to be generated. See `st_sample` in the
  package `sf`.

- outer_poly:

  (object) `sf geometry` storing the "outer map" we want to compute the
  predictions in.

- id_var:

  if `x` is a set of `POLYGONS` (areal data) instead of a set of points,
  the `id_var` is the name (or index) of the unique identifier
  associated to each polygon.

## Value

a `list` of size 4 belonging to the class `spm_pred`. This list contains
the predicted values and the mean and covariance matrix associated with
the conditional distribution used to compute the predictions.

## Examples

``` r
data(liv_lsoa) ## loading the LSOA data
data(liv_msoa) ## loading the MSOA data

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

pred_lsoa <- predict_spm(x = liv_lsoa, spm_obj = fit_msoa, id_var = "lsoa11cd")
#> Warning: st_centroid assumes attributes are constant over geometries
```
