# single `sf` to `spm`

Transforming a `sf` into a `spm` object (Internal use)

## Usage

``` r
single_sf_to_spm(
  sf_obj,
  n_pts,
  type = "regular",
  by_polygon = FALSE,
  poly_ids = NULL,
  var_ids = NULL
)

sf_to_spm(
  sf_obj,
  n_pts,
  type = "regular",
  by_polygon = FALSE,
  poly_ids = NULL,
  var_ids = NULL
)
```

## Arguments

- sf_obj:

  a `sf` object s.t. its geometries are polygons.

- n_pts:

  a `numeric` scalar representing the number of points to create a grid
  in the study region on which the polygons in `sf_obj` is observed.
  Alternatively, it can be a vector of the same length as
  `nrow(sf_obj)`. In this case, it generates the given number of points
  for each polygon in `sf_obj`.

- type:

  a `character` indicating the type of grid to be generated. The options
  are `c("random", "regular", "hexagonal")`. For more details, see
  `st_sample` in the `sf` package.

- by_polygon:

  a `logical` indicating whether we should generate `n_pts` by polygon
  or for the `n_pts` for the whole study region.

- poly_ids:

  a `character` vector informing the name of the variable in `sf_obj`
  that represents the polygons unique identifiers. In case this is not
  informed, we assume the id of the polygons are given by their row
  numbers.

- var_ids:

  a scalar or vector of type `character` indicating the (numerical)
  variables that are going to be analyzed.

## Value

a named `list` of size 6 belonging to the class `spm`. This list stores
all the objects necessary to fit models using the
[`fit_spm`](https://lcgodoy.me/smile/reference/fit_spm.md).

## Examples

``` r
data(liv_lsoa) # loading the LSOA data

msoa_spm <- sf_to_spm(sf_obj = liv_msoa, n_pts = 1000,
                      type = "regular", by_polygon = FALSE,
                      poly_ids = "msoa11cd",
                      var_ids = "leb_est")
```
