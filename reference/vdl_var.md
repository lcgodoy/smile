# Voronoi Data Linkage - Single variable and variance

Reminder, have to create an example. This will be exported after we
submit the paper for publication.

## Usage

``` r
vdl_var(coords_sf, areal_sf, res_var, variance, var_method = "CS", buff)
```

## Arguments

- coords_sf:

  `sf` POINT target dataset.

- areal_sf:

  `sf` POLYGON source dataset.

- res_var:

  a `character` - the name of the variable in the `areal_sf` to be
  estimated in the `coords_sf`.

- variance:

  a `character` - the name of the variable variance in the `areal_sf` to
  be estimated in the `coords_sf`.

- var_method:

  a `character` representing the method to approximate the variance of
  the AI estimates. Possible values are "CS" (Cauchy-Schwartz) or "MI"
  (Moran's I).

- buff:

  scalar `numeric`. Mostly for internal use.

## Value

a `sf` object, containing the `id_coords` variable and the `list_vars`
for the `coords_sf` spatial data set.
