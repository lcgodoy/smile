# Voronoi Data Linkage

Reminder, have to create an example. This will be exported after we
submit the paper for publication.

## Usage

``` r
vdl(coords_sf, areal_sf, vars, buff)
```

## Arguments

- coords_sf:

  `sf` POINT target dataset.

- areal_sf:

  `sf` POLYGON source dataset.

- vars:

  a `character` representing the variables (observed at the source -
  polygon) to be estimated at the target data.

- buff:

  scalar `numeric`. Mostly for internal use.

## Value

a `sf` object for the `coords_sf` spatial data set.
