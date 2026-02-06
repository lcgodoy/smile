# Voronoi Tessellation inside a polygon

voronoi tessellation of a given a set of points inside a polygon. This
is an internal use function.

## Usage

``` r
vor_build(points_sf, poly_sf)
```

## Arguments

- points_sf:

  `sf data frame` containing the points' coordinates

- poly_sf:

  a `sf` polygon

## Value

a `sf` object containing the polygons associated with the voronoi
tessellation of `points_sf` the polygon `poly_sf`
