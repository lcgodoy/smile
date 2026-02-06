# Liverpool Middle Super Output Area.

A dataset containing containing the MSOA's for Liverpool along with
estimates for Life Expectancy at Birth. Data taken from [Johnson et al.
2020](https://ij-healthgeographics.biomedcentral.com/articles/10.1186/s12942-020-00200-w)

## Usage

``` r
liv_msoa
```

## Format

A `sf` data frame with 61 rows and 4 variables:

- msoa11cd:

  MSOA code

- msoa11cd:

  MSOA name

- lev_est:

  Estimated life expectancy at birth, in years

- area:

  MSOA area, in \\km^2\\

## Source

<https://ij-healthgeographics.biomedcentral.com/articles/10.1186/s12942-020-00200-w>

## Details

The data was projected to EPSG 27700 and units changed to km
