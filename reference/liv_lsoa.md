# Liverpool Lower Super Output Area.

A dataset containing containing the LSOA's for Liverpool along with
estimates for Index of Multiple Deprivation. Data taken from [Johnson et
al.
2020](https://ij-healthgeographics.biomedcentral.com/articles/10.1186/s12942-020-00200-w)

## Usage

``` r
liv_lsoa
```

## Format

A `sf` data frame with 298 rows and 6 variables:

- lsoa11cd:

  LSOA code

- lsoa11cd:

  LSOA name

- male:

  Male population

- female:

  Female population

- imdscore:

  Index of Multiple Deprivation

- area:

  LMSOA area, in \\km^2\\

## Source

<https://ij-healthgeographics.biomedcentral.com/articles/10.1186/s12942-020-00200-w>

## Details

The data was projected to EPSG 27700 and units changed to km
