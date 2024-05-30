
<!-- README.md is generated from README.Rmd. Please edit that file -->

# smile

## **S**patial **M**isalignment: **I**nterpolation, **L**inkage, and **E**stimation

<!-- badges: start -->

[![R-CMD-check](https://github.com/lcgodoy/smile/workflows/R-CMD-check/badge.svg)](https://github.com/lcgodoy/smile/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/smile)](https://CRAN.R-project.org/package=smile)
<!-- badges: end -->

The goal of `smile` is to provide an easy to use `R` package that
functions to estimate, predict and interpolate areal data. For
estimation and prediction we use a model-based approach assuming areal
data to be an average of an underlying continuous Gaussian random field.
For interpolation, we use the simple areal interpolation, which is a
nonparametric approach. Also, we provide functions to quantify the
uncertainty associated with the simple areal interpolation when the
variance of the observed variables is known.

The package accompanies a web page (powered by
[pkgdown](https://pkgdown.r-lib.org/)) and 5 vignettes.

### Vignettes

1.  [Converting `sf` to `spm`
    objects](https://lcgodoy.me/smile/articles/sf-to-spm.html);
2.  [Fitting models and making
    predictions](https://lcgodoy.me/smile/articles/fit-and-pred.html);
3.  [Areal Interpolation](https://lcgodoy.me/smile/articles/sai.html);
4.  [Method](https://lcgodoy.me/smile/articles/theory.html);
5.  [Spatial covariance
    functions](https://lcgodoy.me/smile/articles/sp-cov-functions.html);

### Installation

To install the [CRAN](https://cran.r-project.org) version of the, use

``` r
install.packages("smile")
## remotes::install_github("lcgodoy/smile")
```

The installationg of the development version from Github can be done via

``` r
remotes::install_github("lcgodoy/smile")
## or devtools::install_github("lcgodoy/smile")
```
