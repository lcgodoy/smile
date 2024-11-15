---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# smile

## **S**patial **M**isalignment: **I**nterpolation, **L**inkage, and **E**stimation

<!-- badges: start -->
[![R-CMD-check](https://github.com/lcgodoy/smile/workflows/R-CMD-check/badge.svg)](https://github.com/lcgodoy/smile/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/smile)](https://CRAN.R-project.org/package=smile)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/smile)](https://CRAN.R-project.org/package=smile)
<!-- badges: end -->

The `smile` package for R simplifies the analysis of spatial data. It offers
tools to:

* **Estimate and predicte**: Using a model-based approach, smile treats areal
  data as averages derived from a continuous Gaussian random field. It is
  possible to model point-referenced and areal data jointly.

* **Interpolate**: For straightforward interpolation, smile employs a
  non-parametric method known as simple areal interpolation.

* **Quantify uncertainty:** When dealing with simple areal interpolation and
  knowing the variance of your observations, smile helps you quantify the
  uncertainty associated with the interpolation.

Essentially, smile provides a user-friendly way to explore, model, and
interpolate spatial data with R, offering both model-based and non-parametric
approaches.

The package accompanies a web page (powered by
[pkgdown](https://pkgdown.r-lib.org/)) and 5 vignettes.

### Vignettes

1. [Converting `sf` to `spm`
   objects](https://lcgodoy.me/smile/articles/sf-to-spm.html);
2. [Fitting models and making
   predictions](https://lcgodoy.me/smile/articles/fit-and-pred.html);
3. [Areal Interpolation](https://lcgodoy.me/smile/articles/sai.html);
4. [Method](https://lcgodoy.me/smile/articles/theory.html);
5. [Spatial covariance functions](https://lcgodoy.me/smile/articles/sp-cov-functions.html);

### Installation

To install the [CRAN](https://cran.r-project.org) version of the package, use
```r
install.packages("smile")
```

The installation of the development version from GitHub can be done via
```r
remotes::install_github("lcgodoy/smile")
## or devtools::install_github("lcgodoy/smile")
```
