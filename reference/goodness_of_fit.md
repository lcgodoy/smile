# Akaike's (and Bayesian) An Information Criterion for `spm_fit` objects.

Akaike's (and Bayesian) An Information Criterion for `spm_fit` objects.

## Usage

``` r
# S3 method for class 'spm_fit'
AIC(object, ..., k = 2)

# S3 method for class 'spm_fit'
BIC(object, ...)
```

## Arguments

- object:

  a `spm_fit` object.

- ...:

  optionally more fitted model objects.

- k:

  `numeric`, the *penalty* per parameter to be used; the default 'k = 2'
  is the classical AIC. (for compatibility with
  [`stats::AIC`](https://rdrr.io/r/stats/AIC.html).

## Value

a `numeric` scalar corresponding to the goodness of fit measure.
