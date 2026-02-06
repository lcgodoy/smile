# Summarizing `spm_fit`

Provides a `data.frame` with point estimates and confidence intervals
for the parameters of the model fitted using the `spm_fit` function.

## Usage

``` r
summary_spm_fit(x, sig = 0.05)
```

## Arguments

- x:

  a `spm_fit` object.

- sig:

  a real number between 0 and 1 indicating significance level to be used
  to compute the confidence intervals for the parameter estimates.

## Value

a `data.frame` summarizing the parameters estimated by the `fit_spm`
function.
