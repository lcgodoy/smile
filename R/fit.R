##' @name fit_spm
##' @export
fit_spm <- function(x, ...) UseMethod("fit_spm")

##' @name fit_spm
##'
##' @title Fitting an underlying continuous process to areal data
##'
##' @description This function fits a spatial process model to areal data by
##'   estimating the parameters of an underlying continuous process. It
##'   leverages the [stats::optim] function for Maximum Likelihood estimation,
##'   providing flexibility in optimization algorithms and control parameters.
##'
##' @details The function utilizes optimization algorithms from [stats::optim]
##'   to determine Maximum Likelihood estimators (MLEs) and their standard
##'   errors for a model adapted for areal data. Users can customize the
##'   optimization process by providing control parameters through the
##'   `control_opt` argument (a named list, see [stats::optim] for details),
##'   specifying lower and upper bounds for parameters, and selecting the
##'   desired optimization algorithm via `opt_method` (must be available in
##'   [stats::optim]).
##'
##'   Additionally, the function supports various covariance functions,
##'   including Matern, Exponential, and Powered Exponential. The `apply_exp`
##'   argument, a logical scalar, allows for exponentiation of non-negative
##'   parameters, enabling optimization over the entire real line.
##'
##'   The underlying model assumes a point-level process:
##'   \deqn{Y(\mathbf{s}) = \mu + S(\mathbf{s})} where
##'   \eqn{S ~ GP(0, \sigma^2 C(\lVert \mathbf{s} - \mathbf{s}_2 \rVert; \theta))}.
##'   The observed areal data is then assumed to behave as the average of the process
##'   over each area:
##'   \deqn{Y(B) = \lvert B \rvert^{-1} \int_{B} Y(\mathbf{s}) \, \textrm{d} \mathbf{s}}.
##'
##' @param x An object of type \code{spm}. The dimension of \code{theta_st}
##'   depends on the number of variables analyzed and whether the input is an
##'   \code{spm} object.
##' @param model A \code{character} scalar indicating the covariance function
##'   family. Options are \code{c("matern", "pexp")}.
##' @param theta_st A named \code{numeric} vector containing initial parameter
##'   values.
##' @param nu A \code{numeric} value specifying the \eqn{\nu} parameter for the
##'   Matern or Powered Exponential covariance functions. If \code{model} is
##'   "matern" and \code{nu} is not provided, it defaults to 0.5. If
##'   \code{model} is "pexp" and \code{nu} is not provided, it defaults to 1. In
##'   both cases, this results in the exponential covariance function.
##' @param apply_exp A \code{logical} scalar indicating whether to exponentiate
##'   non-negative parameters.
##' @param opt_method A \code{character} scalar specifying the optimization
##'   algorithm (see [stats::optim]).
##' @param control_opt A named \code{list} of control arguments for the
##'   optimization algorithm (see [stats::optim]).
##' @param comp_hess A \code{logical} indicating whether to compute the Hessian
##'   matrix.
##' @param phi_min A \code{numeric} scalar representing the minimum \eqn{phi}
##'   value for grid search.
##' @param phi_max A \code{numeric} scalar representing the maximum \eqn{phi}
##'   value for grid search.
##' @param nphi A \code{numeric} scalar indicating the number of \eqn{phi}
##'   values for grid search.
##' @param cores An \code{integer} scalar indicating the number of cores to use.
##'   Defaults to \code{getOption("mc.cores")}. No effect on Windows.
##' @param ... Additional parameters passed to [stats::optim].
##'
##' @import Matrix
##'
##' @return An \code{spm_fit} object containing the estimated model parameters.
##'
##' @examples
##' data(liv_lsoa) ## loading the LSOA data
##'
##' msoa_spm <- sf_to_spm(sf_obj = liv_msoa, n_pts = 500,
##'                      type = "regular", by_polygon = FALSE,
##'                      poly_ids = "msoa11cd",
##'                      var_ids = "leb_est")
##' ## fitting model
##' theta_st_msoa <- c("phi" = 1) # initial value for the range parameter
##'
##' fit_msoa <-
##'   fit_spm(x = msoa_spm,
##'           theta_st = theta_st_msoa,
##'           model = "matern",
##'           nu = .5,
##'           apply_exp  = TRUE,
##'           opt_method = "L-BFGS-B",
##'           control    = list(maxit = 500))
##'
##' AIC(fit_msoa)
##'
##' summary_spm_fit(fit_msoa, sig = .05)
##'
##' @export
fit_spm.spm <- function(x, model, theta_st,
                        nu = NULL,
                        apply_exp = FALSE,
                        opt_method  = "Nelder-Mead",
                        control_opt = list(),
                        comp_hess = TRUE,
                        ...) {    
  stopifnot(!is.null(names(theta_st)))
  stopifnot(NCOL(x$var) == 1)
  stopifnot(inherits(x, "spm"))
  if (! missing(nu))
    stopifnot(length(nu) == 1)
  stopifnot(model %in% c("matern", "pexp"))
  npar <- length(theta_st)
  p    <- npar + 2L
  if (npar == 2) {
    op_val <-
      stats::optim(par = theta_st,
                   fn  = singl_log_plik,
                   method  = opt_method,
                   control = control_opt,
                   hessian = FALSE,
                   .dt     = x$var,
                   dists   = x$dists,
                   npix    = x$npix,
                   model   = model,
                   nu = nu,
                   apply_exp = apply_exp,
                   ...)
  } else if(npar == 1) {
    op_val <-
      stats::optim(par = theta_st,
                   fn  = singl_log_lik_nn,
                   method  = opt_method,
                   control = control_opt,
                   hessian = FALSE,
                   .dt     = x$var,
                   dists   = x$dists,
                   npix    = x$npix,
                   model   = model,
                   nu = nu,
                   apply_exp = apply_exp,
                   ...)
  }

  estimates <- op_val$par

  if (apply_exp) {
    estimates <- exp(estimates)
    ## Using Delta-Method
    ## https://stats.idre.ucla.edu/r/faq/how-can-i-estimate-the-standard-error-of-transformed-regression-parameters-in-r-using-the-delta-method/
    ## grad_mat <-
    ##     diag(c(1, exp(estimates[2:npar])))
    ## info_mat <- crossprod(grad_mat, info_mat) %*% grad_mat
  }

  .n <- NROW(x$var)

  ## can be turned in to a function to make to code cleaner
  switch(model,
         "matern" = {
           if (is.null(nu))
             nu <- .5
           V <- comp_mat_cov(x$dists,
                             n = .n, n2 = .n,
                             phi   = estimates["phi"],
                             sigsq = 1,
                             nu = nu)
         },
         "pexp" = {
           if (is.null(nu))
             nu <- 1

           V <- comp_pexp_cov(x$dists,
                              n = .n, n2 = .n,
                              phi   = estimates["phi"],
                              sigsq = 1,
                              nu = nu)
         })

  ones_n <- matrix(rep(1, .n), ncol = 1L)
  y <- matrix(x$var, ncol = 1L)

  if (npar == 2) {
    V <- V + diag(estimates["al"] / x$npix,
                  nrow = .n, ncol = .n)
    inv_v <- chol2inv(chol(V))
    mles <- est_mle(x$var, inv_v)
    estimates <- c(mles,
                   "al" = unname(estimates["al"]),
                   "phi" = unname(estimates["phi"]))
    if (comp_hess) {
      info_mat <- solve(
          numDeriv::hessian(func = singl_log_lik,
                            x = estimates,
                            .dt = x$var,
                            dists = x$dists,
                            npix = x$npix,
                            model = model,
                            nu = nu,
                            apply_exp = FALSE)
      )
    } else {
      info_mat <- matrix(NA_real_, ncol = p, nrow = p)
    }
  } else if (npar == 1) {
    inv_v <- chol2inv(chol(V))
    mles  <- est_mle(x$var, inv_v)
    estimates <- c(mles,
                   "phi" = unname(estimates["phi"]))
    if (comp_hess) {
      ## stats::optimHess(par = estimates,
      ##                  fn  = singl_log_lik,
      info_mat <- solve(
          numDeriv::hessian(func = singl_ll_nn_hess,
                            x = estimates,
                            .dt = x$var,
                            dists = x$dists,
                            npix = x$npix,
                            model = model,
                            nu = nu,
                            apply_exp = FALSE)
      )
    } else {
      info_mat <- matrix(NA_real_,
                         ncol = length(estimates),
                         nrow = length(estimates))
    }
  }

  output <- list(
      estimate  = estimates,
      info_mat  = info_mat,
      converged = ifelse(op_val$convergence == 0,
                         "yes", "no"),
      log_lik   = - op_val$value,
      call_data = x,
      model     = model,
      nu        = nu
  )

  class(output) <- append(class(output), "spm_fit")

  return(output)
}

## ##' @name fit_spm
## fit_spm2 <- function(x, model, theta_st,
##                      nu = NULL,
##                      apply_exp = FALSE,
##                      opt_method  = "Nelder-Mead",
##                      control_opt = list(),
##                      comp_hess = TRUE,
##                      ...) {
##     stopifnot(!is.null(names(theta_st)))
##     p    <- NCOL(x$var)
##     npar <- length(theta_st) 
##     op_val <-
##         stats::optim(par = theta_st,
##                      fn  = singl_log_rel,
##                      method  = opt_method,
##                      control = control_opt,
##                      hessian = FALSE,
##                      .dt     = x$var,
##                      dists   = x$dists,
##                      npix    = x$npix,
##                      model   = model,
##                      nu   = nu,
##                      apply_exp = apply_exp,
##                      ...)
##     estimates <- op_val$par
##     if(apply_exp) {
##         ## estimates[2:npar] <- c(expm1(estimates[2]), exp(estimates[3:npar]))
##         estimates <- exp(estimates)
##         ## Using Delta-Method
##         ## https://stats.idre.ucla.edu/r/faq/how-can-i-estimate-the-standard-error-of-transformed-regression-parameters-in-r-using-the-delta-method/
##         ## grad_mat <-
##         ##     diag(c(1, exp(estimates[2:npar])))
##         ## info_mat <- crossprod(grad_mat, info_mat) %*% grad_mat
##     } 
##     .n <- NROW(x$var)
##     ## can be turned in to a function to make to code cleaner
##     switch(model,
##            "matern" = {
##                if(is.null(nu))
##                    nu <- .5
##                V <- comp_mat_cov(x$dists,
##                                  n = .n, n2 = .n,
##                                  phi   = estimates["phi"],
##                                  sigsq = 1,
##                                  nu = nu)
##            },
##            "pexp" = {
##                if(is.null(nu))
##                    nu <- 1
##                V <- comp_pexp_cov(x$dists,
##                                   n = .n, n2 = .n,
##                                   phi   = estimates["phi"],
##                                   sigsq = 1,
##                                   nu = nu)
##            },
##            "gaussian" = {
##                V <- comp_gauss_cov(x$dists,
##                                    n = .n, n2 = .n,
##                                    phi   = estimates["phi"],
##                                    sigsq = 1)
##            },
##            "spherical" = {
##                V <- comp_spher_cov(x$dists,
##                                    n = .n, n2 = .n,
##                                    phi   = estimates["phi"],
##                                    sigsq = 1)
##            })
##     ones_n <- matrix(rep(1, .n), ncol = 1L)
##     y <- matrix(x$var, ncol = 1L)
##     V <- V + diag(estimates["nu"] / x$npix,
##                   nrow = .n, ncol = .n)
##     inv_v <- chol2inv(chol(V))
##     mu_hat <- as.numeric( crossprod(ones_n, inv_v) %*% y / (sum(inv_v)) )
##     estimates <- c("mu" = mu_hat, 
##                    "tausq" = unname(estimates["sigsq"] * estimates["nu"]),
##                    "sigsq" = unname(estimates["sigsq"]),
##                    "phi"   = unname(estimates["phi"]))
##     if(comp_hess) {
##         info_mat <- solve(
##             stats::optimHess(par = estimates,
##                              fn  = singl_log_lik,
##                              .dt = x$var,
##                              dists = x$dists,
##                              npix = x$npix,
##                              model = model,
##                              nu = nu,
##                              apply_exp = FALSE)
##         )
##     } else {
##         info_mat <- matrix(NA_real_, ncol = p, nrow = p)
##     }
##     output <- list(
##         estimate  = estimates,
##         info_mat  = info_mat,
##         converged = ifelse(op_val$convergence == 0,
##                            "yes", "no"),
##         call_data = x,
##         model     = model,
##         nu     = nu
##     )
##     class(output) <- append(class(output), "spm_fit")
##     return(output)
## }

##' @title Summarizing \code{spm_fit}
##'
##' @description Provides a \code{data.frame} with point estimates and
##'     confidence intervals for the parameters of the model fitted using the
##'     \code{spm_fit} function.
##' 
##' @param x a \code{spm_fit} object.
##' @param sig a real number between 0 and 1 indicating significance level to be
##'     used to compute the confidence intervals for the parameter estimates.
##' @return a \code{data.frame} summarizing the parameters estimated by the
##'     \code{fit_spm} function.
##' @export
summary_spm_fit <- function(x, sig = .05) {
  se_est <- diag(x$info_mat)
  est <- x$estimate
  ## cc <- matrix(c(0, 1, 1, 0), ncol = 1)
  ## est <- c(est, as.numeric(crossprod(cc, est)))
  ## se_est <- c(se_est, as.numeric(crossprod(cc, x$info_mat) %*% cc))
  se_est <- sqrt(se_est)
  lb <- est - stats::qnorm(1 - (sig/2)) * se_est
  k  <- 2:length(est)
  lb[k] <- pmax(lb[k], 0)
  ub <- est + stats::qnorm(1 - (sig/2)) * se_est
  ci_est <- sprintf("[%.3f; %.3f]",
                    lb, ub)
  if(is.null(names(x$estimate))) {
    tbl <- data.frame(
        estimate = est,
        se       = se_est,
        ci       = ci_est,
        row.names = NULL
    )
  } else {
    tbl <- data.frame(
        par      = names(x$estimate),
        estimate = est,
        se       = se_est,
        ci       = ci_est,
        row.names = NULL
    )
  }
  cat(sprintf("\n optimization algorithm converged: %s \n \n", x$converged))    
  return(tbl)
}

##' @name fit_spm
##' @export
fit_spm2 <- function(x, model, nu,
                     comp_hess = TRUE,
                     phi_min, phi_max, nphi = 10,
                     cores = getOption("mc.cores", 1L)) {
  stopifnot(NCOL(x$var) == 1)
  stopifnot(inherits(x, "spm"))
  stopifnot(length(nphi) == 1)
  stopifnot(model %in% c("matern", "pexp"))
  if (! missing(nu))
    stopifnot(length(nu) == 1)
  
  my_phi <- seq(from = phi_min,
                to   = phi_max,
                length.out = nphi)

  ## vector to store the profile likelihood value for each phi
  pl <- vector(mode = "numeric",
               length = 2 * length(my_phi))
  
  if(.Platform$OS.type != "unix") {
    for(i in seq_along(my_phi)) {
      pl[i] <- singl_log_lik_nn(my_phi[i], .dt = x$var,
                                dists = x$dists, npix = 1,
                                model = model, nu = nu)
    }
  } else {
    pl <- parallel::mclapply(my_phi,
                             FUN = singl_log_lik_nn,
                             .dt = x$var,
                             dists = x$dists, npix = 1,
                             model = model, nu = nu,
                             mc.cores = cores)
    pl <- unlist(pl)         
  }

  k <- which.min(pl[seq_len(nphi)])

  if(k == 1) {
    my_phi2 <- seq(my_phi[1] - 1e-05, my_phi[1],
                   length.out = length(my_phi))
  } else if(k == length(my_phi)) {
    my_phi2 <- seq(my_phi[nphi], my_phi[nphi] + 1e-05,
                   length.out = length(my_phi))
  } else {
    my_phi2 <- seq(my_phi[k - 1], my_phi[k + 1],
                   length.out = length(my_phi))
  }

  for( i in seq_along(my_phi2) ) {
    pl[i + nphi] <- 
      singl_log_lik_nn(my_phi2[i], .dt = x$var,
                       dists = x$dists, npix = 1,
                       model = model, nu = nu)
  }

  phi_out <- c(my_phi, my_phi2)[which.min(pl)]

  .n <- NROW(x$var)
  
  ## can be turned in to a function to make to code cleaner
  switch(model,
         "matern" = {
           if(is.null(nu))
             nu <- .5

           V <- comp_mat_cov(x$dists,
                             n = .n, n2 = .n,
                             phi   = phi_out,
                             sigsq = 1,
                             nu = nu)
         },
         "pexp" = {
           if(is.null(nu))
             nu <- 1

           V <- comp_pexp_cov(x$dists,
                              n = .n, n2 = .n,
                              phi   = phi_out,
                              sigsq = 1,
                              nu = nu)
         })
  
  ones_n <- matrix(rep(1, .n), ncol = 1L)
  y <- matrix(x$var, ncol = 1L)

  inv_v <- chol2inv(chol(V))
  mles <- est_mle(x$var, inv_v)
  estimates <- c(mles,
                 "phi" = unname(phi_out))

  if(comp_hess) {
    ## stats::optimHess(par = estimates,
    ##                  fn  = singl_log_lik,
    info_mat <- solve(
        numDeriv::hessian(func = singl_ll_nn_hess,
                          x = estimates,
                          .dt = x$var,
                          dists = x$dists,
                          npix = x$npix,
                          model = model,
                          nu = nu,
                          apply_exp = FALSE)
    )
  } else {
    info_mat <- matrix(NA_real_,
                       ncol = length(estimates),
                       nrow = length(estimates))
  }

  output <- list(
      estimate  = estimates,
      info_mat  = info_mat,
      converged = "yes",
      log_lik   = - pl[which.min(pl)],
      call_data = x,
      model     = model,
      nu        = nu
  )
  
  class(output) <- append(class(output), "spm_fit")

  return(output)
}

##' @name goodness_of_fit
##'
##' @title Akaike's (and Bayesian) An Information Criterion for \code{spm_fit}
##'     objects.
##' 
##' @param object a \code{spm_fit} object.
##' @param ... optionally more fitted model objects.
##' @param k \code{numeric}, the _penalty_ per parameter to be used; the default
##'     'k = 2' is the classical AIC. (for compatibility with \code{stats::AIC}.
##'
##' @importFrom stats AIC BIC
##' 
##' @return a \code{numeric} scalar corresponding to the goodness of fit
##'     measure.
##' @export
AIC.spm_fit <- function(object, ..., k = 2) {
  p <- length(object$estimates)
  ll <- object$log_lik
  k * (p - object$log_lik)
}

##' @name goodness_of_fit
##' @export
BIC.spm_fit <- function(object, ...) {
  .n <- length(object$call_data$var)
  ll <- object$log_lik
  log(.n) - 2 * object$log_lik
}
