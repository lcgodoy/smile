##' @name fit_spm
##' @export
fit_spm <- function(x, ...) UseMethod("fit_spm")

##' @name fit_spm
##'
##' @title Fitting an underlying continuous process to areal data
##' 
##' @details This function uses the \code{optim} function optimization
##'     algorithms to find the Maximum Likelihood estimators, and their standard
##'     errors, from a model adapted from. The function allows the user to input
##'     the control parameters from the {optim} function through the argument
##'     \code{control_opt}, which is a named list. Additionally, the one can
##'     input lower and upper boundaries for the optimization problem, as well
##'     as the preferred optimization algorithm (as long as it is available for
##'     \code{optim}). The preferred algorithm is selected by the argument
##'     \code{opt_method}. In addition to the control of the optimization, the
##'     user can select a covariance function among the following: Matern,
##'     Exponential, Powered Exponential, Gaussian, and Spherical. The parameter
##'     \code{apply_exp} is a \code{logical} scalar such that, if set to
##'     \code{TRUE}, the \eqn{\exp} function is applied to the nonnegative
##'     parameters, allowing the optimization algorithm to search for all the
##'     parameters over the real numbers.
##'
##'     The model assumes \deqn{Y(\mathbf{s}) = \mu + S(\mathbf{s})} at the
##'     point level.  Where \eqn{S ~ GP(0, \sigma^2 C(\lVert \mathbf{s} -
##'     \mathbf{s}_2 \rVert; \theta))}.  Further, the observed data is supposed
##'     to be \eqn{Y(B) = \lvert B \rvert^{-1} \int_{B} Y(\mathbf{s}) \,
##'     \textrm{d} \mathbf{s}}.
##' 
##' @param x an object of type \code{spm}. Note that, the dimension of
##'     \code{theta_st} depends on the 2 factors. 1) the number of variables
##'     being analyzed, and 2) if the input is a \code{spm} object.
##' @param model a \code{character} scalar indicating the family of the
##'     covariance function to be used. The options are \code{c("matern",
##'     "pexp", "gaussian", "spherical", "gw")}.
##' @param theta_st a \code{numeric} (named) vector containing the initial
##'     parameters.
##' @param nu a \code{numeric} value indicating either the \eqn{\nu}
##'     paramater from the Matern covariance function (controlling the process
##'     differentiability), or the "pexp" for the Powered Exponential family. If
##'     the \code{model} chosen by the user is Matern and \code{nu} is not
##'     informed, it is automatically set to .5. On the other hand, if the user
##'     choses the Powered Exponential family and do not inform \code{nu},
##'     then it is set to 1. In both cases, the covariance function becomes the
##'     so covalled exponential covariance function.
##' @param tr tapper range
##' @param kappa \eqn{\kappa \in \{0, \ldots, 3 \}} parameter for the GW cov
##'     function.
##' @param mu2 the smoothness parameter \eqn{\mu} for the GW function.
##' @param apply_exp a \code{logical} scalar indicating wheter the parameters
##'     that cannot assume negative values should be exponentiate or not.
##' @param opt_method a \code{character} scalar indicating the optimization
##'     algorithm to be used. For details, see {optim}.
##' @param control_opt a named \code{list} containing the control arguments for
##'     the optimization algorithm to be used. For details, see {optim}.
##' @param comp_hess a \code{boolean} indicating whether the Hessian matrix
##'     should be computed.
##' @param phi_min a \code{numeric} scalar representing the minimum \eqn{phi}
##'     value to look for.
##' @param phi_max a \code{numeric} scalar representing the maximum \eqn{phi}
##'     value to look for.
##' @param nphi a \code{numeric} scalar indicating the number of values to
##'     compute a grid-search over \eqn{phi}.
##' @param ... additionnal parameters, either passed to \code{optim}.
##'
##' @import Matrix
##' 
##' @return a \code{spm_fit} object containing the information about the
##'     estimation of the model parameters.
##'
##' @examples
##' 
##' @examples
##'
##' data(liv_lsoa) ## loading the LSOA data
##' 
##' msoa_spm <- sf_to_spm(sf_obj = liv_msoa, n_pts = 500,
##'                       type = "regular", by_polygon = FALSE,
##'                       poly_ids = "msoa11cd",
##'                       var_ids = "leb_est")
##' ## fitting model
##' theta_st_msoa <- c("phi" = 1) # initial value for the range parameter
##'
##' fit_msoa <-
##'    fit_spm(x = msoa_spm,
##'            theta_st = theta_st_msoa,
##'            model = "matern",
##'            nu = .5,
##'            apply_exp  = TRUE,
##'            opt_method = "L-BFGS-B",
##'            control    = list(maxit = 500))
##'
##' AIC(fit_msoa)
##' 
##' summary_spm_fit(fit_msoa, sig = .05)
##' 
##' @export
fit_spm.spm <- function(x, model, theta_st,
                        nu = NULL,
                        tr = NULL,
                        kappa = 1, mu2 = 1.5,
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
    stopifnot(model %in% c("matern", "pexp", "gaussian",
                           "spherical", "cs", "gw"))
    if (model == "gw")
        stopifnot(mu2 >= 1)
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
                         tr = tr,
                         kappa = kappa,
                         mu2 = mu2,
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
                         tr = tr,
                         kappa = kappa,
                         mu2 = mu2,
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
           },
           "gaussian" = {
               V <- comp_gauss_cov(x$dists,
                                   n = .n, n2 = .n,
                                   phi   = estimates["phi"],
                                   sigsq = 1)
           },
           "spherical" = {
               V <- Matrix(
                   comp_spher_cov(x$dists,
                                  n = .n, n2 = .n,
                                  phi   = estimates["phi"],
                                  sigsq = 1),
                   sparse = TRUE
               )
           },
           "gw" = {
               V <- Matrix(
                   comp_gw_cov(x$dists,
                               n = .n, n2 = .n,
                               phi   = estimates["phi"],
                               sigsq = 1,
                               kappa = kappa,
                               mu    = mu2),
                   sparse = TRUE
               )
           },
           "cs" = {
               V <- Matrix(
                   comp_cs_cov(x$dists,
                               n = .n, n2 = .n,
                               phi   = estimates["phi"],
                               sigsq = 1),
                   sparse = TRUE
               )
           },
           "tapmat" = {
               V <- Matrix(
                   comp_tapmat_cov(x$dists,
                                   n = .n, n2 = .n,
                                   phi   = estimates["phi"],
                                   sigsq = 1,
                                   nu = nu,
                                   theta = tr),
                   sparse = TRUE
               )
           })

    ones_n <- matrix(rep(1, .n), ncol = 1L)
    y <- matrix(x$var, ncol = 1L)

    if (npar == 2) {
        V <- V + diag(estimates["al"] / x$npix,
                      nrow = .n, ncol = .n)
        inv_v <- chol2inv(chol(V))
        mles <- est_mle(x$var, inv_v)
        estimates <- c(mles,
                       "al" = estimates["al"],
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
                                  tr = tr,
                                  kappa = kappa,
                                  mu2 = mu2,
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
                                  tr = tr,
                                  kappa = kappa,
                                  mu2 = mu2,
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
        nu        = nu,
        taper_rg  = tr,
        gw_pars   = c(kappa, mu2)
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
##'     confidence intervals for the paramters of the model fitted using the
##'     \code{spm_fit} function.
##' 
##' @param x a \code{spm_fit} object.
##' @param sig a real number between 0 and 1 indicating significance level to be
##'     used to compute the confidence intervals for the parameter estimates.
##' @return a \code{data.frame} summarising the parameters estimated by the
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
                     tr,
                     kappa = 1, mu2 = 1.5,
                     comp_hess = TRUE, 
                     phi_min, phi_max, nphi = 10) {
    stopifnot(NCOL(x$var) == 1)
    stopifnot(inherits(x, "spm"))
    stopifnot(length(nphi) == 1)
    stopifnot(model %in% c("matern", "pexp", "gaussian",
                           "spherical", "cs", "gw"))
    if(model == "gw")
        stopifnot(mu2 >= 1)
    if(! missing(nu))
        stopifnot(length(nu) == 1)
    
    my_phi <- seq(from = phi_min,
                  to   = phi_max,
                  length.out = nphi)

    ## vector to store the profile likelihood value for each phi
    pl <- vector(mode = "numeric",
                 length = 2 * length(my_phi))

    for( i in seq_along(my_phi) ) {
        pl[i] <- singl_log_lik_nn(my_phi[i], .dt = x$var,
                                  dists = x$dists, npix = 1,
                                  model = model, nu = nu,
                                  kappa = kappa, mu2 = mu2)
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
                             model = model, nu = nu,
                             kappa = kappa, mu2 = mu2)
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
           },
           "gaussian" = {
               V <- comp_gauss_cov(x$dists,
                                   n = .n, n2 = .n,
                                   phi   = phi_out,
                                   sigsq = 1)
           },
           "spherical" = {
               V <- Matrix(
                   comp_spher_cov(x$dists,
                                  n = .n, n2 = .n,
                                  phi   = phi_out,
                                  sigsq = 1),
                   sparse = TRUE
               )
           },
           "gw" = {
               V <- Matrix(
                   comp_gw_cov(x$dists,
                               n = .n, n2 = .n,
                               phi   = phi_out,
                               sigsq = 1,
                               kappa = kappa,
                               mu    = mu2),
                   sparse = TRUE
               )
           },
           "cs" = {
               V <- Matrix(
                   comp_cs_cov(x$dists,
                               n = .n, n2 = .n,
                               phi   = phi_out,
                               sigsq = 1),
                   sparse = TRUE
               )
           },
           "tapmat" = {
               V <- Matrix(
                   comp_tapmat_cov(x$dists,
                                   n = .n, n2 = .n,
                                   phi   = phi_out,
                                   sigsq = 1,
                                   nu = nu,
                                   theta = tr),
                   sparse = TRUE
               )
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
                              tr    = tr,
                              kappa = kappa,
                              mu2   = mu2,
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
        nu        = nu,
        gw_pars   = c(kappa, mu2)
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
