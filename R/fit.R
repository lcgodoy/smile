##' @name fit_spm
##' @export
fit_spm <- function(x, ...) UseMethod("fit_spm")

##' @name fit_spm
##'
##' @title  Fit a model adapted from
##' \href{https://ij-healthgeographics.biomedcentral.com/articles/10.1186/s12942-020-00200-w}{Johnson
##' et al. 2020}.
##'
##' @details This function uses the \code{\link[stats]{optim}} optimization
##'     algorithms to find the Maximum Likelihood estimators, and their standard
##'     errors, from a model adapted from
##'     \href{https://ij-healthgeographics.biomedcentral.com/articles/10.1186/s12942-020-00200-w}{Johnson
##'     et al. 2020}. The function allows the user to input the control
##'     parameters from the \code{\link[stats]{optim}} function through the
##'     argument \code{control_opt}, which is a named list. Additionally, the
##'     user can input lower and upper boundaries for the optimization problem,
##'     as well as the preferred optimization algorithm (as long as it is
##'     available at \code{\link[stats]{optim}}). The preferred algorithm is
##'     selected by the argument \code{opt_method}. In addition to the control
##'     of the optimization, the user can select a covariance function among the
##'     following: Matern, Exponential, Powered Exponential, Gaussian, and
##'     Spherical. The parameter \code{apply_exp} is a \code{logical} scalar
##'     such that, if set to \code{TRUE}, the \eqn{\exp} function is applied to
##'     the nonnegative parameters, allowing the optimization algorithm to
##'     search for all the paraters over the real numbers.
##' 
##' @title Fitting the model based on
##' @param x an object of type \code{sspm} or \code{mspm}. Note that, the
##'     dimension of \code{theta_st} depends on the 2 factors. 1) the number of
##'     variables being analyzed, and 2) if the input is a \code{sspm} or a
##'     \code{mspm} object.
##' @param model a \code{character} scalar indicating the family of the
##'     covariance function to be used. The options are \code{c("matern",
##'     "pexp", "gaussian", "spherical")}.
##' @param kappa a \code{numeric} value indicating either the \eqn{\kappa}
##'     paramater from the Matern covariance function (controlling the process
##'     differentiability), or the "pexp" for the Powered Exponential family. If
##'     the \code{model} chosen by the user is Matern and \code{kappa} is not
##'     informed, it is automatically set to .5. On the other hand, if the user
##'     choses the Powered Exponential family and do not inform \code{kappa},
##'     then it is set to 1. In both cases, the covariance function becomes the
##'     so covalled exponential covariance function.
##' @param theta_st a \code{numeric} (named) vector containing the initial
##'     parameters.
##' @param apply_exp a \code{logical} scalar indicating wheter the parameters
##'     that cannot assume negative values should be exponentiate or not.
##' @param opt_method a \code{character} scalar indicating the optimization
##'     algorithm to be used. For details, see \code{\link[stats]{optim}}.
##' @param control_opt a named \code{list} containing the control arguments for
##'     the optimization algorithm to be used. For details, see
##'     \code{\link[stats]{optim}}.
##' @param comp_hess a \code{boolean} indicating whether the Hessian matrix
##'     should be computed
##' @param ... additionnal parameters, either passed to \code{optim}.
##' 
##' @return a \code{sspm_fit} or a \code{mspm_fit} object.
##' 
##' @export
fit_spm.sspm <- function(x, model, theta_st,
                         kappa = NULL,
                         apply_exp = FALSE,
                         opt_method  = "Nelder-Mead",
                         control_opt = list(),
                         comp_hess   = TRUE,
                         ...) {
    
    stopifnot(!is.null(names(theta_st)))

    p <- NCOL(x$var)
    npar <- length(theta_st) 
    
    stopifnot(npar == (p + 2 + (.5 * p * (p + 1))))
    
    op_val <-
        stats::optim(par = theta_st,
                     fn  = singl_log_lik,
                     method  = opt_method,
                     control = control_opt,
                     hessian = comp_hess,
                     .dt     = x$var,
                     dists   = x$dists,
                     npix    = x$npix,
                     model   = model,
                     kappa   = kappa,
                     apply_exp = apply_exp,
                     ...)

    if(comp_hess) 
        info_mat <- solve(op_val$hessian)
    else
        info_mat <- matrix(NA, ncol = p, nrow = p)

    estimates <- op_val$par
    
    if(apply_exp) {
        ## estimates[2:npar] <- c(expm1(estimates[2]), exp(estimates[3:npar]))
        estimates[2:npar] <- exp(estimates[2:npar])
        ## Using Delta-Method
        ## https://stats.idre.ucla.edu/r/faq/how-can-i-estimate-the-standard-error-of-transformed-regression-parameters-in-r-using-the-delta-method/
        if(comp_hess) {
            grad_mat <-
                diag(c(1, exp(estimates[2:npar])))
            info_mat <- crossprod(grad_mat, info_mat) %*% grad_mat
        }
    } 

    output <- list(
        estimate  = estimates,
        info_mat  = info_mat,
        converged = ifelse(op_val$convergence == 0,
                           "yes", "no"),
        call_data = x,
        model     = model,
        kappa     = kappa
    )

    class(output) <- append(class(output), "sspm_fit")

    return(output)
}

##' @name fit_spm
##' @export
fit_spm2 <- function(x, model, theta_st,
                     kappa = NULL,
                     apply_exp = FALSE,
                     opt_method  = "Nelder-Mead",
                     control_opt = list(),
                     comp_hess = TRUE,
                     ...) {
    
    stopifnot(!is.null(names(theta_st)))

    p    <- NCOL(x$var)
    npar <- length(theta_st) 
    
    ## stopifnot(npar == (p + 2 + (.5 * p * (p + 1))))
    
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
                     kappa   = kappa,
                     apply_exp = apply_exp,
                     ...)
    
    ## info_mat  <- solve(op_val$hessian)
    estimates <- op_val$par
    
    if(apply_exp) {
        ## estimates[2:npar] <- c(expm1(estimates[2]), exp(estimates[3:npar]))
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
               if(is.null(kappa))
                   kappa <- .5

               V <- comp_mat_cov(x$dists,
                                 n = .n, n2 = .n,
                                 phi   = estimates["phi"],
                                 sigsq = 1,
                                 kappa = kappa)
           },
           "pexp" = {
               if(is.null(kappa))
                   kappa <- 1

               V <- comp_pexp_cov(x$dists,
                                  n = .n, n2 = .n,
                                  phi   = estimates["phi"],
                                  sigsq = 1,
                                  kappa = kappa)
           },
           "gaussian" = {
               V <- comp_gauss_cov(x$dists,
                                   n = .n, n2 = .n,
                                   phi   = estimates["phi"],
                                   sigsq = 1)
           },
           "spherical" = {
               V <- comp_spher_cov(x$dists,
                                   n = .n, n2 = .n,
                                   phi   = estimates["phi"],
                                   sigsq = 1)
           })
    
    ones_n <- matrix(rep(1, .n), ncol = 1L)
    y <- matrix(x$var, ncol = 1L)

    V <- V + diag(estimates["nu"] / x$npix,
                  nrow = .n, ncol = .n)

    inv_v <- chol2inv(chol(V))
    
    mu_hat    <- as.numeric( crossprod(ones_n, inv_v) %*% y / (sum(inv_v)) )
    y_mu_hat  <- y - mu_hat
    sigsq_hat <- as.numeric( crossprod(y_mu_hat, inv_v) %*% y_mu_hat / .n )

    estimates <- c("mu" = mu_hat, 
                   "tausq" = unname(sigsq_hat * estimates["nu"]),
                   "sigsq" = sigsq_hat,
                   "phi"   = unname(estimates["phi"]))

    if(comp_hess) {
        info_mat <- solve(
            stats::optimHess(par = estimates,
                             fn  = singl_log_lik,
                             .dt = x$var,
                             dists = x$dists,
                             npix = x$npix,
                             model = model,
                             kappa = kappa,
                             apply_exp = FALSE)
        )
    } else {
        info_mat <- matrix(NA_real_, ncol = p, nrow = p)
    }
    
    output <- list(
        estimate  = estimates,
        info_mat  = info_mat,
        converged = ifelse(op_val$convergence == 0,
                           "yes", "no"),
        call_data = x,
        model     = model,
        kappa     = kappa
    )

    class(output) <- append(class(output), "sspm_fit")

    return(output)
}

##' @name fit_spm
##' @export
fit_spm3 <- function(x, model, theta_st,
                     kappa = NULL,
                     apply_exp = FALSE,
                     opt_method  = "Nelder-Mead",
                     control_opt = list(),
                     comp_hess = TRUE,
                     ...) {
    
    stopifnot(!is.null(names(theta_st)))

    p    <- NCOL(x$var)
    npar <- length(theta_st) 
    
    ## stopifnot(npar == (p + 2 + (.5 * p * (p + 1))))
    
    op_val <-
        stats::optim(par = theta_st,
                     fn  = singl_log_rel,
                     method  = opt_method,
                     control = control_opt,
                     hessian = FALSE,
                     .dt     = x$var,
                     dists   = x$dists,
                     npix    = x$npix,
                     model   = model,
                     kappa   = kappa,
                     apply_exp = apply_exp,
                     ...)
    
    ## info_mat  <- solve(op_val$hessian)
    estimates <- op_val$par
    
    if(apply_exp) {
        ## estimates[2:npar] <- c(expm1(estimates[2]), exp(estimates[3:npar]))
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
               if(is.null(kappa))
                   kappa <- .5

               V <- comp_mat_cov(x$dists,
                                 n = .n, n2 = .n,
                                 phi   = estimates["phi"],
                                 sigsq = 1,
                                 kappa = kappa)
           },
           "pexp" = {
               if(is.null(kappa))
                   kappa <- 1

               V <- comp_pexp_cov(x$dists,
                                  n = .n, n2 = .n,
                                  phi   = estimates["phi"],
                                  sigsq = 1,
                                  kappa = kappa)
           },
           "gaussian" = {
               V <- comp_gauss_cov(x$dists,
                                   n = .n, n2 = .n,
                                   phi   = estimates["phi"],
                                   sigsq = 1)
           },
           "spherical" = {
               V <- comp_spher_cov(x$dists,
                                   n = .n, n2 = .n,
                                   phi   = estimates["phi"],
                                   sigsq = 1)
           })
    
    ones_n <- matrix(rep(1, .n), ncol = 1L)
    y <- matrix(x$var, ncol = 1L)

    V <- V + diag(estimates["nu"] / x$npix,
                  nrow = .n, ncol = .n)

    inv_v <- chol2inv(chol(V))
    
    mu_hat <- as.numeric( crossprod(ones_n, inv_v) %*% y / (sum(inv_v)) )

    estimates <- c("mu" = mu_hat, 
                   "tausq" = unname(estimates["sigsq"] * estimates["nu"]),
                   "sigsq" = unname(estimates["sigsq"]),
                   "phi"   = unname(estimates["phi"]))

    if(comp_hess) {
        info_mat <- solve(
            stats::optimHess(par = estimates,
                             fn  = singl_log_lik,
                             .dt = x$var,
                             dists = x$dists,
                             npix = x$npix,
                             model = model,
                             kappa = kappa,
                             apply_exp = FALSE)
        )
    } else {
        info_mat <- matrix(NA_real_, ncol = p, nrow = p)
    }
    
    output <- list(
        estimate  = estimates,
        info_mat  = info_mat,
        converged = ifelse(op_val$convergence == 0,
                           "yes", "no"),
        call_data = x,
        model     = model,
        kappa     = kappa
    )

    class(output) <- append(class(output), "sspm_fit")

    return(output)
}

##' @name fit_spm
##' @export
fit_spm.mspm <- function(x, model, theta_st,
                         kappa = NULL,
                         apply_exp = FALSE,
                         opt_method  = "Nelder-Mead",
                         control_opt = list(),
                         ...) {
    
    stopifnot(!is.null(names(theta_st)))

    p    <- NCOL(x$y$var)
    .q   <- NCOL(x$x$var)
    npar <- length(theta_st) 

    stopifnot(.q == 1)
    stopifnot(npar == (2 * p + 4 + (.5 * p * (p + 1))))
    
    op_val <-
        stats::optim(par = theta_st,
                     fn  = mult_log_lik,
                     method  = opt_method,
                     control = control_opt,
                     hessian = TRUE,
                     .dt     = list(x$y$var, x$x$var),
                     dists   = list(x$y$dists, x$x$dists, x$cdist),
                     model   = model,
                     kappa   = kappa,
                     apply_exp = apply_exp,
                     ...)
    
    info_mat  <- solve(op_val$hessian)
    estimates <- op_val$par
    
    if(apply_exp) {
        estimates[((2 * p) + 2):npar] <- exp(estimates[((2 * p) + 2):npar])
        ## Using Delta-Method
        ## https://stats.idre.ucla.edu/r/faq/how-can-i-estimate-the-standard-error-of-transformed-regression-parameters-in-r-using-the-delta-method/
        grad_mat <-
            diag(c(rep(1, (2 * p) + 1), estimates[((2 * p) + 2):npar]))
        info_mat <- crossprod(grad_mat, info_mat) %*% info_mat
    } 

    output <- list(
        estimate  = estimates,
        info_mat  = info_mat,
        converged = ifelse(op_val$convergence == 0,
                           "yes", "no"),
        call_data = x,
        model     = model,
        kappa     = kappa
    )

    class(output) <- append(class(output), "mspm_fit")

    return(output)
}

##' @title summarizing \code{spm_fit}
##'
##' @param x either a \code{sspm_fit} or a \code{mspm_fit} object
##' @param sig signigicance level
##' @export
summary_sspm_fit <- function(x, sig = .05) {
    se_est <- sqrt(diag(x$info_mat))
    lb <- x$estimate - stats::qnorm(1 - (sig/2)) * se_est
    k  <- which(grepl("(tausq|phi|sigsq)", names(x$estimate)))
    lb[k] <- pmax(lb[k], 0)
    ub <- x$estimate + stats::qnorm(1 - (sig/2)) * se_est
    ci_est <- sprintf("[%.3f; %.3f]",
                      lb, ub)
    if(is.null(names(x$estimate))) {
        tbl <- data.frame(
            estimate = x$estimate,
            se       = se_est,
            ci       = ci_est,
            row.names = NULL
        )
    } else {
        tbl <- data.frame(
            par      = names(x$estimate),
            estimate = x$estimate,
            se       = se_est,
            ci       = ci_est,
            row.names = NULL
        )
    }

    cat(sprintf("\n optimization algorithm converged: %s \n \n", x$converged))
    
    return(tbl)
}

## summary.mspm_fit <- function(x, sig = .05, ...) {
##     se_est <- sqrt(diag(x$info_mat))
##     ci_est <- sprintf("[%.3f; %.3f]",
##                       x$estimate - stats::qnorm(1 - (sig/2)) * se_est,
##                       x$estimate + stats::qnorm(1 - (sig/2)) * se_est)
##     if(is.null(names(x$estimate))) {
##         tbl <- data.frame(
##             estimate = x$estimate,
##             se       = se_est,
##             ci       = ci_est,
##             row.names = NULL
##         )
##     } else {
##         tbl <- data.frame(
##             par      = names(x$estimate),
##             estimate = x$estimate,
##             se       = se_est,
##             ci       = ci_est,
##             row.names = NULL
##         )
##     }

##     cat(sprintf("\n optimization algorithm converged: %s", x$converged))

##     return(tbl)
## }


##' @export
print.mspm_fit <- function(x, ...) {
    summary_sspm_fit(x, ...)
}

