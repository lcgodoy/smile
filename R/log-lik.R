##' Evaluate the log-likelihood for a given set of parameters
##'
##' Internal use.
##' @title Evaluate log-lik
##' @param theta a \code{numeric} vector of size 4 (\eqn{\mu, \sigma^2, \tau^2,
##'     \phi}) containing the parameters associated with the model.
##' @param .dt a \code{numeric} vector containing the variable \eqn{Y}.
##' @param dists a \code{list} of size distance matrices at the point level.
##' @param npix a \code{integer vector} containing the number of pixels within
##'     each polygon. (Ordered by the id variables for the polygons).
##' @param model a \code{character} indicating which covariance function to
##'     use. Possible values are \code{c("matern", "pexp", "gaussian",
##'     "spherical", "cs", "w1", "tapmat")}.
##' @param kappa \eqn{\kappa} parameter. Not necessary if \code{model} is
##'     \code{"gaussian"} or \code{"spherical"}
##' @param tr \eqn{\theta_r} taper range.
##' @param apply_exp a \code{logical} indicater wheter the exponential
##'     transformation should be applied to variance parameters. This
##'     facilitates the optimization process.
##' 
##' @return a scalar representing \code{-log.lik}.
##' @keywords internal
singl_log_lik <- function(theta, .dt, dists, npix, model,
                          kappa = NULL, tr = NULL,
                          apply_exp = FALSE) {

    if(! apply_exp & any(theta[3:4] < 0 )) {
        return(NA_real_)
    }

    mu    <- theta[1]
    sigsq <- theta[2]
    tausq <- theta[3]
    phi   <- theta[4]

    ## tausq <- matrix(nrow = p, ncol = p)
    ## tausq[upper.tri(tausq, diag = TRUE)] <-
    ##     theta[(p + 1):((p + 1) + .5*(p * ( p  + 1 )) - 1)]
    ## tausq[lower.tri(tausq)] <- tausq[upper.tri(tausq)]
    ## sigsq <- theta[((p + 1) + .5*(p * ( p  + 1 )))]
    ## phi   <- theta[((p + 1) + .5*(p * ( p  + 1 )) + 1)]

    if(apply_exp) {
        tausq <- exp(tausq)
        sigsq <- exp(sigsq)
        phi   <- exp(phi)
    }

    .n <- NROW(.dt)

    switch(model,
           "matern" = {
               varcov_u1 <- comp_mat_cov(cross_dists = dists,
                                         n = .n, n2 = .n,
                                         phi = phi,
                                         sigsq = sigsq,
                                         kappa = kappa)
           },
           "pexp" = {
               varcov_u1 <- comp_pexp_cov(cross_dists = dists,
                                          n = .n, n2 = .n,
                                          phi = phi,
                                          sigsq = sigsq,
                                          kappa = kappa)
           },
           "gaussian" = {
               varcov_u1 <- comp_gauss_cov(cross_dists = dists,
                                           n = .n, n2 = .n,
                                           phi = phi,
                                           sigsq = sigsq)
           },
           "spherical" = {
               varcov_u1 <- comp_spher_cov(cross_dists = dists,
                                           n = .n, n2 = .n,
                                           phi = phi,
                                           sigsq = sigsq)
           },
           "cs" = {
               varcov_u1 <- comp_cs_cov(cross_dists = dists,
                                        n = .n, n2 = .n,
                                        phi = phi,
                                        sigsq = sigsq)
           },
           "w1" = {
               varcov_u1 <- comp_w1_cov(cross_dists = dists,
                                        n = .n, n2 = .n,
                                        phi = phi,
                                        sigsq = sigsq)
           },
           "tapmat" = {
               varcov_u1 <- comp_tapmat_cov(cross_dists = dists,
                                            n = .n, n2 = .n,
                                            phi = phi,
                                            sigsq = sigsq,
                                            kappa = kappa,
                                            theta = tr)
           })
    
    varcov_y  <- varcov_u1 + diag(tausq / npix,
                                  nrow = .n, ncol = .n)
    
    log_lik_y <- mvtnorm::dmvnorm(x = matrix(.dt, nrow = 1),
                                  mean  = matrix(rep(mu, .n),
                                                 ncol = 1),
                                  sigma = varcov_y,
                                  log = TRUE,
                                  checkSymmetry = FALSE)

    ## varcov_y  <- kronecker(tcrossprod(matrix(rep(1, p), ncol = 1)), varcov_u1) +
    ##     kronecker(tausq, diag(1 / npix, nrow = .n, ncol = .n)) 
    ## log_lik_y <- mvtnorm::dmvnorm(x = matrix(c(.dt), nrow = 1),
    ##                               mean  = c(kronecker(mu,
    ##                                                   matrix(rep(1, .n), ncol = 1))),
    ##                               sigma = varcov_y,
    ##                               log = TRUE,
    ##                               checkSymmetry = FALSE)

    return( - log_lik_y )
}

##' Evaluate the log-likelihood for a given set of parameters - New
##' parametrization + profile likelihood
##'
##' Internal use.
##' @title Evaluate log-lik
##' @param theta a \code{vector} of size 2 containing the parameters associated
##'     with the model. These parameters are \eqn{\nu} and \eqn{\phi},
##'     respectively.
##' @param .dt a \code{numeric} vector containing the variable \eqn{Y}.
##' @param dists a \code{list} of size three. The first containing the distance
##'     matrices associated with the regions where \eqn{Y} was measured, the
##'     second for the distance matrices associated with \eqn{X}, and the last
##'     containing the cross-distance matrices.
##' @param npix a \code{integer vector} containing the number of pixels within
##'     each polygon. (Ordered by the id variables for the polygons).
##' @param model a \code{character} indicating which covariance function to
##'     use. Possible values are \code{c("matern", "pexp", "gaussian",
##'     "spherical", "cs", "w1", "tapmat")}.
##' @param kappa \eqn{\kappa} parameter. Not necessary if \code{mode} is
##'     \code{"gaussian"} or \code{"spherical"}
##' @param tr \eqn{\theta_r} taper range.
##' @param apply_exp a \code{logical} indicater wheter the exponential
##'     transformation should be applied to variance parameters. This
##'     facilitates the optimization process.
##' 
##' @return a scalar representing \code{-log.lik}.
##' @keywords internal
singl_log_plik <- function(theta, .dt, dists, npix, model,
                           kappa = NULL, tr = NULL,
                           apply_exp = FALSE) {
    
    if(! apply_exp & any(theta < 0 )) {
        return(NA_real_)
    }
    
    nu  <- theta[1]
    phi <- theta[2]
    if(apply_exp) {
        nu  <- exp(nu)
        phi <- exp(phi)
    }
    .n <- NROW(.dt)
    switch(model,
           "matern" = {
               varcov_u1 <- comp_mat_cov(cross_dists = dists,
                                         n = .n, n2 = .n,
                                         phi = phi,
                                         sigsq = 1,
                                         kappa = kappa)
           },
           "pexp" = {
               varcov_u1 <- comp_pexp_cov(cross_dists = dists,
                                          n = .n, n2 = .n,
                                          phi = phi,
                                          sigsq = 1,
                                          kappa = kappa)
           },
           "gaussian" = {
               varcov_u1 <- comp_gauss_cov(cross_dists = dists,
                                           n = .n, n2 = .n,
                                           phi = phi,
                                           sigsq = 1)
           },
           "spherical" = {
               varcov_u1 <- Matrix(
                   comp_spher_cov(cross_dists = dists,
                                  n = .n, n2 = .n,
                                  phi = phi,
                                  sigsq = 1),
                   sparse = TRUE
               )
           },
           "cs" = {
               varcov_u1 <- Matrix(
                   comp_cs_cov(cross_dists = dists,
                               n = .n, n2 = .n,
                               phi = phi,
                               sigsq = 1),
                   sparse = TRUE
               )
           },
           "w1" = {
               varcov_u1 <- Matrix(
                   comp_w1_cov(cross_dists = dists,
                               n = .n, n2 = .n,
                               phi = phi,
                               sigsq = 1),
                   sparse = TRUE
               )
           },
           "tapmat" = {
               varcov_u1 <- Matrix(
                   comp_tapmat_cov(cross_dists = dists,
                                   n = .n, n2 = .n,
                                   phi = phi,
                                   sigsq = 1,
                                   kappa = kappa,
                                   theta = tr),
                   sparse = TRUE
               )
           })
    
    
    V <- varcov_u1 + diag(nu / npix,
                          nrow = .n, ncol = .n)
    chol_v <- try(chol(V))
    if(inherits(chol_v, "try-error")) {
        inv_v <- solve(V)
        mles   <- est_mle(.dt, inv_v)
        log_lik_y <- .5 * (.n * log(2 * pi) + .n * log(mles[length(mles)]) +
                           2 * log(det(V)) + .n)
    } else {
        inv_v  <- chol2inv(chol_v)
        mles   <- est_mle(.dt, inv_v)
        log_lik_y <- .5 * (.n * log(2 * pi) + .n * log(mles[length(mles)]) +
                           2 * sum(log(diag(chol_v))) + .n)
    }
    return( log_lik_y )
}

##' Evaluate the log-likelihood for a given set of parameters - No nugget +
##' profile likelihood
##'
##' Internal use.
##' @title Evaluate log-lik
##' @param theta a scalar for the \eqn{\phi} parameter.
##' @param .dt a \code{numeric} vector containing the variable \eqn{Y}.
##' @param dists a \code{list} of size three. The first containing the distance
##'     matrices associated with the regions where \eqn{Y} was measured, the
##'     second for the distance matrices associated with \eqn{X}, and the last
##'     containing the cross-distance matrices.
##' @param npix a \code{integer vector} containing the number of pixels within
##'     each polygon. (Ordered by the id variables for the polygons).
##' @param model a \code{character} indicating which covariance function to
##'     use. Possible values are \code{c("matern", "pexp", "gaussian",
##'     "spherical", "cs", "w1", "tapmat")}.
##' @param kappa \eqn{\kappa} parameter. Not necessary if \code{mode} is
##'     \code{"gaussian"} or \code{"spherical"}
##' @param apply_exp a \code{logical} indicater wheter the exponential
##'     transformation should be applied to variance parameters. This
##'     facilitates the optimization process.
##' 
##' @return a scalar representing \code{-log.lik}.
##' @keywords internal
singl_log_lik_nn <- function(theta, .dt, dists, npix, model,
                             kappa = NULL,
                             tr = NULL, apply_exp = FALSE) {
    
    if(! apply_exp & theta < 0 ) {
        return(NA_real_)
    }
    
    phi <- theta

    if(apply_exp) {
        phi <- exp(phi)
    }

    .n <- NROW(.dt)

    switch(model,
           "matern" = {
               varcov_u1 <- comp_mat_cov(cross_dists = dists,
                                         n = .n, n2 = .n,
                                         phi = phi,
                                         sigsq = 1,
                                         kappa = kappa)
           },
           "pexp" = {
               varcov_u1 <- comp_pexp_cov(cross_dists = dists,
                                          n = .n, n2 = .n,
                                          phi = phi,
                                          sigsq = 1,
                                          kappa = kappa)
           },
           "gaussian" = {
               varcov_u1 <- comp_gauss_cov(cross_dists = dists,
                                           n = .n, n2 = .n,
                                           phi = phi,
                                           sigsq = 1)
           },
           "spherical" = {
               varcov_u1 <- Matrix(
                   comp_spher_cov(cross_dists = dists,
                                  n = .n, n2 = .n,
                                  phi = phi,
                                  sigsq = 1),
                   sparse = TRUE
               )
           },
           "cs" = {
               varcov_u1 <- Matrix(
                   comp_cs_cov(cross_dists = dists,
                               n = .n, n2 = .n,
                               phi = phi,
                               sigsq = 1),
                   sparse = TRUE
               )
           },
           "w1" = {
               varcov_u1 <- Matrix(
                   comp_w1_cov(cross_dists = dists,
                               n = .n, n2 = .n,
                               phi = phi,
                               sigsq = 1),
                   sparse = TRUE
               )
           },
           "tapmat" = {
               varcov_u1 <- Matrix(
                   comp_tapmat_cov(cross_dists = dists,
                                   n = .n, n2 = .n,
                                   phi = phi,
                                   sigsq = 1,
                                   kappa = kappa,
                                   theta = tr),
                   sparse = TRUE
               )
           })
    
    V <- varcov_u1
    chol_v <- try(chol(V))
    if(inherits(chol_v, "try-error")) {
        inv_v <- solve(V)
        mles <- est_mle(.dt, inv_v)
        log_lik_y <- .5 * (.n * log(2 * pi) + .n * log(mles[length(mles)]) +
                           2 * log(det(V)) + .n)
    } else {
        inv_v <- chol2inv(chol_v)
        mles <- est_mle(.dt, inv_v)
        log_lik_y <- .5 * (.n * log(2 * pi) + .n * log(mles[length(mles)]) +
                           2 * sum(log(diag(chol_v))) + .n)
    }
    
    return( log_lik_y )
}

## ##' Evaluate the log-likelihood for a given set of parameters - New
## ##' parametrization + restricted likelihood
## ##'
## ##' Internal use.
## ##' @title Evaluate log-lik
## ##' @param theta a \code{list} of size 3 containing the parameters associated
## ##'     with the model. (Explain why)
## ##' @param .dt a \code{numeric} vector containing the variable \eqn{Y}.
## ##' @param X a \code{numeric} design matrix \eqn{X} of dimension \eqn{n \times
## ##'     p}.
## ##' @param dists a \code{list} of size three. The first containing the distance
## ##'     matrices associated with the regions where \eqn{Y} was measured, the
## ##'     second for the distance matrices associated with \eqn{X}, and the last
## ##'     containing the cross-distance matrices.
## ##' @param npix a \code{integer vector} containing the number of pixels within
## ##'     each polygon. (Ordered by the id variables for the polygons).
## ##' @param model a \code{character} indicating which covariance function to
## ##'     use. Possible values are \code{c("matern", "pexp", "gaussian",
## ##'     "spherical")}.
## ##' @param kappa \eqn{\kappa} parameter. Not necessary if \code{mode} is
## ##'     \code{"gaussian"} or \code{"spherical"}
## ##' @param apply_exp a \code{logical} indicater wheter the exponential
## ##'     transformation should be applied to variance parameters. This
## ##'     facilitates the optimization process.
## ##' 
## ##' @return a scalar representing \code{-log.lik}.
## singl_log_rel <- function(theta, .dt, dists, npix, model,
##                           kappa = NULL, apply_exp = FALSE) {    
##     if(! apply_exp & any(theta[(p + 1):length(theta)] < 0 )) {
##         return(NA_real_)
##     }
##     nu    <- theta[1]
##     sigsq <- theta[2]
##     phi   <- theta[3]
##     ## nu <- matrix(nrow = p, ncol = p)
##     ## nu[upper.tri(nu, diag = TRUE)] <-
##     ##     theta[1:(1 + .5*(p * ( p  + 1 )) - 1)]
##     ## nu[lower.tri(nu)] <- nu[upper.tri(nu)]
##     ## sigsq <- theta[( .5 * (p * ( p  + 1 )) + 1)]
##     ## phi   <- theta[( .5 * (p * ( p  + 1 )) + 2)]
##     if(apply_exp) {
##         nu    <- exp(nu)
##         sigsq <- exp(sigsq)
##         phi   <- exp(phi)
##     }
##     .n <- NROW(.dt)
##     switch(model,
##            "matern" = {
##                varcov_u1 <- comp_mat_cov(cross_dists = dists,
##                                          n = .n, n2 = .n,
##                                          phi = phi,
##                                          sigsq = 1,
##                                          kappa = kappa)
##            },
##            "pexp" = {
##                varcov_u1 <- comp_pexp_cov(cross_dists = dists,
##                                           n = .n, n2 = .n,
##                                           phi = phi,
##                                           sigsq = 1,
##                                           kappa = kappa)
##            },
##            "gaussian" = {
##                varcov_u1 <- comp_gauss_cov(cross_dists = dists,
##                                            n = .n, n2 = .n,
##                                            phi = phi,
##                                            sigsq = 1)
##            },
##            "spherical" = {
##                varcov_u1 <- comp_spher_cov(cross_dists = dists,
##                                            n = .n, n2 = .n,
##                                            phi = phi,
##                                            sigsq = 1)
##            })
##     if(p == 1) {
##         V <- varcov_u1 + diag(nu / npix,
##                               nrow = .n, ncol = .n)
##         chol_v <- chol(V)
##         inv_v  <- chol2inv(chol_v)
##         ones_n <- matrix(rep(1, .n), ncol = 1L)
##         y <- matrix(.dt, ncol = 1L)
##         mu_hat    <- as.numeric( crossprod(ones_n, inv_v) %*% y / (sum(inv_v)) )
##         y_mu_hat  <- y - mu_hat
##         log_lik_y <- .5 * (.n * log(2 * pi) + .n * log(sigsq) + 2 * sum(log(diag(chol_v))) +
##                            (crossprod(y_mu_hat, inv_v) %*% y_mu_hat / sigsq) +
##                            (sum(inv_v) / sigsq))
##     } else {
##         stop("to be implemented.")
##     }
##     return( log_lik_y )
## }

##' Evaluate the log-likelihood for a given set of parameters
##'
##' Internal use.
##' @title Evaluate log-lik
##' @param theta a \code{numeric} vector of size \eqn{3} containing the
##'     parameters values associated with \eqn{\mu}, \eqn{\sigma^2}, and
##'     \eqn{\phi}, respectively.
##' @param .dt a \code{numeric} vector containing the variable \eqn{Y}.
##' @param dists a \code{list} of size three. The first containing the distance
##'     matrices associated with the regions where \eqn{Y} was measured, the
##'     second for the distance matrices associated with \eqn{X}, and the last
##'     containing the cross-distance matrices.
##' @param npix a \code{integer vector} containing the number of pixels within
##'     each polygon. (Ordered by the id variables for the polygons).
##' @param model a \code{character} indicating which covariance function to
##'     use. Possible values are \code{c("matern", "pexp", "gaussian",
##'     "spherical")}.
##' @param kappa \eqn{\kappa} parameter. Not necessary if \code{mode} is
##'     \code{"gaussian"} or \code{"spherical"}
##' @param tr taper range
##' @param apply_exp a \code{logical} indicater wheter the exponential
##'     transformation should be applied to variance parameters. This
##'     facilitates the optimization process.
##' 
##' @return a scalar representing \code{-log.lik}.
##' @keywords internal
singl_ll_nn_hess <- function(theta, .dt, dists, npix, model,
                             kappa = NULL, tr = NULL,
                             apply_exp = FALSE) {
    npar <- length(theta)
    
    mu    <- theta[1]
    sigsq <- theta[2]
    phi   <- theta[3]
    ## mu    <- matrix(theta[1:p], ncol = 1)
    ## sigsq <- theta[p + 1]
    ## phi   <- theta[p + 2]
    if(apply_exp) {
        sigsq <- exp(sigsq)
        phi   <- exp(phi)
    }

    .n <- NROW(.dt)

    switch(model,
           "matern" = {
               varcov_u1 <- comp_mat_cov(cross_dists = dists,
                                         n = .n, n2 = .n,
                                         phi = phi,
                                         sigsq = sigsq,
                                         kappa = kappa)
           },
           "pexp" = {
               varcov_u1 <- comp_pexp_cov(cross_dists = dists,
                                          n = .n, n2 = .n,
                                          phi = phi,
                                          sigsq = sigsq,
                                          kappa = kappa)
           },
           "gaussian" = {
               varcov_u1 <- comp_gauss_cov(cross_dists = dists,
                                           n = .n, n2 = .n,
                                           phi = phi,
                                           sigsq = sigsq)
           },
           "spherical" = {
               varcov_u1 <- comp_spher_cov(cross_dists = dists,
                                           n = .n, n2 = .n,
                                           phi = phi,
                                           sigsq = sigsq)
           },
           "w1" = {
               varcov_u1 <- comp_w1_cov(cross_dists = dists,
                                        n = .n, n2 = .n,
                                        phi = phi,
                                        sigsq = sigsq)
           },
           "cs" = {
               varcov_u1 <- comp_cs_cov(cross_dists = dists,
                                        n = .n, n2 = .n,
                                        phi = phi,
                                        sigsq = sigsq)
           },
           "tapmat" = {
               varcov_u1 <- comp_tapmat_cov(cross_dists = dists,
                                            n = .n, n2 = .n,
                                            phi = phi,
                                            sigsq = sigsq,
                                            kappa = kappa,
                                            theta = tr)
           })
    log_lik_y <- mvtnorm::dmvnorm(x = matrix(.dt, nrow = 1),
                                  mean  = matrix(rep(mu, .n),
                                                 ncol = 1),
                                  sigma = varcov_u1,
                                  log = TRUE,
                                  checkSymmetry = FALSE)
    return( - log_lik_y )
}
