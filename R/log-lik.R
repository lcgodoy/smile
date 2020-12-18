##' Evaluate the log-likelihood for a given set of parameters
##'
##' Internal use.
##' @title Evaluate log-lik
##' @param theta a \code{list} of size \eqn{2 p + \frac{p(p + 1)}{2} + 3}
##'     containing the parameters associated with the model. (Explain why)
##' @param .dt a \code{list} with two positions. The first containing the
##'     numerical data for the variable \eqn{Y}, and the second for \eqn{X}.
##' @param dists a \code{list} of size three. The first containing the distance
##'     matrices associated with the regions where \eqn{Y} was measured, the
##'     second for the distance matrices associated with \eqn{X}, and the last
##'     containing the cross-distance matrices.
##' @param model a \code{character} indicating which covariance function to
##'     use. Possible values are \code{c("matern", "pexp", "gaussian",
##'     "spherical")}.
##' @param kappa \eqn{\kappa} parameter. Not necessary if \code{mode} is
##'     \code{"gaussian"} or \code{"spherical"}
##' @param apply_exp a \code{logical} indicater wheter the exponential
##'     transformation should be applied to variance parameters. This
##'     facilitates the optimization process.
##' 
##' @return a scalar representing \code{-log.lik}.
mult_log_lik <- function(theta, .dt, dists, model, kappa = NULL,
                         apply_exp = FALSE) {
    p <- NCOL(.dt[[1]])

    if(p == 1) {
        alpha <- theta[1]
        .beta  <- theta[2]
        omega <- theta[3]
        mu_x  <- theta[4]
        sigsq <- theta[5]
        phi   <- theta[6]
        nusq  <- theta[7]
    } else {
        alpha <- matrix(theta[1:p], ncol = 1)
        .beta <- matrix(theta[(p + 1):(2 * p)], ncol = 1)
        mu_x  <- theta[(2 * p) + 1]
        omega <- matrix(nrow = p, ncol = p)
        omega[upper.tri(omega, diag = TRUE)] <-
            theta[((2 * p) + 2):((2 * p) + 1 + .5*(p * (p  + 1 )))]
        omega[lower.tri(omega)] <- omega[upper.tri(omega)]
        sigsq <- theta[((2 * p) + 2 + .5*(p * (p  + 1 )))]
        phi   <- theta[((2 * p) + 3 + .5*(p * (p  + 1 )))]
        nusq  <- theta[((2 * p) + 4 + .5*(p * (p  + 1 )))]
    }
    if(apply_exp) {
        omega <- exp(omega)
        sigsq <- exp(sigsq)
        phi   <- exp(phi)
        nusq  <- exp(nusq)
    }
    
    .n    <- sapply(.dt, NROW)

    switch(model,
           "matern" = {
               if(is.null(kappa))
                   kappa <- .5
               varcov_u1 <- comp_mat_cov(cross_dists = dists[[1]],
                                         n = .n[[1]], n2 = .n[[1]],
                                         phi = phi,
                                         sigsq = sigsq,
                                         kappa = kappa)

               varcov_u2 <- comp_mat_cov(cross_dists = dists[[2]],
                                         n = .n[[2]], n2 = .n[[2]],
                                         phi = phi,
                                         sigsq = sigsq,
                                         kappa = kappa)
               
               varcov_u12 <- comp_mat_cov(cross_dists = dists[[3]],
                                          n = .n[[1]], n2 = .n[[2]],
                                          phi = phi, sigsq = sigsq,
                                          kappa  = kappa)
           },
           "pexp" = {
               if(is.null(kappa))
                   kappa <- 1
               varcov_u1 <- comp_pexp_cov(cross_dists = dists[[1]],
                                          n = .n[[1]], n2 = .n[[1]],
                                          phi = phi,
                                          sigsq = sigsq,
                                          kappa = kappa)

               varcov_u2 <- comp_pexp_cov(cross_dists = dists[[2]],
                                          n = .n[[2]], n2 = .n[[2]],
                                          phi = phi,
                                          sigsq = sigsq,
                                          kappa = kappa)
               
               varcov_u12 <- comp_pexp_cov(cross_dists = dists[[3]],
                                           n = .n[[1]], n2 = .n[[2]],
                                           phi = phi, sigsq = sigsq,
                                           kappa  = kappa)
           },
           "gaussian" = {
               varcov_u1 <- comp_gauss_cov(cross_dists = dists[[1]],
                                           n = .n[[1]], n2 = .n[[1]],
                                           phi = phi,
                                           sigsq = sigsq)

               varcov_u2 <- comp_gauss_cov(cross_dists = dists[[2]],
                                           n = .n[[2]], n2 = .n[[2]],
                                           phi = phi,
                                           sigsq = sigsq)
               
               varcov_u12 <- comp_gauss_cov(cross_dists = dists[[3]],
                                            n = .n[[1]], n2 = .n[[2]],
                                            phi = phi, sigsq = sigsq)
           },
           "spherical" = {
               varcov_u1 <- comp_spher_cov(cross_dists = dists[[1]],
                                           n = .n[[1]], n2 = .n[[1]],
                                           phi = phi,
                                           sigsq = sigsq)

               varcov_u2 <- comp_spher_cov(cross_dists = dists[[2]],
                                           n = .n[[2]], n2 = .n[[2]],
                                           phi = phi,
                                           sigsq = sigsq)
               
               varcov_u12 <- comp_spher_cov(cross_dists = dists[[3]],
                                            n = .n[[1]], n2 = .n[[2]],
                                            phi = phi, sigsq = sigsq)
           })

    if(p == 1) {
        varcov_y  <- (.beta * .beta * varcov_u1) + diag(omega, nrow = .n[[1]], ncol = .n[[1]])
        crossv_yx <- .beta * varcov_u12
    } else {
        varcov_y  <- kronecker(tcrossprod(.beta), varcov_u1) +
            kronecker(omega, diag(1, nrow = .n[[1]],
                                  ncol = .n[[1]]))
        crossv_yx <- kronecker(t(.beta), varcov_u12)
    }

    varcov_x  <- varcov_u2 + diag(nusq, nrow = .n[[2]], ncol = .n[[2]])

    ## dealing with possible inconsistencies during the numerical optimization
    inv_u2 <- tryCatch(
        chol2inv(chol(varcov_u2)),
        error = tryCatch(
            solve(varcov_u2),
            error = chol2inv(chol(sigsq*Matrix::nearPD(varcov_u2/sigsq)$mat))
        )
    )

    if(NCOL(.dt[[1]]) == 1) {
        mu_y_x <- matrix(rep(alpha, .n[[1]]), ncol = 1) +
            crossprod(crossv_yx, inv_u2) %*% matrix(.dt[[2]] - rep(mu_x, .n[[2]]), ncol = 1)
        varcov_y_x <- varcov_y - tcrossprod(crossv_yx %*% inv_u2, crossv_yx)
    } else {
        c_inv_u2 <- crossprod(crossv_yx, inv_u2)
        mu_y_x <- kronecker(matrix(alpha, ncol = 1), rep(1, NROW(.dt[[1]]))) +
            c_inv_u2 %*% matrix(.dt[[2]] - rep(mu_x, .n[[2]]), ncol = 1)
        varcov_y_x <- varcov_y - c_inv_u2 %*% crossv_yx
    }

    log_lik_x <- mvtnorm::dmvnorm(x = matrix(.dt[[2]], nrow = 1),
                                  mean = matrix(rep(mu_x, .n[[2]]), ncol = 1),
                                  sigma = varcov_x,
                                  log = TRUE,
                                  checkSymmetry = FALSE)
    
    log_lik_yx <- mvtnorm::dmvnorm(x = matrix(c(.dt[[1]]), nrow = 1),
                                   mean  = mu_y_x,
                                   sigma = varcov_y_x,
                                   log   = TRUE,
                                   checkSymmetry = FALSE)

    return(-(log_lik_x + log_lik_yx))
}

##' Evaluate the log-likelihood for a given set of parameters
##'
##' Internal use.
##' @title Evaluate log-lik
##' @param theta a \code{list} of size \eqn{2 p + \frac{p(p + 1)}{2} + 3}
##'     containing the parameters associated with the model. (Explain why)
##' @param .dt a \code{list} with two positions. The first containing the
##'     numerical data for the variable \eqn{Y}, and the second for \eqn{X}.
##' @param dists a \code{list} of size three. The first containing the distance
##'     matrices associated with the regions where \eqn{Y} was measured, the
##'     second for the distance matrices associated with \eqn{X}, and the last
##'     containing the cross-distance matrices.
##' @param model a \code{character} indicating which covariance function to
##'     use. Possible values are \code{c("matern", "pexp", "gaussian", "spherical")}.
##' @param kappa \eqn{\kappa} parameter. Not necessary if \code{mode} is
##'     \code{"gaussian"} or \code{"spherical"}
##' @param apply_exp a \code{logical} indicater wheter the exponential
##'     transformation should be applied to variance parameters. This
##'     facilitates the optimization process.
##' 
##' @return a scalar representing \code{-log.lik}.
singl_log_lik <- function(theta, .dt, dists, model, kappa = NULL,
                          apply_exp = FALSE) {

    p <- NCOL(.dt)

    if(p == 1) {
        alpha <- theta[1]
        omega <- theta[2]
        sigsq <- theta[3]
        phi   <- theta[4]
    } else {
        alpha <- matrix(theta[1:p], ncol = 1)
        omega <- matrix(nrow = p, ncol = p)
        omega[upper.tri(omega, diag = TRUE)] <-
            theta[(p + 1):((p + 1) + .5*(p * ( p  + 1 )) - 1)]
        omega[lower.tri(omega)] <- omega[upper.tri(omega)]
        sigsq <- theta[((p + 1) + .5*(p * ( p  + 1 )))]
        phi   <- theta[((p + 1) + .5*(p * ( p  + 1 )) + 1)]
    }
    if(apply_exp) {
        omega <- exp(omega)
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
           })
    
    if(p == 1) {
        varcov_y  <- varcov_u1 + diag(omega, nrow = .n, ncol = .n)
        
        log_lik_y <- mvtnorm::dmvnorm(x = matrix(.dt, nrow = 1),
                                      mean  = rep(alpha, .n),
                                      sigma = varcov_y,
                                      log = TRUE,
                                      checkSymmetry = FALSE)
    } else {
        varcov_y  <- kronecker(tcrossprod(matrix(rep(1, p), ncol = 1)), varcov_u1) +
            kronecker(omega, diag(1, nrow = .n, ncol = .n))
        
        log_lik_y <- mvtnorm::dmvnorm(x = matrix(c(.dt), nrow = 1),
                                      mean  = c(kronecker(alpha,
                                                          matrix(rep(1, .n), ncol = 1))),
                                      sigma = varcov_y,
                                      log = TRUE,
                                      checkSymmetry = FALSE)
    }

    return( - log_lik_y )
}
