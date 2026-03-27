##' Unified function for log-likelihood evaluation
##'
##' Internal use.
##' @title Evaluate log-lik
##' @param theta a named \code{numeric} vector containing parameters (\eqn{\mu,
##'   \sigma^2, \alpha, \phi}).
##' @param .dt a \code{numeric} vector containing the variable \eqn{Y}.
##' @param dists a \code{list} of distance matrices.
##' @param npix a \code{integer vector} containing the number of pixels within
##'   each polygon.
##' @param model a \code{character} indicating the covariance function family.
##' @param nu \eqn{\nu} parameter for the Matern or Powered Exponential
##'   covariance functions.
##' @param apply_exp a \code{logical} indicating whether to exponentiate
##'   non-negative parameters.
##' @param type a \code{character} specifying the likelihood type: `"full"` or
##'   `"profile"`.
##' 
##' @return a scalar representing \code{-log.lik}.
##' @keywords internal
log_lik_spm <- function(theta, .dt, dists, npix, model,
                        nu = NULL,
                        apply_exp = FALSE,
                        type = c("full", "profile")) {
  type <- match.arg(type)
  
  if (apply_exp) {
    to_exp <- names(theta) %in% c("sigsq", "phi", "al")
    theta[to_exp] <- exp(theta[to_exp])
  }
  
  if (!apply_exp) {
    if (any(theta[names(theta) %in% c("sigsq", "phi", "al")] < 0))
      return(NA_real_)
  }

  phi   <- theta["phi"]
  al    <- ifelse("al" %in% names(theta), theta["al"], 0)
  .n    <- NROW(.dt)

  V <- switch(model,
              "matern" = comp_mat_cov(cross_dists = dists,
                                      n = .n, n2 = .n,
                                      phi = phi,
                                      sigsq = 1,
                                      nu = ifelse(is.null(nu), 0.5, nu)),
              "pexp"   = comp_pexp_cov(cross_dists = dists,
                                       n = .n, n2 = .n,
                                       phi = phi,
                                       sigsq = 1,
                                       nu = ifelse(is.null(nu), 1.0, nu)),
              stop("Model '", model, "' not supported."))

  if (al > 0) {
    V <- V + diag(al / npix, nrow = .n, ncol = .n)
  }

  if (type == "profile") {
    chol_v <- try(chol(V), silent = TRUE)
    if (inherits(chol_v, "try-error")) {
      inv_v <- solve(V)
      mles <- est_mle(.dt, inv_v)
      log_lik_y <- .5 * (.n * log(2 * pi) + .n * log(mles["sigsq"]) + log(det(V)) + .n)
    } else {
      inv_v <- chol2inv(chol_v)
      mles <- est_mle(.dt, inv_v)
      log_lik_y <- .5 * (.n * log(2 * pi) + .n * log(mles["sigsq"]) + 2 * sum(log(diag(chol_v))) + .n)
    }
    return(as.numeric(log_lik_y))
  } else {
    mu <- theta["mu"]
    sigsq <- theta["sigsq"]

    log_lik_y <- mvtnorm::dmvnorm(x = matrix(.dt, nrow = 1),
                                  mean = matrix(rep(mu, .n), ncol = 1),
                                  sigma = sigsq * V,
                                  log = TRUE,
                                  checkSymmetry = FALSE)
    return(-as.numeric(log_lik_y))
  }
}

##' Backward compatibility wrapper for full log-likelihood
##' @keywords internal
singl_log_lik <- function(theta, .dt,
                          dists, npix,
                          model,
                          nu = NULL,
                          apply_exp = FALSE) {
  if (is.null(names(theta))) {
    if (length(theta) == 4) {
      names(theta) <- c("mu", "sigsq", "al", "phi")
    } else if (length(theta) == 3) {
      names(theta) <- c("mu", "sigsq", "phi")
    }
  }
  log_lik_spm(theta, .dt, dists, npix, model, nu, apply_exp, type = "full")
}

##' Backward compatibility wrapper for profile log-likelihood
##' @keywords internal
singl_log_plik <- function(theta, .dt, dists, npix, model, nu = NULL, apply_exp = FALSE) {
  if (is.null(names(theta))) {
    if (length(theta) == 2) {
      names(theta) <- c("al", "phi")
    } else if {
      (length(theta) == 1) names(theta) <- c("phi")
    }
  }
  log_lik_spm(theta, .dt, dists, npix, model, nu, apply_exp, type = "profile")
}

##' Backward compatibility wrapper for profile log-likelihood (no nugget)
##' @keywords internal
singl_log_lik_nn <- function(theta, .dt, dists, npix, model, nu = NULL, apply_exp = FALSE) {
  if (is.null(names(theta))) {
    names(theta) <- "phi"
  }
  log_lik_spm(theta, .dt, dists, npix, model, nu, apply_exp, type = "profile")
}

##' Backward compatibility wrapper for full log-likelihood (Hessian)
##' @keywords internal
singl_ll_nn_hess <- function(theta, .dt, dists, npix, model, nu = NULL, apply_exp = FALSE) {
  if (is.null(names(theta))) {
    if (length(theta) == 4) {
      names(theta) <- c("mu", "sigsq", "al", "phi")
    }
    else if (length(theta) == 3) {
      names(theta) <- c("mu", "sigsq", "phi")
    }
  }
  log_lik_spm(theta, .dt, dists, npix, model, nu, apply_exp, type = "full")
}
