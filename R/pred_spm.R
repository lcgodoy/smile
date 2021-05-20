##' @name predict_spm
##' @export
predict_spm <- function(x, ...) UseMethod("predict_spm")

##' @title Computing the matrices to make predictions.
##'
##' @description Computes the matrices needed to make predictions.
##' 
##' @param dists either a \code{numeric matrix} representing the (hausdorff)
##'     distance matrix between polygons or a \code{list} used to compute the
##'     "grid covariances".
##' @param cdists either a \code{numeric matrix} representing the (hausdorff)
##'     distance matrix between observed polygons and predicted locations or a
##'     \code{list} used to compute the "grid cross-covariances".
##' @param predists a \code{numeric matrix} representing the (hausdorff)
##'     distance matrix between the prediction locations.
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
##' @param phi a \code{numeric} value indicating the \eqn{\theta} parameter
##'     associated with the "decay" speed of the spatial dependence.
##' @param sigsq a \code{numeric} value indicating the \eqn{\sigma^2} parameter
##'     associated with the variance of a Gaussian Process.
##' @param n_obs number of observed polygons (used only when \code{method ==
##'     "grid"}).
##' @param n_pred number of prediction locations (used only when \code{method ==
##'     "grid"}).
##'
##' @name compute_mats
##' 
##' @return a \code{list} containing three covariance matrices associated with
##'     the observed data, cross covariance between observed data and locations
##'     on which we want to make the predictions, and the predicted locations.
comp_haus_mats <- function(dists, cdists, predists,
                           model, phi, sigsq, kappa) {

    out <- vector(mode = "list", length = 3L)
    names(out) <- c("sig_y", "d_mat", "sig_pred")
    
    if( model == "matern" ) {
        if(is.null(kappa))
            kappa <- .5
        
        out$sig_y <- mat_cov(dists,
                             phi   = phi,
                             sigsq = sigsq,
                             kappa = kappa)
        
        out$d_mat <- mat_cov(cdists,
                             phi   = phi,
                             sigsq = sigsq,
                             kappa = kappa)
        
        out$sig_pred <- mat_cov(predists,
                                phi   = phi,
                                sigsq = sigsq,
                                kappa = kappa)
        
    } else if( model == "pexp" ) {
        if(is.null(kappa))
            kappa <- 1
        
        out$sig_y <- pexp_cov(dists,
                              phi   = phi,
                              sigsq = sigsq,
                              kappa = kappa)
        
        out$d_mat <- pexp_cov(cdists,
                              phi   = phi,
                              sigsq = sigsq,
                              kappa = kappa)
        
        out$sig_pred <- pexp_cov(predists,
                                 phi   = phi,
                                 sigsq = sigsq,
                                 kappa = kappa)
        
    } else if( model == "gaussian" ) {
        out$sig_y <- gauss_cov(dists,
                               phi   = phi,
                               sigsq = sigsq)
        
        out$d_mat <- gauss_cov(cdists,
                               phi   = phi,
                               sigsq = sigsq)
        
        out$sig_pred <- gauss_cov(predists,
                                  phi   = phi,
                                  sigsq = sigsq)
    } else {
        out$sig_y <- spher_cov(dists, 
                               phi   = phi,
                               sigsq = sigsq)
        
        out$d_mat <- spher_cov(cdists,
                               phi   = phi,
                               sigsq = sigsq)
        
        out$sig_pred <- spher_cov(dists = predists,
                                  phi   = phi,
                                  sigsq = sigsq)
    }

    return(out)
}

##' @name compute_mats
comp_grid_mats <- function(dists, cdists, predists,
                           model, phi, sigsq,
                           kappa, n_obs, n_pred) {

    out <- vector(mode = "list", length = 3L)
    names(out) <- c("sig_y", "d_mat", "sig_pred")
    
    if( model == "matern" ) {
        if(is.null(kappa))
            kappa <- .5
        
        out$sig_y <- comp_mat_cov(dists,
                                  n = n_obs, n2 = n_obs,
                                  phi   = phi,
                                  sigsq = sigsq,
                                  kappa = kappa)
        
        
        out$d_mat <- comp_mat_cov(cross_dists = cdists,
                                  n = n_obs, n2 = n_pred,
                                  phi   = phi,
                                  sigsq = sigsq,
                                  kappa = kappa)
        
        out$sig_pred <- mat_cov(dists = predists,
                                phi   = phi,
                                sigsq = sigsq,
                                kappa = kappa)
        
    } else if( model == "pexp" ) {
        if(is.null(kappa))
            kappa <- 1
        
        out$sig_y <- comp_pexp_cov(dists,
                                   n = n_obs, n2 = n_obs,
                                   phi   = phi,
                                   sigsq = sigsq,
                                   kappa = kappa)
        
        out$d_mat <- comp_pexp_cov(cross_dists = cdists,
                               n = n_obs, n2 = n_pred,
                               phi   = phi,
                               sigsq = sigsq,
                               kappa = kappa)
        
        out$sig_pred <- pexp_cov(dists = predists,
                                 phi   = phi,
                                 sigsq = sigsq,
                                 kappa = kappa)
        
    } else if( model == "gaussian" ) {
        out$sig_y <- comp_gauss_cov(dists,
                                    n = n_obs, n2 = n_obs,
                                    phi   = phi,
                                    sigsq = sigsq)
        
        out$d_mat <- comp_gauss_cov(cross_dists = cdists,
                                    n = n_obs, n2 = n_pred,
                                    phi   = phi,
                                    sigsq = sigsq)
        
        out$sig_pred <- gauss_cov(dists = predists,
                                  phi   = phi,
                                  sigsq = sigsq)
    } else {
        out$sig_y <- comp_spher_cov(dists, 
                                    n = n_obs, n2 = n_obs,
                                    phi   = phi,
                                    sigsq = sigsq)
        
        out$d_mat <- comp_spher_cov(cross_dists = cdists,
                                    n = n_obs, n2 = n_pred,
                                    phi   = phi,
                                    sigsq = sigsq)
        
        out$sig_pred <- spher_cov(dists = predists,
                                  phi   = phi,
                                  sigsq = sigsq)
    }

    return(out)
}

##' @name predict_spm
##' @export
predict_spm.sspm_fit <- function(x, .aggregate = TRUE, ...) {
    p     <- NCOL(x$call_data$var)
    n_obs <- NROW(x$call_data$var)

    if(p == 1) {        
        mean_y <- matrix(rep(x$estimate["alpha"], n_obs), ncol = 1)
        if( x$call_data$method == "grid" ) {
            pred_grid <- x$call_data$grid
            ## create the distance matrix of the predictive location
            coords_pred <- sf::st_coordinates(pred_grid)
            ## u_pred      <- as.matrix(stats::dist(coords_pred))
            u_pred      <- distmat(coords_pred)
            n_pred      <- nrow(coords_pred)  # number of predicted location

            u_res_pred <- pred_cdist(get_grid_list(x_to_list = x$call_data$grid,
                                                   by = x$call_data$ids_var),
                                     coords_pred)
            
            var_lists <- comp_grid_mats(x$call_data$dists,
                                        u_res_pred,
                                        u_pred,
                                        x$model,
                                        x$estimate["phi"],
                                        x$estimate["sigsq"],
                                        x$kappa,
                                        n_obs, n_pred)
        } else {
            warning("when {method == 'hausdorff}' and a sf object is not specified for making the predictions, the predictions are done in a 3000 points grid.")
            ## create the distance matrix of the predictive location
            pred_grid  <- 
                sf::st_sample(x = sf::st_cast(sf::st_geometry(x$call_data$sf_poly),
                                              "POLYGON"),
                              size = 3000, 
                              type = "regular",
                              by_polygon = FALSE)
            
            coords_pred <- sf::st_coordinates(pred_grid)
            ## u_pred      <- as.matrix(stats::dist(coords_pred))
            u_pred      <- distmat(coords_pred)
            n_pred      <- nrow(coords_pred)  # number of predicted location
            
            u_res_pred <- cdist_haus(poly2coords(x$call_data$sf_poly),
                                     coords_pred)
            
            var_lists <- comp_haus_mats(x$call_data$dists,
                                        u_res_pred,
                                        u_pred,
                                        x$model,
                                        x$estimate["phi"],
                                        x$estimate["sigsq"],
                                        x$kappa)
        }

        var_lists$sig_y <- var_lists$sig_y +
            diag(x$estimate["omega"] / x$call_data$nap,
                 nrow = n_obs, ncol = n_obs)

        var_lists$sig_pred <- var_lists$sig_pred +
            diag(x$estimate["omega"],
                 nrow = nrow(var_lists$sig_pred),
                 ncol = ncol(var_lists$sig_pred))
        
        sig_y_inv <- chol2inv( chol(var_lists$sig_y) )

        mean_pred <- matrix(rep(x$estimate["alpha"], n_pred),
                            ncol = 1)

        dt_yinv  <- crossprod(var_lists$d_mat, sig_y_inv)

        sig_pred_y <- var_lists$sig_pred - (dt_yinv %*% var_lists$d_mat)

        mean_pred_y <- mean_pred +
            dt_yinv %*% (matrix(x$call_data$var - x$estimate["alpha"],
                                ncol = 1))

        if(any(diag(sig_pred_y) < 0)) {
            warning("Negative variance for at least one predicted region. Taking absolute value.")
            se_pred_y <- sqrt(abs(diag(sig_pred_y)))
        } else {
            se_pred_y <- sqrt(diag(sig_pred_y))
        }

        pred_grid <- sf::st_sf(data.frame("mu_pred"  = mean_pred_y,
                                          "se_pred"  = se_pred_y,
                                          "geometry" = pred_grid),
                               crs = sf::st_crs(pred_grid))

        if(.aggregate) {
            out_poly <- aggregate_aux(x = pred_grid[c("mu_pred", "se_pred")],
                                      by = x$call_data$sf_poly,
                                      FUN = mean,
                                      join = sf::st_intersects)
            
            output <-
                list(
                    mu_pred   = mean_pred_y,
                    sig_pred  = sig_pred_y,
                    pred_grid = pred_grid,
                    pred_agg  = out_poly
                )
        } else {
            output <-
                list(
                    mu_pred   = mean_pred_y,
                    sig_pred  = sig_pred_y,
                    pred_grid = pred_grid,
                    pred_agg  = NA
                )
        }

        
    } else {
        ## URGENT - not really
        stop("still has to be implemented.")
    }

    class(output) <- append(class(output), "spm_pred")
    
    return(output)
}

##' @name predict_spm
##' @export
predict_spm.mspm_fit <- function(x, ...) {
    stop("still has to be implemented.")
}

##' @name predict_spm
##'
##' @title Prediction over the same or a different set of regions (or points).
##'
##' @description Realizes predictions that can be useful when researchers are
##'     interested in predict a variable observed in one political division of a
##'     city (or state) on another division of the same region.
##' 
##' @param x a \code{sf} object such that its geometris are either points or
##'     polygons.
##' @param spm_obj an object of either class \code{sspm_fit} or \code{mspm_fit}
##' @param .aggregate \code{logical}. Should the predictions be aggregated? In
##'     case the input is only a "fit" object, the aggregation is made over the
##'     polygons on which the original data was observed. In case the input
##'     \code{x} is composed by \code{sf POLYGONS}, the aggregation is made over
##'     this new partition of the study region.
##' @param n_pts a \code{numeric} scalar standing for number of points to form a
##'     grid over the whole region to make the predictions
##' @param type \code{character} type of grid to be generated. See
##'     \code{\link[sf]{st_sample}}.
##' @param ... additional parameters
##'
##' @return an object of class \code{spm_pred}
##' 
##' @export
predict_spm.sf <- function(x, spm_obj, .aggregate = TRUE,
                           n_pts, type,
                           ...) {
    if(inherits(spm_obj, "mspm_fit"))
        stop("yet to be implemented.")

    if(sf::st_crs(x) != sf::st_crs(spm_obj$call_data$sf_poly)) {
        warning("`x` and the data on which the model was ajdusted are not in the same CRS. Reprojecting `x`")
        x <- sf::st_transform(x, sf::st_crs(spm_obj$call_data$sf_poly))
    }

    p     <- NCOL(spm_obj$call_data$var)
    n_obs <- NROW(spm_obj$call_data$var)

    ## HIGH PRIORITY AND EASY TO BE FIXED
    if(p > 1)
        stop("predictions for more than one random variable still have to be implemented.")

    stopifnot(all(grepl("(POLYGON|POINT)", sf::st_geometry_type(x))))

    if( all(grepl("POINT", sf::st_geometry_type(x))) ) {
        coords_pred <- sf::st_coordinates(x)
        pred_grid   <- sf::st_geometry(x)
        if( ! missing(n_pts) | ! missing(type) )
            warning("The arguments 'n_pts' and 'type' are ignored when the sf geometry type is of type POINT.")
    } else {
        if( spm_obj$call_data$method == "grid" ) {
            pred_grid   <- sf::st_sample(x    = sf::st_union(spm_obj$call_data$sf_poly),
                                         size = n_pts, 
                                         by_polygon = FALSE,
                                         type = type)
            coords_pred <- sf::st_coordinates(pred_grid)
        }
    }

    if( spm_obj$call_data$method == "grid" ) {
        u_pred <- distmat(coords_pred)
        n_pred <- nrow(coords_pred)  # number of locations to make predictions
        
        u_res_pred <- pred_cdist(get_grid_list(x_to_list = spm_obj$call_data$grid,
                                               by = spm_obj$call_data$ids_var),
                                 coords_pred)

        var_lists <- comp_grid_mats(spm_obj$call_data$dists,
                                    u_res_pred,
                                    u_pred,
                                    spm_obj$model,
                                    spm_obj$estimate["phi"],
                                    spm_obj$estimate["sigsq"],
                                    spm_obj$kappa,
                                    n_obs, n_pred)
    } else {
        if( all(grepl("POINT", sf::st_geometry_type(x))) ) {
            warning("when {method == 'hausdorff}' and a sf object is not specified for making the predictions, the predictions are done in a 3000 points grid.")
            ## create the distance matrix of the predictive location
            pred_grid   <- 
                sf::st_sample(x = sf::st_cast(sf::st_geometry(spm_obj$call_data$sf_poly),
                                              "POLYGON"),
                              size = 3000, 
                              type = "regular",
                              by_polygon = FALSE)
            
            coords_pred <- sf::st_coordinates(pred_grid)
            u_pred      <- distmat(coords_pred)
            n_pred      <- nrow(coords_pred)  # number of predicted location
            
            u_res_pred <- cdist_haus(poly2coords(spm_obj$call_data$sf_poly),
                                     coords_pred)
        } else {
            pred_poly <- poly2coords(x)
            u_pred    <- dist_haus(pred_poly)
            n_pred    <- nrow(x)  # number of predicted location
            
            u_res_pred <- cdist_haus_lists(poly2coords(spm_obj$call_data$sf_poly),
                                           pred_poly)
        }
        var_lists <- comp_haus_mats(spm_obj$call_data$dists,
                                    u_res_pred,
                                    u_pred,
                                    spm_obj$model,
                                    spm_obj$estimate["phi"],
                                    spm_obj$estimate["sigsq"],
                                    spm_obj$kappa)
    }
    
    if(all(grepl("POINT", sf::st_geometry_type(x))) & .aggregate) {
        warning("If you want to make predictions only for a set of locations, it does not make sense to use `.aggregate`.")
    }

    var_lists$sig_y <- var_lists$sig_y +
        diag(spm_obj$estimate["omega"] / spm_obj$call_data$nap,
             nrow = n_obs, ncol = n_obs)

    if( spm_obj$call_data$method == "hausdorff" & ! any(grepl("POINT", sf::st_geometry_type(x)))) {
        nap_pred <- as.numeric( sf::st_area(x) )
        var_lists$sig_pred <- var_lists$sig_pred +
            diag(spm_obj$estimate["omega"] / nap_pred,
                 nrow = nrow(var_lists$sig_pred),
                 ncol = ncol(var_lists$sig_pred))
    } else {
        var_lists$sig_pred <- var_lists$sig_pred +
            diag(spm_obj$estimate["omega"],
                 nrow = nrow(var_lists$sig_pred),
                 ncol = ncol(var_lists$sig_pred))
    }
    
    ## sig_y_inv <- solve(sig_y)
    sig_y_inv <- chol2inv( chol(var_lists$sig_y) )
    
    mean_y <- matrix(rep(spm_obj$estimate["alpha"], n_obs), ncol = 1)

    mean_pred <- matrix(rep(spm_obj$estimate["alpha"], n_pred),
                        ncol = 1)
    
    dt_yinv  <- crossprod(var_lists$d_mat, sig_y_inv)
    
    sig_pred_y <- var_lists$sig_pred - (dt_yinv %*% var_lists$d_mat)
    
    mean_pred_y <- mean_pred +
        dt_yinv %*% (spm_obj$call_data$var - mean_y)

    if(any(diag(sig_pred_y) < 0)) {
        warning("Negative variance for at least one predicted region. Taking absolute value.")
        se_pred_y <- sqrt(abs(diag(sig_pred_y)))
    } else {
        se_pred_y <- sqrt(diag(sig_pred_y))
    }

    if( spm_obj$call_data$method == "hausdorff" & ! any(grepl("POINT", sf::st_geometry_type(x)))) {
        warning(".aggregate is ignored when x's geometry type is POLYGON and method = 'hausdorff'.")        
        pred_grid <- sf::st_sf(data.frame("mu_pred"  = mean_pred_y,
                                          "se_pred"  = se_pred_y,
                                          "geometry" = sf::st_geometry(x)),
                               crs = sf::st_crs(x))
        output <-
            list(
                mu_pred   = mean_pred_y,
                sig_pred  = sig_pred_y,
                pred_grid = NA,
                pred_agg  = pred_grid
            )
        
    } else {
        pred_grid <- sf::st_sf(data.frame("mu_pred"  = mean_pred_y,
                                          "se_pred"  = se_pred_y,
                                          "geometry" = pred_grid),
                               crs = sf::st_crs(pred_grid))
        if(.aggregate) {
            out_poly <- aggregate_aux(x = pred_grid[c("mu_pred", "se_pred")],
                                      by = sf::st_geometry(x),
                                      FUN = mean,
                                      join = sf::st_intersects)
            
            output <-
                list(
                    mu_pred   = mean_pred_y,
                    sig_pred  = sig_pred_y,
                    pred_grid = pred_grid,
                    pred_agg  = out_poly
                )
        } else {
            output <-
                list(
                    mu_pred   = mean_pred_y,
                    sig_pred  = sig_pred_y,
                    pred_grid = pred_grid,
                    pred_agg  = NA
                )
        }
    }
    
    
    class(output) <- append(class(output), "spm_pred")

    return(output)
}

##' @name predict_spm
##' @export
predict_spm.matrix <- function(x, ...) {
    stop("yet to be implemented.")
}

##' @title Internal use only
##' @param x internal use
##' @param by internal use
##' @param FUN internal use
##' @param ... internal use
##' @param do_union internal use
##' @param simplify internal use
##' @param join internal use
##' @description taken from \code{\link[sf]{aggregate.sf}}.
aggregate_aux <- utils::getFromNamespace("aggregate.sf", "sf")
