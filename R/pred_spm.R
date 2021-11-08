##' @name predict_spm
##' @export
predict_spm <- function(x, ...) UseMethod("predict_spm")

##' @name predict_spm
##' @export
predict_spm.spm_fit <- function(x, .aggregate = TRUE, ...) {
    n_obs <- NROW(x$call_data$var)
    ids_betas <- which(grepl("^beta", names(x$estimate)))
    
    ## create the distance matrix of the predictive location
    coords_pred <- sf::st_coordinates(x$call_data$grid)
    ## u_pred      <- as.matrix(stats::dist(coords_pred))
    u_pred      <- distmat(coords_pred)
    n_pred      <- nrow(coords_pred)  # number of predicted location

    u_res_pred <- pred_cdist(get_grid_list(x_to_list = x$call_data$grid,
                                           by = x$call_data$ids_var),
                             coords_pred)

    ## can be turned in to a function to make to code cleaner
    switch(x$model,
           "matern" = {
               if(is.null(x$kappa))
                   x$kappa <- .5

               sig_y <- comp_mat_cov(x$call_data$dists,
                                     n = n_obs, n2 = n_obs,
                                     phi   = x$estimate["phi"],
                                     ## sigsq = x$estimate["sigsq"],
                                     sigsq = 1,
                                     kappa = x$kappa)
               d_mat <- comp_mat_cov(cross_dists = u_res_pred,
                                     n = n_obs, n2 = n_pred,
                                     phi   = x$estimate["phi"],
                                     sigsq = x$estimate["sigsq"],
                                     kappa = x$kappa)
               sig_pred <- mat_cov(dists = u_pred,
                                   phi   = x$estimate["phi"],
                                   ## sigsq = x$estimate["sigsq"],
                                   sigsq = 1,
                                   kappa = x$kappa)
           },
           "pexp" = {
               if(is.null(x$kappa))
                   x$kappa <- 1

               sig_y <- comp_pexp_cov(x$call_data$dists,
                                      n = n_obs, n2 = n_obs,
                                      phi   = x$estimate["phi"],
                                      ## sigsq = x$estimate["sigsq"],
                                      sigsq = 1,
                                      kappa = x$kappa)
               
               d_mat <- comp_pexp_cov(cross_dists = u_res_pred,
                                      n = n_obs, n2 = n_pred,
                                      phi   = x$estimate["phi"],
                                      sigsq = x$estimate["sigsq"],
                                      kappa = x$kappa)
               
               sig_pred <- pexp_cov(dists = u_pred,
                                    phi   = x$estimate["phi"],
                                    ## sigsq = x$estimate["sigsq"],
                                    sigsq = 1,
                                    kappa = x$kappa)
           },
           "gaussian" = {
               sig_y <- comp_gauss_cov(x$call_data$dists,
                                       n = n_obs, n2 = n_obs,
                                       phi   = x$estimate["phi"],
                                       ## sigsq = x$estimate["sigsq"]
                                       sigsq = 1)
               
               d_mat <- comp_gauss_cov(cross_dists = u_res_pred,
                                       n = n_obs, n2 = n_pred,
                                       phi   = x$estimate["phi"],
                                       sigsq = x$estimate["sigsq"])
               
               sig_pred <- gauss_cov(dists = u_pred,
                                     phi   = x$estimate["phi"],
                                     ## sigsq = x$estimate["sigsq"]
                                     sigsq = 1)
           },
           "spherical" = {
               sig_y <- comp_spher_cov(x$call_data$dists, 
                                       n = n_obs, n2 = n_obs,
                                       phi   = x$estimate["phi"],
                                       ## sigsq = x$estimate["sigsq"]
                                       sigsq = 1)
               
               d_mat <- comp_spher_cov(cross_dists = u_res_pred,
                                       n = n_obs, n2 = n_pred,
                                       phi   = x$estimate["phi"],
                                       sigsq = x$estimate["sigsq"])
               
               sig_pred <- spher_cov(dists = u_pred,
                                     phi   = x$estimate["phi"],
                                     ## sigsq = x$estimate["sigsq"]
                                     sigsq = 1)
           })

    if(length(x$estimate) > (length(ids_betas) + 2)) {
        if("tausq" %in% names(x$estimate)) {
            sig_y <- (x$estimate["sigsq"] * sig_y) +
                diag(x$estimate["tausq"] / x$call_data$npix,
                     nrow = n_obs, ncol = n_obs)
            
            sig_pred <- (x$estimate["sigsq"] * sig_pred) +
                diag(x$estimate["tausq"],
                     nrow = nrow(sig_pred),
                     ncol = ncol(sig_pred))
        } else if("nu" %in% names(x$estimate)) {
            sig_y <- x$estimate["sigsq"] *
                (sig_y +
                 diag(x$estimate["nu"] / x$call_data$npix,
                      nrow = n_obs, ncol = n_obs))
            sig_pred <- x$estimate["sigsq"] *
                (sig_pred +
                 diag(x$estimate["nu"],
                      nrow = nrow(sig_pred),
                      ncol = ncol(sig_pred)))
        } 
    } else {
        sig_y <- (x$estimate["sigsq"] * sig_y)
        sig_pred <- (x$estimate["sigsq"] * sig_pred) 
    }

    mean_y <- x$call_data$X %*% x$estimate[ids_betas]
    mean_pred <- x$call_data$X0 %*% x$estimate[ids_betas]

    sig_y_inv <- chol2inv(chol(sig_y))

    dt_yinv  <- crossprod(d_mat, sig_y_inv)

    sig_pred_y <- sig_pred - (dt_yinv %*% d_mat)

    mean_pred_y <- mean_pred +
        dt_yinv %*% (x$call_data$var - mean_y)

    if(any(diag(sig_pred_y) < 0)) {
        warning("Negative variance for at least one predicted region. Taking absolute value.")
        var_pred_y <- abs(diag(sig_pred_y))
    } else {
        var_pred_y <- diag(sig_pred_y)
    }
    
    pred_grid <- transform(x$call_data$grid,
                           mu_pred = mean_pred_y,
                           se_pred = var_pred_y)

    if(.aggregate) {
        out_poly <- aggregate_aux(x = pred_grid[c("mu_pred", "se_pred")],
                                  by = x$call_data$sf_poly,
                                  FUN = mean,
                                  join = sf::st_intersects)

        out_poly[["se_pred"]]  <- sqrt(out_poly[["se_pred"]])
        pred_grid[["se_pred"]] <- sqrt(pred_grid[["se_pred"]])
        
        output <-
            list(
                mu_pred   = mean_pred_y,
                sig_pred  = sig_pred_y,
                pred_grid = pred_grid,
                pred_agg  = out_poly
            )
    } else {
        pred_grid[["se_pred"]] <- sqrt(pred_grid[["se_pred"]])

        output <-
            list(
                mu_pred   = mean_pred_y,
                sig_pred  = sig_pred_y,
                pred_grid = pred_grid,
                pred_agg  = NA
            )
    }

    class(output) <- append(class(output), "spm_pred")
    
    return(output)
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
##' @param X0 a \code{matrix} representing the covariates observed at \code{x}.
##' @param spm_obj an object of either class \code{spm_fit} or \code{mspm_fit}
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
predict_spm.sf <- function(x, spm_obj, X0,
                           .aggregate = TRUE,
                           n_pts, type,
                           ...) {
    if(sf::st_crs(x) != sf::st_crs(spm_obj$call_data$sf_poly)) {
        warning("`x` and the data on which the model was ajdusted are not in the same CRS. Reprojecting `x`")
        x <- sf::st_transform(x, sf::st_crs(spm_obj$call_data$sf_poly))
    }

    n_obs <- NROW(spm_obj$call_data$var)

    stopifnot(all(grepl("(POLYGON|POINT)", sf::st_geometry_type(x))))

    if( all(grepl("POINT", sf::st_geometry_type(x))) ) {
        coords_pred <- sf::st_coordinates(x)
        pred_grid   <- sf::st_geometry(x)
        if( ! missing(n_pts) | ! missing(type) )
            warning("The arguments 'n_pts' and 'type' are ignored when the sf geometry type is POINT.")
    } else {
        pred_grid   <- sf::st_sample(x    = sf::st_union(spm_obj$call_data$sf_poly),
                                     size = n_pts, 
                                     by_polygon = FALSE,
                                     type = type)
        coords_pred <- sf::st_coordinates(pred_grid)
    }

    ids_betas <- which(grepl("^beta", names(spm_obj$estimate)))
    mean_y <- spm_obj$call_data$X %*% spm_obj$estimate[ids_betas]
    
    ## this part can be improved
    ## u_pred <- as.matrix(stats::dist(coords_pred))
    u_pred <- distmat(coords_pred)
    n_pred <- nrow(coords_pred)  # number of locations to make predictions

    u_res_pred <- pred_cdist(get_grid_list(x_to_list = spm_obj$call_data$grid,
                                           by = spm_obj$call_data$ids_var),
                             coords_pred)
    
    switch(spm_obj$model,
           "matern" = {
               if(is.null(spm_obj$kappa))
                   spm_obj$kappa <- .5

               sig_y <- comp_mat_cov(spm_obj$call_data$dists,
                                     n = n_obs, n2 = n_obs,
                                     phi   = spm_obj$estimate["phi"],
                                     ## sigsq = spm_obj$estimate["sigsq"],
                                     sigsq = 1,
                                     kappa = spm_obj$kappa)


               d_mat <- comp_mat_cov(cross_dists = u_res_pred,
                                     n = n_obs, n2 = n_pred,
                                     phi   = spm_obj$estimate["phi"],
                                     ## sigsq = spm_obj$estimate["sigsq"],
                                     sigsq = 1,
                                     kappa = spm_obj$kappa)
               
               sig_pred <- mat_cov(dists = u_pred,
                                   phi   = spm_obj$estimate["phi"],
                                   ## sigsq = spm_obj$estimate["sigsq"],
                                   sigsq = 1,
                                   kappa = spm_obj$kappa)
           },
           "pexp" = {
               if(is.null(spm_obj$kappa))
                   spm_obj$kappa <- 1

               sig_y <- comp_pexp_cov(spm_obj$call_data$dists,
                                      n = n_obs, n2 = n_obs,
                                      phi   = spm_obj$estimate["phi"],
                                      ## sigsq = spm_obj$estimate["sigsq"],
                                      sigsq = 1,
                                      kappa = spm_obj$kappa)
               
               d_mat <- comp_pexp_cov(cross_dists = u_res_pred,
                                      n = n_obs, n2 = n_pred,
                                      phi   = spm_obj$estimate["phi"],
                                      ## sigsq = spm_obj$estimate["sigsq"],
                                      sigsq = 1,
                                      kappa = spm_obj$kappa)
               
               sig_pred <- pexp_cov(dists = u_pred,
                                    phi   = spm_obj$estimate["phi"],
                                    ## sigsq = spm_obj$estimate["sigsq"],
                                    sigsq = 1,
                                    kappa = spm_obj$kappa)
           },
           "gaussian" = {
               sig_y <- comp_gauss_cov(spm_obj$call_data$dists,
                                       n = n_obs, n2 = n_obs,
                                       phi   = spm_obj$estimate["phi"],
                                       ## sigsq = spm_obj$estimate["sigsq"]
                                       sigsq = 1)
               
               d_mat <- comp_gauss_cov(cross_dists = u_res_pred,
                                       n = n_obs, n2 = n_pred,
                                       phi   = spm_obj$estimate["phi"],
                                       ## sigsq = spm_obj$estimate["sigsq"]
                                       sigsq = 1)
               
               sig_pred <- gauss_cov(dists = u_pred,
                                     phi   = spm_obj$estimate["phi"],
                                     ## sigsq = spm_obj$estimate["sigsq"]
                                     sigsq = 1)
           },
           "spherical" = {
               sig_y <- comp_spher_cov(spm_obj$call_data$dists, 
                                       n = n_obs, n2 = n_obs,
                                       phi   = spm_obj$estimate["phi"],
                                       ## sigsq = spm_obj$estimate["sigsq"]
                                       sigsq = 1)
               
               d_mat <- comp_spher_cov(cross_dists = u_res_pred,
                                       n = n_obs, n2 = n_pred,
                                       phi   = spm_obj$estimate["phi"],
                                       ## sigsq = spm_obj$estimate["sigsq"]
                                       sigsq = 1)
               
               sig_pred <- spher_cov(dists = u_pred,
                                     phi   = spm_obj$estimate["phi"],
                                     ## sigsq = spm_obj$estimate["sigsq"]
                                     sigsq = 1)
           })

    if(all(grepl("POINT", sf::st_geometry_type(x))) & .aggregate) {
        warning("If you want to make predictions only for a set of locations, it does not make sense to use `.aggregate`.")
    }

    if(length(spm_obj$estimate) > length(ids_betas) + 2) {
        if("tausq" %in% names(spm_obj$estimate)) {
            sig_y <- (spm_obj$estimate["sigsq"] * sig_y) +
                diag(spm_obj$estimate["tausq"] / spm_obj$call_data$npix,
                     nrow = n_obs, ncol = n_obs)
            
            sig_pred <- (spm_obj$estimate["sigsq"] * sig_pred) +
                diag(spm_obj$estimate["tausq"],
                     nrow = nrow(sig_pred),
                     ncol = ncol(sig_pred))
        } else if("nu" %in% names(spm_obj$estimate)) {
            sig_y <- spm_obj$estimate["sigsq"] *
                (sig_y +
                 diag(spm_obj$estimate["nu"] / spm_obj$call_data$npix,
                      nrow = n_obs, ncol = n_obs))
            
            sig_pred <- spm_obj$estimate["sigsq"] *
                (sig_pred +
                 diag(spm_obj$estimate["nu"],
                      nrow = nrow(sig_pred),
                      ncol = ncol(sig_pred)))
        } else {
            ## this part will be deleted soon
            sig_y <- (spm_obj$estimate["sigsq"] * sig_y) +
                diag(spm_obj$estimate["omega"] / spm_obj$call_data$npix,
                     nrow = n_obs, ncol = n_obs)
            
            sig_pred <- (spm_obj$estimate["sigsq"] * sig_pred) +
                diag(spm_obj$estimate["omega"],
                     nrow = nrow(sig_pred),
                     ncol = ncol(sig_pred))
        }
    } else {
        sig_y <- spm_obj$estimate["sigsq"] * sig_y
            
        sig_pred <- spm_obj$estimate["sigsq"] * sig_pred
    }
    
    ## sig_y_inv <- solve(sig_y)
    sig_y_inv <- chol2inv(chol(sig_y))
    
    ## mean_y <- matrix(rep(spm_obj$estimate["mu"], n_obs), ncol = 1)

    if(! is.null(X0) ) {
        if(NROW(X0) != n_pred)
            warning("X0 and prediction region are not conformable")
    } else {
        X0 <- matrix(1, nrow = n_pred, ncol = 1)
    }
    
    mean_pred <- X0 %*% spm_obj$estimate[ids_betas]
    
    dt_yinv  <- crossprod(d_mat, sig_y_inv)
    
    sig_pred_y <- sig_pred - (dt_yinv %*% d_mat)
    
    mean_pred_y <- mean_pred +
        dt_yinv %*% (spm_obj$call_data$var - mean_y)

    if(any(diag(sig_pred_y) < 0)) {
        warning("Negative variance for at least one predicted region. Taking absolute value.")
        var_pred_y <- abs(diag(sig_pred_y))
    } else {
        var_pred_y <- diag(sig_pred_y)
    }
    
    pred_grid <- transform(sf::st_sf(geometry = pred_grid),
                           mu_pred = mean_pred_y,
                           se_pred = var_pred_y)
    
    if(.aggregate) {
        out_poly <- aggregate_aux(x = pred_grid[c("mu_pred", "se_pred")],
                                  by = sf::st_geometry(x),
                                  FUN = mean,
                                  join = sf::st_intersects)

        out_poly[["se_pred"]]  <- sqrt(out_poly[["se_pred"]])
        pred_grid[["se_pred"]] <- sqrt(pred_grid[["se_pred"]])
        
        output <-
            list(
                mu_pred   = mean_pred_y,
                sig_pred  = sig_pred_y,
                pred_grid = pred_grid,
                pred_agg  = out_poly
            )
    } else {
        pred_grid[["se_pred"]] <- sqrt(pred_grid[["se_pred"]])
        
        output <-
            list(
                mu_pred   = mean_pred_y,
                sig_pred  = sig_pred_y,
                pred_grid = pred_grid,
                pred_agg  = NA
            )
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
