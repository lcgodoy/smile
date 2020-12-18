##' @name predict_spm
##' @export
predict_spm <- function(x, ...) UseMethod("predict_spm")

##' @name predict_spm
##' @export
predict_spm.sspm_fit <- function(x, .aggregate = TRUE, ...) {
    p     <- NCOL(x$call_data$var)
    n_obs <- NROW(x$call_data$var)

    if(p == 1) {        
        mean_y <- matrix(rep(x$estimate["alpha"], n_obs), ncol = 1)
        
        ## create the distance matrix of the predictive location
        coords_pred <- sf::st_coordinates(x$call_data$grid)
        u_pred      <- as.matrix(stats::dist(coords_pred))
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
                                         sigsq = x$estimate["sigsq"],
                                         kappa = x$kappa)


                   d_mat <- t(comp_mat_cov(cross_dists = u_res_pred,
                                           n = n_obs, n2 = n_pred,
                                           phi   = x$estimate["phi"],
                                           sigsq = x$estimate["sigsq"],
                                           kappa = x$kappa))
                   
                   sig_pred <- mat_cov(dists = u_pred,
                                       phi   = x$estimate["phi"],
                                       sigsq = x$estimate["sigsq"],
                                       kappa = x$kappa)
               },
               "pexp" = {
                   if(is.null(x$kappa))
                       x$kappa <- 1

                   sig_y <- comp_pexp_cov(x$call_data$dists,
                                          n = n_obs, n2 = n_obs,
                                          phi   = x$estimate["phi"],
                                          sigsq = x$estimate["sigsq"],
                                          kappa = x$kappa)
                   
                   d_mat <- t(comp_pexp_cov(cross_dists = u_res_pred,
                                            n = n_obs, n2 = n_pred,
                                            phi   = x$estimate["phi"],
                                            sigsq = x$estimate["sigsq"],
                                            kappa = x$kappa))
                   
                   sig_pred <- pexp_cov(dists = u_pred,
                                        phi   = x$estimate["phi"],
                                        sigsq = x$estimate["sigsq"],
                                        kappa = x$kappa)
               },
               "gaussian" = {
                   sig_y <- comp_gauss_cov(x$call_data$dists,
                                           n = n_obs, n2 = n_obs,
                                           phi   = x$estimate["phi"],
                                           sigsq = x$estimate["sigsq"])
                   
                   d_mat <- t(comp_gauss_cov(cross_dists = u_res_pred,
                                             n = n_obs, n2 = n_pred,
                                             phi   = x$estimate["phi"],
                                             sigsq = x$estimate["sigsq"]))
                   
                   sig_pred <- gauss_cov(dists = u_pred,
                                         phi   = x$estimate["phi"],
                                         sigsq = x$estimate["sigsq"])
               },
               "spherical" = {
                   sig_y <- comp_spher_cov(x$call_data$dists, 
                                           n = n_obs, n2 = n_obs,
                                           phi   = x$estimate["phi"],
                                           sigsq = x$estimate["sigsq"])
                                      
                   d_mat <- t(comp_spher_cov(cross_dists = u_res_pred,
                                             n = n_obs, n2 = n_pred,
                                             phi   = x$estimate["phi"],
                                             sigsq = x$estimate["sigsq"]))
                   
                   sig_pred <- spher_cov(dists = u_pred,
                                         phi   = x$estimate["phi"],
                                         sigsq = x$estimate["sigsq"])
               })

        sig_y_inv <- chol(chol2inv(sig_y))

        mean_pred <- matrix(rep(x$estimate["alpha"], n_pred),
                            ncol = 1)

        dt_yinv  <- crossprod(d_mat, sig_y_inv)

        sig_pred_y <- sig_pred - (dt_yinv %*% d_mat)

        mean_pred_y <- mean_pred +
            dt_yinv %*% (matrix(x$call_data$var - x$estimate["alpha"],
                                ncol = 1))

        if(any(diag(sig_pred_y) < 0)) {
            warning("Negative variance for at least one predicted region. Taking absolute value.")
            se_pred_y <- sqrt(abs(diag(sig_pred_y)))
        } else {
            se_pred_y <- sqrt(diag(sig_pred_y))
        }
        
        pred_grid <- transform(x$call_data$grid,
                               mu_pred = mean_pred_y,
                               se_pred = se_pred_y)

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
        ## URGENT
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

    if(all(grepl("POINT", sf::st_geometry_type(x)))) {
        coords_pred <- sf::st_coordinates(x)
        pred_grid   <- sf::st_geometry(x)
        warning("The arguments `n_pts`, `type`, and `by_polygon` are ignored when the sf geometry type is POINT")
    } else {
        pred_grid   <- sf::st_sample(x    = sf::st_union(spm_obj$call_data$sf_poly),
                                     size = n_pts, 
                                     by_polygon = FALSE,
                                     type = type)
        coords_pred <- sf::st_coordinates(pred_grid)
    }
    
    u_pred <- as.matrix(stats::dist(coords_pred))
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
                                     sigsq = spm_obj$estimate["sigsq"],
                                     kappa = spm_obj$kappa)


               d_mat <- t(comp_mat_cov(cross_dists = u_res_pred,
                                       n = n_obs, n2 = n_pred,
                                       phi   = spm_obj$estimate["phi"],
                                       sigsq = spm_obj$estimate["sigsq"],
                                       kappa = spm_obj$kappa))
               
               sig_pred <- mat_cov(dists = u_pred,
                                   phi   = spm_obj$estimate["phi"],
                                   sigsq = spm_obj$estimate["sigsq"],
                                   kappa = spm_obj$kappa)
           },
           "pexp" = {
               if(is.null(spm_obj$kappa))
                   spm_obj$kappa <- 1

               sig_y <- comp_pexp_cov(spm_obj$call_data$dists,
                                      n = n_obs, n2 = n_obs,
                                      phi   = spm_obj$estimate["phi"],
                                      sigsq = spm_obj$estimate["sigsq"],
                                      kappa = spm_obj$kappa)
               
               d_mat <- t(comp_pexp_cov(cross_dists = u_res_pred,
                                        n = n_obs, n2 = n_pred,
                                        phi   = spm_obj$estimate["phi"],
                                        sigsq = spm_obj$estimate["sigsq"],
                                        kappa = spm_obj$kappa))
               
               sig_pred <- pexp_cov(dists = u_pred,
                                    phi   = spm_obj$estimate["phi"],
                                    sigsq = spm_obj$estimate["sigsq"],
                                    kappa = spm_obj$kappa)
           },
           "gaussian" = {
               sig_y <- comp_gauss_cov(spm_obj$call_data$dists,
                                       n = n_obs, n2 = n_obs,
                                       phi   = spm_obj$estimate["phi"],
                                       sigsq = spm_obj$estimate["sigsq"])
               
               d_mat <- t(comp_gauss_cov(cross_dists = u_res_pred,
                                         n = n_obs, n2 = n_pred,
                                         phi   = spm_obj$estimate["phi"],
                                         sigsq = spm_obj$estimate["sigsq"]))
               
               sig_pred <- gauss_cov(dists = u_pred,
                                     phi   = spm_obj$estimate["phi"],
                                     sigsq = spm_obj$estimate["sigsq"])
           },
           "spherical" = {
               sig_y <- comp_spher_cov(spm_obj$call_data$dists, 
                                       n = n_obs, n2 = n_obs,
                                       phi   = spm_obj$estimate["phi"],
                                       sigsq = spm_obj$estimate["sigsq"])
               
               d_mat <- t(comp_spher_cov(cross_dists = u_res_pred,
                                         n = n_obs, n2 = n_pred,
                                         phi   = spm_obj$estimate["phi"],
                                         sigsq = spm_obj$estimate["sigsq"]))
               
               sig_pred <- spher_cov(dists = u_pred,
                                     phi   = spm_obj$estimate["phi"],
                                     sigsq = spm_obj$estimate["sigsq"])
           })

    if(all(grepl("POINT", sf::st_geometry_type(x))) & .aggregate) {
        warning("If you want to make predictions only for a set of locations, it does not make sense to use `.aggregate`.")
    }
    
    sig_y_inv <- chol(chol2inv(sig_y))
    
    mean_y <- matrix(rep(spm_obj$estimate["alpha"], n_obs), ncol = 1)
            mean_pred <- matrix(rep(spm_obj$estimate["alpha"], n_pred),
                            ncol = 1)
    
    dt_yinv  <- crossprod(d_mat, sig_y_inv)
    
    sig_pred_y <- sig_pred - (dt_yinv %*% d_mat)
    
    mean_pred_y <- mean_pred +
        dt_yinv %*% (matrix(spm_obj$call_data$var - spm_obj$estimate["alpha"],
                            ncol = 1))

    if(any(diag(sig_pred_y) < 0)) {
        warning("Negative variance for at least one predicted region. Taking absolute value.")
        se_pred_y <- sqrt(abs(diag(sig_pred_y)))
    } else {
        se_pred_y <- sqrt(diag(sig_pred_y))
    }
    
    pred_grid <- transform(sf::st_sf(pred_grid),
                           mu_pred = mean_pred_y,
                           se_pred = se_pred_y)
    
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
