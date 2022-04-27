##' @name predict_spm
##' @export
predict_spm <- function(x, ...) UseMethod("predict_spm")

##' @name predict_spm
##' @export
predict_spm.spm_fit <- function(x, .aggregate = TRUE, ...) {
    n_obs <- NROW(x$call_data$var)
    ## ids_betas <- which(grepl("^beta", names(x$estimate)))

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
               if(is.null(x$nu))
                   x$nu <- .5

               sig_y <- comp_mat_cov(x$call_data$dists,
                                     n = n_obs, n2 = n_obs,
                                     phi   = x$estimate["phi"],
                                     ## sigsq = x$estimate["sigsq"],
                                     sigsq = 1,
                                     nu = x$nu)
               d_mat <- comp_mat_cov(cross_dists = u_res_pred,
                                     n = n_obs, n2 = n_pred,
                                     phi   = x$estimate["phi"],
                                     sigsq = x$estimate["sigsq"],
                                     nu = x$nu)
               sig_pred <- mat_cov(dists = u_pred,
                                   phi   = x$estimate["phi"],
                                   ## sigsq = x$estimate["sigsq"],
                                   sigsq = 1,
                                   nu = x$nu)
           },
           "pexp" = {
               if(is.null(x$nu))
                   x$nu <- 1

               sig_y <- comp_pexp_cov(x$call_data$dists,
                                      n = n_obs, n2 = n_obs,
                                      phi   = x$estimate["phi"],
                                      ## sigsq = x$estimate["sigsq"],
                                      sigsq = 1,
                                      nu = x$nu)
               
               d_mat <- comp_pexp_cov(cross_dists = u_res_pred,
                                      n = n_obs, n2 = n_pred,
                                      phi   = x$estimate["phi"],
                                      sigsq = x$estimate["sigsq"],
                                      nu = x$nu)
               
               sig_pred <- pexp_cov(dists = u_pred,
                                    phi   = x$estimate["phi"],
                                    ## sigsq = x$estimate["sigsq"],
                                    sigsq = 1,
                                    nu = x$nu)
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
               sig_y <- Matrix(
                   comp_spher_cov(x$call_data$dists, 
                                  n = n_obs, n2 = n_obs,
                                  phi   = x$estimate["phi"],
                                  ## sigsq = x$estimate["sigsq"]
                                  sigsq = 1),
                   sparse = TRUE
               )
               
               d_mat <- Matrix(
                   comp_spher_cov(cross_dists = u_res_pred,
                                  n = n_obs, n2 = n_pred,
                                  phi   = x$estimate["phi"],
                                  sigsq = x$estimate["sigsq"]),
                   sparse = TRUE
               )
               
               sig_pred <- Matrix(
                   spher_cov(dists = u_pred,
                             phi   = x$estimate["phi"],
                             ## sigsq = x$estimate["sigsq"]
                             sigsq = 1),
                   sparse = TRUE
               )
           },
           "cs" = {
               sig_y <- Matrix(
                   comp_cs_cov(x$call_data$dists, 
                               n = n_obs, n2 = n_obs,
                               phi   = x$estimate["phi"],
                               ## sigsq = x$estimate["sigsq"]
                               sigsq = 1),
                   sparse = TRUE
               )
               
               d_mat <- Matrix(
                   comp_cs_cov(cross_dists = u_res_pred,
                               n = n_obs, n2 = n_pred,
                               phi   = x$estimate["phi"],
                               sigsq = x$estimate["sigsq"]),
                   sparse = TRUE
               )
               
               sig_pred <- Matrix(
                   cs_cov(dists = u_pred,
                          phi   = x$estimate["phi"],
                          ## sigsq = x$estimate["sigsq"]
                          sigsq = 1),
                   sparse = TRUE
               )
           },
           "gw" = {
               sig_y <- Matrix(
                   comp_gw_cov(x$call_data$dists, 
                               n = n_obs, n2 = n_obs,
                               phi   = x$estimate["phi"],
                               ## sigsq = x$estimate["sigsq"]
                               sigsq = 1,
                               kappa = x$gw_pars[1],
                               mu    = x$gw_pars[2]),
                   sparse = TRUE
               )
               
               d_mat <- Matrix(
                   comp_gw_cov(cross_dists = u_res_pred,
                               n = n_obs, n2 = n_pred,
                               phi   = x$estimate["phi"],
                               sigsq = x$estimate["sigsq"],
                               kappa = x$gw_pars[1],
                               mu    = x$gw_pars[2]),
                   sparse = TRUE
               )
               
               sig_pred <- Matrix(
                   gw_cov(dists = u_pred,
                          phi   = x$estimate["phi"],
                          ## sigsq = x$estimate["sigsq"]
                          sigsq = 1,
                          kappa = x$gw_pars[1],
                          mu    = x$gw_pars[2]),
                   sparse = TRUE
               )
           },
           "tapmat" = {
               sig_y <- Matrix(
                   comp_tapmat_cov(x$call_data$dists, 
                                   n = n_obs, n2 = n_obs,
                                   phi   = x$estimate["phi"],
                                   ## sigsq = x$estimate["sigsq"]
                                   sigsq = 1,
                                   nu = x$nu,
                                   theta = x$taper_rg),
                   sparse = TRUE
               )
               
               d_mat <- Matrix(
                   comp_tapmat_cov(cross_dists = u_res_pred,
                                   n = n_obs, n2 = n_pred,
                                   phi   = x$estimate["phi"],
                                   sigsq = x$estimate["sigsq"],
                                   nu = x$nu,
                                   theta = x$taper_rg),
                   sparse = TRUE
               )
               
               sig_pred <- Matrix(
                   tapmat_cov(dists = u_pred,
                              phi   = x$estimate["phi"],
                              ## sigsq = x$estimate["sigsq"]
                              sigsq = 1,
                              nu = x$nu,
                              theta = x$taper_rg),
                   sparse = TRUE
               )
           })

    if(length(x$estimate) > 3) {
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

    ## mean_y <- x$call_data$X %*% x$estimate[ids_betas]
    ## mean_pred <- x$call_data$X0 %*% x$estimate[ids_betas]

    sig_y_inv <- chol2inv(chol(sig_y))
    dt_yinv  <- crossprod(d_mat, sig_y_inv)
    sig_pred_y <- sig_pred - (dt_yinv %*% d_mat)

    mean_pred_y <- x$estimate["mu"] +
        dt_yinv %*% (x$call_data$var - x$estimate["mu"])

    if(any(diag(sig_pred_y) < 0)) {
        warning("Negative variance for at least one predicted region. Taking absolute value.")
        var_pred_y <- abs(diag(sig_pred_y))
    } else {
        var_pred_y <- diag(sig_pred_y)
    }
    
    pred_grid <- transform(x$call_data$grid,
                           mu_pred = as.numeric(mean_pred_y),
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
##' @param spm_obj an object of either class \code{spm_fit} or \code{mspm_fit}
##' @param .aggregate \code{logical}. Should the predictions be aggregated? In
##'     case the input is only a "fit" object, the aggregation is made over the
##'     polygons on which the original data was observed. In case the input
##'     \code{x} is composed by \code{sf POLYGONS}, the aggregation is made over
##'     this new partition of the study region.
##' @param n_pts a \code{numeric} scalar standing for number of points to form a
##'     grid over the whole region to make the predictions
##' @param type \code{character} type of grid to be generated. See
##'     \code{st_sample} in the package \code{sf}.
##' @param outer_poly (object) \code{sf geometry} storing the "outer map" we
##'     want to compute the predictions in.
##' @param id_var if \code{x} is a set of \code{POLYGONS} (areal data) instead
##'     of a set of points, the \code{id_var} is the name (or index) of the
##'     unique identifier associated to each polygon.
##' @param ... additional parameters
##'
##' @return a \code{list} of size 4 belonging to the class \code{spm_pred}. This
##'     list contains the predicted values and the mean and covariance matrix
##'     associated with the conditional distribution used to compute the
##'     predictions.
##'
##' @examples
##'
##' data(liv_lsoa) ## loading the LSOA data
##' data(liv_msoa) ## loading the MSOA data
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
##' pred_lsoa <- predict_spm(x = liv_lsoa, spm_obj = fit_msoa, id_var = "lsoa11cd")
##' 
##' @export
predict_spm.sf <- function(x, spm_obj, 
                           n_pts, type,
                           outer_poly = NULL,
                           id_var, ...) {
    if(sf::st_crs(x) != sf::st_crs(spm_obj$call_data$sf_poly)) {
        warning("`x` and the data on which the model was ajdusted are not in the same CRS. Reprojecting `x`")
        x <- sf::st_transform(x, sf::st_crs(spm_obj$call_data$sf_poly))
    }

    n_obs <- NROW(spm_obj$call_data$var)

    stopifnot(all(grepl("(POLYGON|POINT)", sf::st_geometry_type(x))))

    if( all(grepl("POINT", sf::st_geometry_type(x))) ) {
        coords_pred <- sf::st_coordinates(x)
        pred_grid   <- sf::st_geometry(x)
        u_pred      <- distmat(coords_pred)
        n_pred      <- nrow(coords_pred)  # number of locations to make
                                          # predictions
        if( ! missing(n_pts) | ! missing(type) )
            warning("The arguments 'n_pts' and 'type' are ignored when the sf geometry type is POINT.")
    } else {
        if( missing(n_pts) | missing(type) ) {
            pred_grid <- spm_obj$call_data$grid["geometry"]
            pred_grid <- sf::st_join(pred_grid, x, join = sf::st_within)
        } else if(length(n_pts) == NROW(x)) {
            pred_grid <-
                Map(function(geom, .sz, .tp) {
                    sf::st_sample(x = geom,
                                  size = .sz,
                                  type = .tp)
                },
                geom = sf::st_geometry(x),
                .sz  = n_pts,
                .tp  = type)
            pred_grid <- sf::st_set_crs(do.call("c", pred_grid),
                                        sf::st_crs(x))
            pred_grid <- sf::st_join(sf::st_as_sf(pred_grid), x,
                                     join = sf::st_within)
        } else {
            if( is.null(outer_poly) ) {
                outer_poly <- sf::st_union(spm_obj$call_data$sf_poly)
            }
            pred_grid <- sf::st_sample(x    = outer_poly,
                                       size = n_pts, 
                                       by_polygon = FALSE,
                                       type = type)
            pred_grid <- sf::st_join(sf::st_as_sf(pred_grid), x,
                                     join = sf::st_within)
        }
        
        if(any(! x[[id_var]] %in% unique(pred_grid[[id_var]]))) {
            empty_polys <- which(! x[[id_var]] %in% pred_grid[[id_var]])
            grid_aux <-
                sf::st_centroid(x = x[empty_polys, id_var])
            pred_grid <- rbind(pred_grid[, id_var], sf::st_sf(grid_aux))
        }
        pred_grid[["x"]] <- sf::st_coordinates(pred_grid[["geometry"]])[,1]
        pred_grid[["y"]] <- sf::st_coordinates(pred_grid[["geometry"]])[,2]
        ## pred_grid <- pred_grid[, c(3, 1, 2)]
        pred_grid <- pred_grid[order(pred_grid[[id_var]]),]
        rownames(pred_grid) <- NULL
        u_pred <- dist_from_grids(pred_grid, id_var)
        n_pred <- nrow(x)
        coords_pred <- sf::st_coordinates(pred_grid)
    }

    if(all(grepl("POINT", sf::st_geometry_type(x)))) {
        u_res_pred <- pred_cdist(
            get_grid_list(
                x_to_list = spm_obj$call_data$grid,
                by = spm_obj$call_data$ids_var),
            coords_pred
        )
    } else {
        u_res_pred <- mult_dists(mat_list1 =
                                     get_grid_list(
                                         x_to_list = spm_obj$call_data$grid,
                                         by = spm_obj$call_data$ids_var
                                     ),
                                 mat_list2 =
                                     get_grid_list(
                                         x_to_list = pred_grid,
                                         by = id_var
                                     ),
                                 return_single = FALSE)
    }
    
    switch(spm_obj$model,
           "matern" = {
               if(is.null(spm_obj$nu))
                   spm_obj$nu <- .5

               sig_y <- comp_mat_cov(spm_obj$call_data$dists,
                                     n = n_obs, n2 = n_obs,
                                     phi   = spm_obj$estimate["phi"],
                                     ## sigsq = spm_obj$estimate["sigsq"],
                                     sigsq = 1,
                                     nu    = spm_obj$nu)
               d_mat <- comp_mat_cov(cross_dists = u_res_pred,
                                     n = n_obs, n2 = n_pred,
                                     phi   = spm_obj$estimate["phi"],
                                     sigsq = spm_obj$estimate["sigsq"],
                                     nu    = spm_obj$nu)
               if(all(grepl("POINT", sf::st_geometry_type(x)))) {
                   sig_pred <-
                       mat_cov(dists = u_pred,
                               phi   = spm_obj$estimate["phi"],
                               ## sigsq = spm_obj$estimate["sigsq"],
                               sigsq = 1,
                               nu   = spm_obj$nu)
               } else {
                   sig_pred <-
                       comp_mat_cov(cross_dists = u_pred,
                                    n = n_pred, n2 = n_pred,
                                    phi   = spm_obj$estimate["phi"],
                                    sigsq = 1,
                                    nu    = spm_obj$nu)
               }
           },
           "pexp" = {
               if(is.null(spm_obj$nu))
                   spm_obj$nu <- 1

               sig_y <- comp_pexp_cov(spm_obj$call_data$dists,
                                      n = n_obs, n2 = n_obs,
                                      phi   = spm_obj$estimate["phi"],
                                      sigsq = 1,
                                      nu    = spm_obj$nu)
               d_mat <-
                   comp_pexp_cov(cross_dists = u_res_pred,
                                 n = n_obs, n2 = n_pred,
                                 phi   = spm_obj$estimate["phi"],
                                 sigsq = spm_obj$estimate["sigsq"],
                                 nu    = spm_obj$nu)
               if(all(grepl("POINT", sf::st_geometry_type(x)))) {
                   sig_pred <-
                       pexp_cov(dists = u_pred,
                                phi   = spm_obj$estimate["phi"],
                                sigsq = 1,
                                nu    = spm_obj$nu)
               } else {
                   sig_pred <-
                       comp_pexp_cov(cross_dists = u_pred,
                                     n = n_pred, n2 = n_pred,
                                     phi   = spm_obj$estimate["phi"],
                                     sigsq = 1,
                                     nu    = spm_obj$nu)
               }
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
                                       sigsq = spm_obj$estimate["sigsq"])
               if(all(grepl("POINT", sf::st_geometry_type(x)))) {
                   sig_pred <-
                       gauss_cov(dists = u_pred,
                                 phi   = spm_obj$estimate["phi"],
                                 ## sigsq = spm_obj$estimate["sigsq"]
                                 sigsq = 1)
               } else {
                   sig_pred <-
                       comp_gauss_cov(cross_dists = u_pred,
                                      n = n_pred, n2 = n_pred,
                                      phi   = spm_obj$estimate["phi"],
                                      ## sigsq = spm_obj$estimate["sigsq"],
                                      sigsq = 1)
               }
           },
           "spherical" = {
               sig_y <- Matrix(
                   comp_spher_cov(spm_obj$call_data$dists, 
                                  n = n_obs, n2 = n_obs,
                                  phi   = spm_obj$estimate["phi"],
                                  ## sigsq = spm_obj$estimate["sigsq"]
                                  sigsq = 1),
                   sparse = TRUE
               )
               d_mat <- Matrix(
                   comp_spher_cov(cross_dists = u_res_pred,
                                  n = n_obs, n2 = n_pred,
                                  phi   = spm_obj$estimate["phi"],
                                  sigsq = spm_obj$estimate["sigsq"]),
                   sparse = TRUE
               )
               if(all(grepl("POINT", sf::st_geometry_type(x)))) {
                   sig_pred <- Matrix(
                       spher_cov(dists = u_pred,
                                 phi   = spm_obj$estimate["phi"],
                                 ## sigsq = spm_obj$estimate["sigsq"]
                                 sigsq = 1),
                       sparse = TRUE
                   )
               } else {
                   sig_pred <- Matrix(
                       comp_spher_cov(cross_dists = u_pred,
                                      n = n_pred, n2 = n_pred,
                                      phi   = spm_obj$estimate["phi"],
                                      sigsq = spm_obj$estimate["sigsq"]),
                       sparse = TRUE
                   )
               }
           },
           "cs" = {
               sig_y <- Matrix(
                   comp_cs_cov(spm_obj$call_data$dists, 
                               n = n_obs, n2 = n_obs,
                               phi   = spm_obj$estimate["phi"],
                               ## sigsq = spm_obj$estimate["sigsq"]
                               sigsq = 1),
                   sparse = TRUE
               )
               d_mat <- Matrix(
                   comp_cs_cov(cross_dists = u_res_pred,
                               n = n_obs, n2 = n_pred,
                               phi   = spm_obj$estimate["phi"],
                               sigsq = spm_obj$estimate["sigsq"]),
                   sparse = TRUE
               )
               if(all(grepl("POINT", sf::st_geometry_type(x)))) {
                   sig_pred <- Matrix(
                       cs_cov(dists = u_pred,
                              phi   = spm_obj$estimate["phi"],
                              ## sigsq = spm_obj$estimate["sigsq"]
                              sigsq = 1),
                       sparse = TRUE
                   )
               } else {
                   sig_pred <- Matrix(
                       comp_cs_cov(cross_dists = u_pred,
                                   n = n_pred, n2 = n_pred,
                                   phi   = spm_obj$estimate["phi"],
                                   ## sigsq = spm_obj$estimate["sigsq"],
                                   sigsq = 1),
                       sparse = TRUE
                   )
               }
           },
           "gw" = {
               sig_y <- Matrix(
                   comp_gw_cov(spm_obj$call_data$dists, 
                               n = n_obs, n2 = n_obs,
                               phi   = spm_obj$estimate["phi"],
                               sigsq = 1,
                               kappa = spm_obj$gw_pars[1],
                               mu    = spm_obj$gw_pars[2]),
                   sparse = TRUE
               )
               d_mat <- Matrix(
                   comp_gw_cov(cross_dists = u_res_pred,
                               n = n_obs, n2 = n_pred,
                               phi   = spm_obj$estimate["phi"],
                               sigsq = spm_obj$estimate["sigsq"],
                               kappa = spm_obj$gw_pars[1],
                               mu    = spm_obj$gw_pars[2]),
                   sparse = TRUE
               )
               if(all(grepl("POINT", sf::st_geometry_type(x)))) {
                   sig_pred <- Matrix(
                       gw_cov(dists = u_pred,
                              phi   = spm_obj$estimate["phi"],
                              sigsq = 1,
                              kappa = spm_obj$gw_pars[1],
                              mu    = spm_obj$gw_pars[2]),
                       sparse = TRUE
                   )
               } else {
                   sig_pred <- Matrix(
                       comp_gw_cov(cross_dists = u_pred,
                                   n = n_pred, n2 = n_pred,
                                   phi   = spm_obj$estimate["phi"],
                                   ## sigsq = spm_obj$estimate["sigsq"],
                                   sigsq = 1,
                                   kappa = spm_obj$gw_pars[1],
                                   mu    = spm_obj$gw_pars[2]),
                       sparse = TRUE
                   )
               }
           },
           "tapmat" = {
               sig_y <- Matrix(
                   comp_tapmat_cov(spm_obj$call_data$dists, 
                                   n = n_obs, n2 = n_obs,
                                   phi   = spm_obj$estimate["phi"],
                                   sigsq = 1,
                                   nu    = spm_obj$nu,
                                   theta = spm_obj$taper_rg),
                   sparse = TRUE
               )
               d_mat <- Matrix(
                   comp_tapmat_cov(cross_dists = u_res_pred,
                                   n = n_obs, n2 = n_pred,
                                   phi   = spm_obj$estimate["phi"],
                                   sigsq = spm_obj$estimate["sigsq"],
                                   nu    = spm_obj$nu,
                                   theta = spm_obj$taper_rg),
                   sparse = TRUE
               )
               if(all(grepl("POINT", sf::st_geometry_type(x)))) {
                   sig_pred <- Matrix(
                       tapmat_cov(dists = u_pred,
                                  phi   = spm_obj$estimate["phi"],
                                  sigsq = 1,
                                  nu = spm_obj$nu,
                                  theta = spm_obj$taper_rg),
                       sparse = TRUE
                   )
               } else {
                   sig_pred <- Matrix(
                       comp_tapmat_cov(cross_dists = u_pred,
                                       phi   = spm_obj$estimate["phi"],
                                       sigsq = 1,
                                       nu = spm_obj$nu,
                                       theta = spm_obj$taper_rg),
                       sparse = TRUE
                   )
               }
           })

    if(length(spm_obj$estimate) > 3) {
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
    mean_y <- matrix(rep(spm_obj$estimate["mu"], n_obs),
                     ncol = 1)
    mean_pred <- matrix(rep(spm_obj$estimate["mu"], n_pred),
                        ncol = 1)
    dt_yinv <- crossprod(d_mat, sig_y_inv)
    sig_pred_y <- sig_pred - (dt_yinv %*% d_mat)
    mean_pred_y <- mean_pred +
        dt_yinv %*% (spm_obj$call_data$var - mean_y)

    if(any(diag(sig_pred_y) < 0)) {
        warning("Negative variance for at least one predicted region. Taking absolute value.")
        var_pred_y <- abs(diag(sig_pred_y))
    } else {
        var_pred_y <- diag(sig_pred_y)
    }
    
    if(all(grepl("POINT", sf::st_geometry_type(x)))) {
        if(inherits(pred_grid, "sfc"))
            pred_grid <- transform(sf::st_sf(geometry = pred_grid),
                                   mu_pred = as.numeric(mean_pred_y),
                                   se_pred = var_pred_y)
        else 
            pred_grid <- transform(pred_grid,
                                   mu_pred = as.numeric(mean_pred_y),
                                   se_pred = var_pred_y)
        pred_grid[["se_pred"]] <- sqrt(pred_grid[["se_pred"]])
        output <-
            list(
                mu_pred   = mean_pred_y,
                sig_pred  = sig_pred_y,
                pred_grid = pred_grid,
                pred_agg  = NA
            )
    }
    else {
        pred_grid <- transform(x[order(x[[id_var]]), ],
                               mu_pred = as.numeric(mean_pred_y),
                               se_pred = sqrt(var_pred_y))
        output <-
            list(
                mu_pred   = mean_pred_y,
                sig_pred  = sig_pred_y,
                pred_grid = NA,
                pred_agg  = pred_grid
            )
    }
    class(output) <- append(class(output), "spm_pred")
    return(output)
}

##' @title Internal use only
##' @param x internal use
##' @param by internal use
##' @param FUN internal use
##' @param ... internal use
##' @param do_union internal use
##' @param simplify internal use
##' @param join internal use
##' @description {aggregate.sf} taken from \code{sf}.
##' @keywords internal
aggregate_aux <- utils::getFromNamespace("aggregate.sf", "sf")
