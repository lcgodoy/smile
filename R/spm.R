##' Transforming a \code{sf} into a \code{spm} object (Internal use)
##'
##' @name sf_to_spm
##' 
##' @title single \code{sf} to \code{spm}
##' 
##' @param sf_obj a \code{sf} object s.t. its geometries are polygons.
##' @param n_pts a \code{numeric} scalar representing the number of points to
##'     create a grid in the study region on which the polygons in \code{sf_obj}
##'     is observed. Alternatively, it can be a vector of the same length as
##'     \code{nrow(sf_obj)}. In this case, it generates the given number of
##'     points for each polygon in \code{sf_obj}.
##' @param type a \code{character} indicating the type of grid to be
##'     generated. The options are \code{c("random", "regular",
##'     "hexagonal")}. For more details, see \code{st_sample} in the \code{sf}
##'     package.
##' @param by_polygon a \code{logical} indicating wheter we should generate
##'     \code{n_pts} by polygon or for the \code{n_pts} for the whole study
##'     region.
##' @param poly_ids a \code{character} vector informing the name of the variable
##'     in \code{sf_obj} that represents the polygons unique identifiers. In
##'     case this is not informed, we assume the id of the polygons are given by
##'     their row numbers.
##' @param var_ids a scalar or vector of type \code{character} indicating the
##'     (numerical) variables that are going to be analyzed.
##' @return a \code{list}.
##' @export
single_sf_to_spm <- function(sf_obj,
                             n_pts, type = "regular",
                             by_polygon = FALSE,
                             poly_ids = NULL,
                             var_ids  = NULL) {
    stopifnot(inherits(sf_obj, "sf"))
    stopifnot(all(grepl("POLYGON", sf::st_geometry_type(sf_obj))))
    stopifnot(! is.null(var_ids) )
    stopifnot(inherits(sf_obj, "sf"))
    stopifnot(length(n_pts) == 1 | length(n_pts) == nrow(sf_obj))
    stopifnot(type %in% c("random", "regular", "hexagonal"))
    
    if(is.null(poly_ids) | ( ! poly_ids %in% names(sf_obj))) {
        sf_obj <- transform(sf_obj,
                            poly_ids = seq_len(nrow(sf_obj)))
        poly_ids <- "poly_ids"
    }
    
    if(by_polygon & length(n_pts) == 1) {
        out_grid <-
            lapply(sf::st_cast(sf::st_geometry(sf_obj),
                               "POLYGON"),
                   function(x, .sz, .tp) {
                       sf::st_as_sf(
                               sf::st_sample(x    = x,
                                             size = .sz, 
                                             type = .tp)
                           )
                   },
                   .sz = n_pts, 
                   .tp = type)
    } else if(length(n_pts) > 1) {
        out_grid <-
            Map(function(x, .sz, .tp) {
                sf::st_as_sf(
                        sf::st_sample(x    = x,
                                      size = .sz, 
                                      type = .tp)
                    )
            },
            x = sf::st_cast(sf::st_geometry(sf_obj),
                            "POLYGON"),
            .sz = n_pts, 
            .tp = type)
    } else {
        out_grid <-
            sf::st_sample(x = sf::st_cast(sf::st_geometry(sf_obj),
                                          "POLYGON"),
                          size = n_pts, 
                          type = type,
                          by_polygon = by_polygon)
    }

    if(by_polygon) {
        out_grid <- sf::st_set_crs(do.call("rbind", out_grid),
                                   sf::st_crs(sf_obj))
        
    }
    
    out_grid_pt <-
        sf::st_join(x = sf::st_sf(out_grid),
                    y = sf_obj[poly_ids],
                    join = sf::st_within)
    if(any(! sf_obj[[poly_ids]] %in% out_grid_pt[[poly_ids]])) {
        empty_polys <- which(! sf_obj[[poly_ids]] %in% out_grid_pt[[poly_ids]])
 
        out_grid_aux <-
            sf::st_centroid(x = sf_obj[empty_polys, poly_ids])
        out_grid_pt <- rbind(out_grid_pt, sf::st_sf(out_grid_aux))
    }
    
    out_grid_pt <- out_grid_pt[order(out_grid_pt[[poly_ids]]), ]

    out_grid_pt <- transform(out_grid_pt,
                             x = sf::st_coordinates(out_grid_pt)[, 1],
                             y = sf::st_coordinates(out_grid_pt)[, 2])

    out_dists <- dist_from_grids(out_grid_pt, poly_ids)
    
    if(length(var_ids) == 1) {
        out_var <- sf_obj[[var_ids]]
    } else {
        stop("var_ids must be a single variable")
    }
    
    npix <- as.numeric( sf::st_area(sf_obj) )

    output <- list(
        var     = out_var,
        dists   = out_dists,
        ids_var = poly_ids,
        grid    = out_grid_pt[poly_ids],
        npix    = npix,
        sf_poly = sf::st_geometry(sf_obj)
    )
    
    class(output) <- append(class(output), "spm")
    return(output)
}

##' @name sf_to_spm
##' @export
sf_to_spm <- single_sf_to_spm

print.spm <- function(x, ...) {    
    cat(sprintf("\n Number of variables to be analyzed: %d \n", NCOL(x$var)))
    cat(sprintf("\n ID variable: %ds\n", x$poly_ids))
    if(NCOL(x$var) == 1) {
        cat("\n Variable summary: ")
        summary(x$var)
    } else {
        apply(x$var, 2, summary)
    }
    cat(sprintf("\n Number of polygons in the partition: %d \n", nrow(x$sf_poly)))
    cat(sprintf("\n Number of points in the grid over: %d \n", nrow(x$grid)))
}


## ##' @export
## summary.spm <- function(object, ...) {
##     print.spm(object, ...)
## }
