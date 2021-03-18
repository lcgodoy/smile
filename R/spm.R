##' Transforming a \code{sf} into a \code{spm} object (Internal use)
##'
##' @title single \code{sf} to \code{spm}
##' @param sf_obj a \code{sf} object s.t. its geometries are polygons.
##' @param n_pts a \code{numeric} scalar representing the number of points to
##'     create a grid in the study region on which the polygons in \code{sf_obj}
##'     is observed. Alternatively, it can be a vector of the same length as
##'     \code{nrow(sf_obj)}. In this case, it generates the given number of
##'     points for each polygon in \code{sf_obj}.
##' @param type a \code{character} indicating the type of grid to be
##'     generated. The options are \code{c("random", "regular",
##'     "hexagonal")}. For more details, see \code{\link[sf]{st_sample}}.
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
        out_grid <- st_set_crs(do.call("rbind", out_grid),
                               st_crs(sf_obj))
        
    }
    
    out_grid_pt <-
        sf::st_join(x = sf::st_sf(out_grid),
                    y = sf_obj[poly_ids],
                    join = sf::st_within)

    out_grid_pt <- transform(out_grid_pt,
                             x = sf::st_coordinates(out_grid_pt)[, 1],
                             y = sf::st_coordinates(out_grid_pt)[, 2])
    
    out_dists <- dist_from_grids(out_grid_pt, poly_ids)

    if(length(var_ids) == 1) {
        out_var <- sf_obj[[var_ids]]
    } else {
        out_var <-
            do.call(
                "cbind",
                lapply(var_ids, function(id, y) y[[id]])
            )
    }
    
    output <- list(
        var     = out_var,
        dists   = out_dists,
        ids_var = poly_ids,
        grid    = out_grid_pt[poly_ids],
        sf_poly = sf::st_geometry(sf_obj)
    )

    return(output)
}

##' Transforming one (or two) \code{sf} objects into a \code{sspm} (or \code{mspm}) object
##' 
##' @title \code{sf} to \code{spm}
##' @param sf_obj1 a \code{sf} object s.t. its geometries are polygons.
##' @param sf_obj2 a \code{sf} object s.t. its geometries are polygons.
##' @param n_pts a \code{numeric} scalar representing the number of points to
##'     create a grid in the study region on which the polygons in \code{sf_obj}
##'     is observed.  Alternatively, it can be a vector of the same length as
##'     \code{nrow(sf_obj)}.  In this case, it generates the given number of
##'     points for each polygon in \code{sf_obj}. Considering the scenario on
##'     which you also input a \code{sf_obj2}, this argument can be set to be a
##'     list with two positions such that the first position contains the
##'     \code{n_pts} argument for \code{sf_obj1}, and the second position
##'     contains informations about the same parameter but for the
##'     \code{sf_obj2}.
##' @param type a \code{character} indicating the type of grid to be generated.
##'     The options are \code{c("random", "regular", "hexagonal")}. For more
##'     details, see \code{\link[sf]{st_sample}}. Similarly to \code{n_pts}, you
##'     can input a two positions list in here too.
##' @param by_polygon a \code{logical} indicating wheter we should generate
##'     \code{n_pts} by polygon or for the \code{n_pts} for the whole study
##'     region. Similarly to \code{n_pts} (and \code{by_polygon}), you can input
##'     a two positions list in here too.
##' @param poly_ids a \code{list} with two positions, being the first and the
##'     second one a \code{character} vector informing the name of the variable
##'     in \code{sf_obj} and \code{sf_obj2}, respectively, that represents the
##'     polygons unique identifiers. In case this is not informed, we assume the
##'     id of the polygons are given by their row numbers. Note that, if this is
##'     informed, it must be a two positions list.
##' @param var_ids a two positions \code{list} containing in each position
##'     either a scalar or vector of type \code{character}, indicating the
##'     (numerical) variables that are going to be analyzed for \code{sf_obj1}
##'     and \code{sf_obj2}, respectively.
##'
##' @return a \code{sspm} or \code{mspm} object.
##' 
##' @export
sf_to_spm <- function(sf_obj1,
                      sf_obj2,
                      n_pts, type = "regular",
                      by_polygon = FALSE,
                      poly_ids = NULL,
                      var_ids) {
    stopifnot(!missing(sf_obj1))
    stopifnot(!missing(n_pts))
    stopifnot(!missing(var_ids))

    if(missing(sf_obj2)) {
        stopifnot(length(poly_ids) == 1)
        output <-
            single_sf_to_spm(sf_obj1, n_pts, type,
                             by_polygon,
                             poly_ids, var_ids)
        class(output) <- append(class(output), "sspm")
    } else {
        stopifnot(length(var_ids) == 2)
        stopifnot(length(poly_ids) == 2)
        stopifnot(all(lengths(poly_ids) == 1))

        y <-
            single_sf_to_spm(sf_obj1,
                             n_pts = ifelse(length(n_pts) == 1,
                                            n_pts, n_pts[[1]]),
                             type = ifelse(length(type) == 1,
                                           type, type[1]),
                             by_polygon = ifelse(length(by_polygon) == 1,
                                                 by_polygon, by_polygon[1]),
                             poly_ids = poly_ids[[1]],
                             var_ids = var_ids[[1]])

        x <-
            single_sf_to_spm(sf_obj2,
                             n_pts = ifelse(length(n_pts) == 1,
                                            n_pts, n_pts[[2]]),
                             type = ifelse(length(type) == 1,
                                           type, type[2]),
                             by_polygon = ifelse(length(by_polygon) == 1,
                                                 by_polygon, by_polygon[2]),
                             poly_ids = poly_ids[[2]],
                             var_ids = var_ids[[2]])

        output <-
            list(
                y = y,
                x = x,
                cdist = mult_dist_from_grids(y$grid, x$grid, poly_ids)
            )
        
        class(output) <- append(class(output), "mspm")
    }

    return(output)
}

##' @export
print.sspm <- function(x, ...) {
    
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


##' @export
summary.sspm <- function(object, ...) {
    print.sspm(object, ...)
}

##' @export
print.mspm <- function(x, ...) {
    cat(sprintf("\n Number of variables to be analyzed in the partition A: %d \n", NCOL(x$y$var)))
    cat(sprintf("ID for the polygons observed in the partition A: %ds\n", x$y$poly_ids))
    cat(sprintf("Number of polygons in the partition A: %d \n", nrow(x$y$sf_poly)))
    cat(sprintf("Number of points in the grid over A: %d \n \n", nrow(x$y$grid)))
    cat(sprintf("Number of variables to be analyzed in the partition B: %d \n", NCOL(x$x$var)))
    cat(sprintf("ID for the polygons observed in the partition B: %ds\n", x$x$poly_ids))
    cat(sprintf("Number of polygons in the partition B: %d \n", nrow(x$x$sf_poly)))
    cat(sprintf("Number of points in the grid over B: %d \n", nrow(x$x$grid)))
}


##' @export
summary.mspm <- function(object, ...) {
    print.mspm(object, ...)
}
