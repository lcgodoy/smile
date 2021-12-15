##' @title Voronoi Tesselation inside a polygon
##'
##' @description voronoi tesselation of a given a set of points inside a polygon.
##' This is an internal use function.
##'
##' @param points_sf `sf data frame` containing the points' coordinates
##' @param poly_sf a `sf` polygon
##'
##' @return a `sf` polygon
##' @keywords internal
vor_build <- function(points_sf, poly_sf) {
    if(nrow(points_sf) > 1) {
        voronoi <- do.call(c, sf::st_geometry(points_sf))
        voronoi <- sf::st_voronoi(voronoi, envelope = poly_sf)
        voronoi <- sf::st_collection_extract(voronoi, type = 'POLYGON')
        voronoi <- sf::st_set_crs(voronoi, sf::st_crs(points_sf))
        voronoi <- sf::st_intersection(voronoi, poly_sf)
    } else {
        voronoi <- poly_sf
    }
    voronoi <- sf::st_as_sf(voronoi)

    output <- sf::st_join(x    = voronoi,
                          y    = points_sf,
                          join = sf::st_nearest_feature)

    return(output)
}

##' @title Voronoi Data Linkage
##'
##' @description Reminder, have to create an example.
##'
##' @param coords_sf \code{sf} POINT target dataset.
##' @param areal_sf \code{sf} POLYGON source dataset.
##' @param vars a \code{character} representing the variables (observed at the
##'     source - polygon) to be estimated at the target data.
##' @param buff scalar `numeric`. Mostly for internal use.
##'
##'
##' @return a `sf` object for the `coords_sf` spatial data set.
##' @export
vdl <- function(coords_sf, areal_sf,
                vars,
                buff) {
    if(!all(inherits(coords_sf, "sf"), inherits(areal_sf, "sf")))
        stop("coords_sf and areal_sf must be sf objects.")

    if(missing(vars))
        stop("input at least one pop_vars.")
    ## verifying crs
    crs_coords <- sf::st_crs(coords_sf)
    crs_areal  <- sf::st_crs(areal_sf)
    if(crs_coords != crs_areal) {
        coords_sf <- sf::st_transform(coords_sf, crs_areal)
        warning("Making coordinates CRS compatible to areal data CRS.")
    }

    ## if(is.list(list_vars)) {
    ##     all_vars <- c(names(list_vars),
    ##                   unlist(list_vars, use.names = FALSE))
    ## } else {
    ##     all_vars <- list_vars
    ## }

    if(!all(vars %in% names(areal_sf)))
        stop("all variables inputed in pop_vars and avg_vars must be in areal_sf.")

    sf_border <-
        tryCatch(
            expr = sf::st_union(sf::st_geometry(areal_sf)),
            error = function(e) {
                sf::st_union(
                        sf::st_buffer(
                                sf::st_geometry(areal_sf),
                                if(missing(buff)) .5 else buff
                            )
                    )
            }
        )

    vor_sf <- vor_build(points_sf = coords_sf, poly_sf = sf_border)    
    
    output <- ai(source = areal_sf, target = vor_sf,
                 vars = vars)

    return(output)
}

##' @title Voronoi Data Linkage - Single variable and variance
##'
##' @description Reminder, have to create an example.
##'
##' @param coords_sf \code{sf} POINT target dataset.
##' @param areal_sf \code{sf} POLYGON source dataset.
##' @param res_var a \code{character} - the name of the variable in the
##'     \code{areal_sf} to be estimated in the \code{coords_sf}.
##' @param variance a \code{character} - the name of the variable varinace in
##'     the \code{areal_sf} to be estimated in the \code{coords_sf}.
##' @param var_method a \code{character} representing the method to approximate
##'     the variance of the AI estimates. Possible values are "CS"
##'     (Cauchy-Schwartz) or "MI" (Moran's I).
##' @param buff scalar `numeric`. Mostly for internal use.
##'
##'
##' @return a `sf` object, contaning the `id_coords` variable and the
##' `list_vars` for the `coords_sf` spatial data set.
##' @export
vdl_var <- function(coords_sf, areal_sf,
                    res_var,
                    variance,
                    var_method = "CS",
                    buff) {
    if(!all(inherits(coords_sf, "sf"), inherits(areal_sf, "sf")))
        stop("coords_sf and areal_sf must be sf objects.")

    ## verifying crs
    crs_coords <- sf::st_crs(coords_sf)
    crs_areal  <- sf::st_crs(areal_sf)
    if(crs_coords != crs_areal) {
        coords_sf <- sf::st_transform(coords_sf, crs_areal)
        warning("Making coordinates CRS compatible to areal data CRS.")
    }

    if(!all(c(res_var, variance) %in% names(areal_sf)))
        stop(sprintf("%s and %s are not in areal_sf.", res_var, variance))

    sf_border <-
        tryCatch(
            expr = sf::st_union(sf::st_geometry(areal_sf)),
            error = function(e) {
                sf::st_union(
                        sf::st_buffer(
                                sf::st_geometry(areal_sf),
                                if(missing(buff)) .5 else buff
                            )
                    )
            }
        )
    vor_sf <- vor_build(points_sf = coords_sf, poly_sf = sf_border)
    output <- ai_var(source = areal_sf, target = vor_sf,
                     vars = res_var, vars_var = variance,
                     var_method = var_method)
    return(output)
}

##' Transform method for `sf` objects
##'
##' Internal usage.
##'
##' @param _data a `sf` object
##' @param ... additional options
##'
##' @return an `sf` object.
##' @keywords internal
transform.sf <- function (`_data`, ...) {
    e <- eval(substitute(list(...)), `_data`, parent.frame())
    tags <- names(e)
    inx <- match(tags, names(`_data`))
    matched <- !is.na(inx)
    crs <- sf::st_crs(`_data`)
    if (any(matched)) {
        `_data`[inx[matched]] <- e[matched]
        `_data` <- sf::st_sf(`_data`,
                             crs = crs)
    }
    if (!all(matched))
        sf::st_sf(do.call("data.frame", c(list(`_data`), e[!matched])),
                  crs = crs)
    else `_data`
}
