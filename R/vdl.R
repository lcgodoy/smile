##' Voronoi Tesselation inside a polygon
##'
##' @description voronoi tesselation of a given a set of points inside a polygon.
##' This is an internal use function.
##'
##' @param points_sf `sf data frame` containing the points' coordinates
##' @param poly_sf a `sf` polygon
##'
##' @return a `sf` polygon
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

##' Voronoi Data Linkage
##'
##' Reminder, have to create an example.
##'
##' @param coords_sf `sf` object representing point coordinates
##' @param areal_sf `sf` object representing the areal data (e.g. spatial polygons)
##' @param id_coords `character` indicating the _id_ variable for `coords_sf`. If
##' not informed, the _row number_ is assumed to by the id.
##' @param list_vars a named `list` representing the variables from `areal_sf`
##' that need to be estimated in `coords_sf`. The `list` position names represent
##' the populational variables, while the contend within each position represent
##' the average variavles.
##' @param show_warning `logical`, indicating wheter a warning for the crs
##' transformation should be printed or not. (default is `TRUE`)
##' @param buff scalar `numeric`. Mostly for internal use.
##'
##' @importFrom data.table .I .SD .N ':='
##'
##' @return a `sf` object, contaning the `id_coords` variable and the
##' `list_vars` for the `coords_sf` spatial data set.
##' @export
##'
vdl <- function(coords_sf, areal_sf,
                id_coords = NULL,
                list_vars = NULL,
                show_warning = TRUE,
                buff) {
    if(!all(inherits(coords_sf, "sf"), inherits(areal_sf, "sf")))
        stop("coords_sf and areal_sf must be sf objects.")

    if(is.null(list_vars))
        stop("input at least one pop_vars.")
    if(is.null(id_coords)) {
        id_coords <- "id_coords"
        coords_sf$id_coords <- seq_len(nrow(coords_sf))
    } else {
        if(! (id_coords %in% names(coords_sf))) {
            coords_sf[[id_coords]] <- seq_len(nrow(coords_sf))
        }
    }

                                        # verifying crs
    crs_coords <- sf::st_crs(coords_sf)
    crs_areal  <- sf::st_crs(areal_sf)
    ## if(!(is.na(crs_coords) | crs_coords$epsg == 3857)) {
    ##   coords_sf <- sf::st_transform(x   = coords_sf,
    ##                                 crs = 3857)
    ##   if(show_warning) {
    ##     message("The coords_sf CRS was changed to EPSG:3857")
    ##   }
    ## }

    ## if(!(is.na(crs_areal) | crs_areal$epsg == 3857)) {
    ##   areal_sf <- sf::st_transform(x   = areal_sf,
    ##                                crs = 3857)
    ##   if(show_warning) {
    ##     message("The areal_sf CRS was changed to EPSG:3857")
    ##   }
    ## }

    if(sf::st_crs(areal_sf) != sf::st_crs(coords_sf))
        stop("Error - different CRSs.")

    if(is.list(list_vars)) {
        all_vars <- c(names(list_vars),
                      unlist(list_vars, use.names = F))
    } else {
        all_vars <- list_vars
    }

    if(!all(all_vars %in% names(areal_sf)))
        stop("all variables inputed in pop_vars and avg_vars must be in areal_sf.")

    areal_sf <- areal_sf[, all_vars]
    areal_sf$geom_2 <- sf::st_geometry(areal_sf) # auxiliar geometry
    ## areal_sf <- transform(areal_sf, geom_2 = sf::st_geometry(areal_sf))

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
    vor_join <- vor_sf
    vor_join <- transform(vor_join, cell_area = as.numeric(sf::st_area(geometry)))
    vor_join <- sf::st_join(
                        x = vor_join,
                        y = areal_sf,
                        join = sf::st_intersects
                    )
    vor_join <- transform(vor_join,
                          inter_area = mapply(function(x, y) {
                              as.numeric(sf::st_area(
                                                 sf::st_intersection(x, y)
                                             ))
                          }, x = geometry, y = geom_2))
    vor_join <- transform(vor_join, prop_inter = inter_area/cell_area)
    vor_join <- sf::st_set_geometry(x = vor_join, value = NULL)
    vor_join <- transform(vor_join, geom_2 = NULL)
    data.table::setDT(vor_join)

    output <- split(x = vor_join, f = vor_join[[id_coords]])

    output <- lapply(output, function(x, l_vars) {
        if(is.list(l_vars)) {
            x <- transform(x, prop_inter = prop_inter/sum(prop_inter, na.rm = TRUE))
            for(i in seq_along(l_vars)) {
                pop_var <- names(l_vars)[i]
                mu_var  <- l_vars[[i]]
                if(! is.null(mu_var) )
                    x[, c(mu_var) := Map(f = function(y, prop) {
                        (get(pop_var)*y*prop)/sum(get(pop_var)*prop, na.rm = TRUE)
                    }, y = .SD, prop = list(prop_inter)), .SDcols = mu_var]
            }
            x[, (names(l_vars)) := Map(function(y, prop) {y*prop}, y = .SD, prop = list(prop_inter)),
              .SDcols = (names(l_vars))]
            num_idx <- names(which(sapply(x, is.numeric), useNames = T))
            x[ , lapply(.SD, sum), .SDcols = num_idx, by = eval(id_coords)]
        } else {
            x <- transform(x, prop_inter = prop_inter/sum(prop_inter, na.rm = TRUE))
            x[, (names(l_vars)) := Map(function(y, prop) {y*prop}, y = .SD, prop = list(prop_inter)),
              .SDcols = (names(l_vars))]
            num_idx <- names(which(sapply(x, is.numeric), useNames = T))
            x[ , lapply(.SD, sum), .SDcols = num_idx, by = eval(id_coords)]
        }
    },
    l_vars = list_vars)

    vor_sf <- vor_sf[c(id_coords, 'geometry')]

    output <- data.table::rbindlist(output)
    output <- merge(x = vor_sf, y = output, by = id_coords)

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
