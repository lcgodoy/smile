##' @title MLEs for fixed V.
##' @details internal use.
##' @param y response variable.
##' @param Vinv inverse of \eqn{V}
##' @keywords internal
est_mle <- function(y, Vinv) {
    n <- length(y)
    ones  <- matrix(1, nrow = n)
    IVinv <- crossprod(ones, Vinv)
    mu    <- as.numeric((IVinv %*% y) / (IVinv %*% ones))
    y     <- matrix(y - mu, ncol = 1)
    sigsq <- as.numeric((crossprod(y, Vinv) %*% y) / n)
    out   <- c("mu" = mu, "sigsq" = sigsq)
    return(out)
}

##' @name aux_mat
get_grid_list <- function(x_to_list, by) {
    x_list <- sf::st_coordinates(x_to_list)
    
    x_list <- split(x = as.data.frame(x_list),
                    f = factor(x_to_list[[by]],
                               levels = unique(x_to_list[[by]])))
    
    return(lapply(x_list, as.matrix))
}

##' @name aux_mat
dist_from_grids <- function(y_grid,  by) {
    out_list <- split(x = sf::st_set_geometry(y_grid, NULL),
                      f = factor(y_grid[[by]],
                                 levels = unique(y_grid[[by]])))
    return(
        single_dists(
            lapply(out_list,
                   function(x) as.matrix(x[ , c("x", "y")]))
        )
    )
}

##' @name aux_mat
mult_dist_from_grids <- function(y_grid, x_grid, by) {
    y_list <- sf::st_coordinates(y_grid)
    
    y_list <- split(x = as.data.frame(y_list),
                    f = factor(y_grid[[by[2]]],
                               levels = unique(y_grid[[by[1]]])))

    x_list <- sf::st_coordinates(x_grid)
    
    x_list <- split(x = as.data.frame(x_list),
                      f = factor(x_grid[[by[1]]],
                                 levels = unique(x_grid[[by[2]]])))
    
    return(mult_dists(lapply(y_list, as.matrix),
                      lapply(x_list, as.matrix),
                      FALSE))
}

##' @title Remove holes from a \code{sfc} POLYGON
##' @description internal use. Taken from
##'     \url{https://github.com/michaeldorman/nngeo/}
##' @param x a \code{sf} or \code{sfc} polygon.
##' @return a \code{sf} or \code{sfc} polygon.
##' @keywords internal
st_remove_holes <- function(x) {
    stopifnot(all(sf::st_is(x, "POLYGON") | sf::st_is(x, "MULTIPOLYGON")))
    geometry_is_polygon <- all(sf::st_is(x, "POLYGON"))
    type_is_sfg <- any(inherits(x, "sfg"))
    type_is_sf <- any(inherits(x, "sf"))
    geom <- sf::st_geometry(x)
    if (type_is_sf) 
        dat <- sf::st_set_geometry(x, NULL)
    for (i in 1:length(geom)) {
        if (sf::st_is(geom[i], "POLYGON")) {
            if (length(geom[i][[1]]) > 1) {
                geom[i] <- sf::st_multipolygon(lapply(geom[i], function(p) p[1]))
            }
        }
        if (sf::st_is(geom[i], "MULTIPOLYGON")) {
            tmp <- sf::st_cast(geom[i], "POLYGON")
            for (j in 1:length(tmp)) {
                if (length(tmp[j][[1]]) > 1) {
                    tmp[j] <- sf::st_multipolygon(lapply(tmp[j], function(p) p[1]))
                }
            }
            geom[i] = sf::st_combine(tmp)
        }
    }
    if (geometry_is_polygon) 
        geom = sf::st_cast(geom, "POLYGON")
    if (type_is_sfg) 
        geom = geom[[1]]
    if (type_is_sf) 
        geom = sf::st_sf(dat, geom)
    return(geom)
}

##' @title Find phi parameter for the Exponential spatial auto-correlation
##'     function
##'
##' @description Function designed to find the phi parameter such that the
##'     correlation between points within a given distance \code{d} is at most a
##'     given value.
##'
##' @param d maximum distance for spatial dependence equal to \code{cut}.
##' @param nu smoothness parameter associated with the Matern cov. function.
##' @param kappa one of the smoothness parameters associated with the Generalized
##'     Wendland covariance function
##' @param mu2 one of the smoothness parameters associated with the Generalized
##'     Wendland covariance function
##' @param range Minimum and maximum distance to be considered. The default is
##'     \code{range = c(1e-04, 1000)}.
##' @param family covariance function family, the options are \code{c("matern",
##'     "gw", "cs", "spher", "pexp", "gaussian")}.
##' @param cut desired spatial correlation at a distance \code{d}, the default is
##'     \code{cut = .05}.
##'
##' @return a \code{numeric} value indicating the range parameter such that the
##'     spatial correlation between two points at distance \code{d} is
##'     \code{cut}.
##' @export
find_phi <- function(d, nu, kappa, mu2, family = "matern",
                     range = c(1e-04, 1000), cut = 0.05) {
    if(family %in% c("matern", "gw", "cs",
                     "spher", "pexp",
                     "gaussian")) {
        out <- stats::uniroot(f = function(x, d, nu, cut) {
            if(family == "matern") {
                out <- single_matern(d, 1, x, nu)
            } else if(family == "gw") {
                out <- single_gw(d, 1, x, kappa, mu2)
            } else if(family == "cs") {
                out <- single_cs(d, 1, x)
            } else if(family == "pexp") {
                out <- single_pexp(d, 1, x, nu)
            } else if(family == "gaussian") {
                out <- single_gauss(d, 1, x)
            } else {
                out <- single_spher(d, 1, x)
            }
            out - cut
        }, interval = range, d = d, nu = nu, cut = cut)
        out$root
    } else {
        - d / log(cut)
    }
}

##' @title Calculate Smallest Eigenvalue for Matern and Power Exponential
##'   Correlation Matrices
##'
##' @description This function computes the smallest eigenvalue of a correlation
##'   matrix derived from the power exponential correlation function. It
##'   evaluates this across a grid of values for the power parameter (`nu`) and
##'   the practical range parameter (`rho`), based on a provided distance
##'   matrix.
##'
##' @details The practical range `rho` is defined here as the distance at which
##'   the correlation is 0.1.
##'
##' @param range_nu A numeric vector of length 1 or 2 especifying the smoothness
##'   parameter.
##' @param range_rho A numeric vector of length 2, specifying the minimum and
##'   maximum values for the practical range parameter `rho`. `rho` must be
##'   positive.
##' @param grid_len An integer specifying the number of points to create for
##'   both `nu` and `rho` sequences. The total number of grid combinations will
##'   be `grid_len^2`. Default is 50.
##' @param dmat A numeric matrix representing the distance matrix between
##'   locations.  The distances should be non-negative.
##'
##' @return A [tibble][dplyr::tibble] with three columns:
##'   \item{rho}{The practical range parameter value.}
##'   \item{nu}{The power parameter value.}
##'   \item{lambda}{The smallest eigenvalue of the power exponential correlation
##'                 matrix corresponding to the `rho` and `nu` pair.}
##'
##' @details The function first creates a grid of `nu` and `rho` parameters. For
##'   each pair of (`rho`, `nu`) in the grid: 1. It calculates the scale
##'   parameter `phi` for the power exponential correlation function, where `phi
##'   = rho / (log(10)^(1/nu))`. This definition implies that the correlation is
##'   0.1 at the distance `rho`.  2. It computes the power exponential
##'   correlation matrix using `smile:::pexp_cov(dists = dmat, sill = 1, range =
##'   phi, smooth = nu)`.  Note the use of an internal function from the `smile`
##'   package.  3. It calculates the eigenvalues of this correlation matrix.
##'   4. The minimum eigenvalue is extracted.  The final output is a tibble
##'   containing all parameter combinations and their corresponding minimum
##'   eigenvalues.
##'
##' @name sev
##' @export
sev_pexp <- function(range_nu, range_rho,
                     grid_len = 50,
                     dmat) {
  rhos <- unique(seq(from = range_rho[1], range_rho[2],
                     length.out = grid_len))
  nus  <- unique(seq(from = range_nu[1], range_nu[2],
                     length.out = grid_len))
  pars_mat <- expand.grid(rho = rhos, nu = nus)
  min_lbd <-
    sapply(seq_len(nrow(pars_mat)),
           \(i, pmts, dst) {
             phi <-
               pmts$rho[i] /
               (log(10) ^ (1 / pmts$nu[i]))
             out <-
               eigen(pexp_cov(dists = dst,
                              1,
                              phi,
                              pmts$nu[i]))
             min(out$values)
           }, pmts = pars_mat,
           dst = dmat)
  transform(as.data.frame(pars_mat), lambda = min_lbd)
}

##' @rdname sev
##' @export
sev_mat <- function(range_nu, range_rho,
                    grid_len = 50,
                    dmat) {
  rhos <- unique(seq(from = range_rho[1], range_rho[2],
                     length.out = grid_len))
  nus  <- unique(seq(from = range_nu[1], range_nu[2],
                     length.out = grid_len))
  pars_mat <- expand.grid(rho = rhos, nu = nus)
  min_lbd <-
    sapply(seq_len(nrow(pars_mat)),
           \(i, pmts, dst) {
             phi <-
               pmts$rho[i] / sqrt(8 * pmts$nu[i])
             out <-
               eigen(mat_cov(dists = dst,
                             1,
                             phi,
                             pmts$nu[i]))
             min(out$values)
           }, pmts = pars_mat,
           dst = dmat)
  transform(as.data.frame(pars_mat), lambda = min_lbd)
}
