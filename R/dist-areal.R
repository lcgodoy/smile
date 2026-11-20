## Hausdorff and center-radii (CR) distances between spatial units.
##
## The C++ code is in src/haus-r2.cpp (planar Hausdorff), src/haus-s2.cpp
## (spherical Hausdorff), and src/cr-dist.cpp (CR and Chebyshev centers).

## ---- Internal helpers -------------------------------------------------------

## Stops if `x` and `y` do not have the same CRS.
.check_same_crs <- function(x, y) {
  if (sf::st_crs(x) != sf::st_crs(y)) {
    stop("x and y must have the same CRS.", call. = FALSE)
  }
  invisible(TRUE)
}

## Returns the geometry of `x`, ready for the C++ code. On the sphere, the
## coordinates must be longitude and latitude (EPSG:4326). In the plane, they
## are used as they are. `arg` is the argument name for the warnings.
.prep_geom <- function(x, use_sphere = FALSE, arg = "x") {
  x <- sf::st_geometry(x)
  if (use_sphere) {
    if (is.na(sf::st_crs(x))) {
      warning(arg, " has no CRS. Its coordinates are read as longitude and ",
              "latitude (EPSG:4326).", call. = FALSE)
      x <- sf::st_set_crs(x, "EPSG:4326")
    } else if (!isTRUE(sf::st_is_longlat(x))) {
      x <- sf::st_transform(x, "EPSG:4326")
    }
  } else if (isTRUE(sf::st_is_longlat(x))) {
    warning(arg, " has longitude and latitude coordinates. With use_sphere = ",
            "FALSE, they are read as planar coordinates, so the results are ",
            "in degrees. Use use_sphere = TRUE for meters on the sphere.",
            call. = FALSE)
  }
  x
}

## Converts an sfc object to a list with one matrix per unit, for the
## Hausdorff distances in src/haus-r2.cpp and src/haus-s2.cpp. The columns are
## x, y, and a part id (plus an `is_poly` column when `include_poly = TRUE`).
## Each ring, line, and point is a separate part, so the matrices keep every
## polygon, hole, and island of a unit. An empty unit gives a matrix with no
## rows.
.to_coord_list <- function(x, include_poly = FALSE) {
  n <- length(x)
  ncols <- if (include_poly) 4L else 3L
  if (n == 0L) return(list())
  crd <- wk::wk_coords(x)
  if (nrow(crd) == 0L) {
    return(replicate(n, matrix(numeric(0), nrow = 0L, ncol = ncols),
                     simplify = FALSE))
  }
  new_part <- c(TRUE, diff(crd$feature_id) != 0 | diff(crd$part_id) != 0 |
                  diff(crd$ring_id) != 0)
  part <- as.numeric(cumsum(new_part))
  rows <- split(seq_len(nrow(crd)),
                factor(crd$feature_id, levels = seq_len(n)))
  if (!include_poly) {
    return(unname(lapply(rows, function(i) {
      cbind(crd$x[i], crd$y[i], part[i])
    })))
  }
  ## wk::wk_coords() increments ring_id whenever entering a polygon ring and
  ## leaves ring_id unchanged across points and linestrings.
  part_start <- which(new_part)
  is_poly_part <- c(crd$ring_id[1L] > 0L, diff(crd$ring_id) > 0L)[part_start]
  is_poly <- as.numeric(is_poly_part[part])
  unname(lapply(rows, function(i) {
    cbind(crd$x[i], crd$y[i], part[i], is_poly[i])
  }))
}

.to_list_mat <- function(x) {
  .to_coord_list(x, include_poly = FALSE)
}

.to_polygon_list <- function(x) {
  .to_coord_list(x, include_poly = TRUE)
}

## `bad` holds the indices of the units of argument `arg`.
.warn_non_hemisphere <- function(bad, arg = "x") {
  if (length(bad)) {
    warning("These units of ", arg, " are not inside an open hemisphere: ",
            paste(bad, collapse = ", "), ". Their Chebyshev radius is at ",
            "least pi * R / 2. The algorithm needs a radius below ",
            "pi * R / 2, so the function returns NA for them.", call. = FALSE)
  }
}

## Warns about the units of argument `arg` that are not inside an open
## hemisphere. `crd` is the output of wk::wk_coords() for `n` units. The check
## runs only if `has_na` is TRUE, because only these units give NA for a
## non-empty unit.
.warn_cr_sphere_if_needed <- function(crd, n, has_na, arg = "x") {
  if (!has_na || n == 0L) return(invisible(NULL))
  check <- mbc_sphere_cpp(crd$x, crd$y, crd$feature_id, n,
                          earth_radius = s2::s2_earth_radius_meters())
  .warn_non_hemisphere(which(is.na(check$radius) & check$non_empty), arg)
}

## ---- Hausdorff distance -----------------------------------------------------

##' @title Hausdorff distance between spatial units
##'
##' @description Computes the Hausdorff distance between all pairs of spatial
##'   units in \code{x}, or between each unit in \code{x} and each unit in
##'   \code{y}.
##'
##' @details The \code{method} argument selects how polygons, lines, and points
##'   are represented:
##'   \itemize{
##'     \item \code{"discrete"} (the default): computes the discrete
##'       vertex-to-edge Hausdorff distance, matching GEOS and
##'       \code{sf::st_distance(which = "Hausdorff")} in the plane. For each
##'       vertex of a unit, it finds the distance to the nearest edge or point
##'       of the other unit and takes the maximum across both directions. This
##'       matches continuous Hausdorff for point sets (\code{POINT} and
##'       \code{MULTIPOINT}) and for a single planar \code{POINT} against a
##'       \code{POLYGON}, but only approximates continuous Hausdorff for
##'       polygon-to-polygon or \code{MULTIPOINT}-to-polygon pairs because it
##'       does not sample edge interiors or polygon interiors. On the sphere,
##'       even a single \code{POINT} against a \code{POLYGON} can differ if the
##'       point's antipode enters or projects onto the polygon.
##'     \item \code{"boundary"}: computes the continuous Hausdorff distance
##'       between unit boundaries (all rings of polygons including holes and
##'       islands, linestrings, and points) using branch-and-bound segment or
##'       arc subdivision within tolerance \code{tol}.
##'     \item \code{"polygon"}: computes the continuous Hausdorff distance
##'       between spatial units treated as closed sets, where polygons include
##'       their 2D interiors (excluding holes) and points and linestrings retain
##'       their 0D and 1D geometry, within tolerance \code{tol}.
##'   }
##'
##'   For \code{"discrete"}, the C++ code uses the early break of Taha and
##'   Hanbury (2015), and the result is exact for the discrete vertex-to-edge
##'   definition. Continuous methods (\code{"boundary"} and \code{"polygon"})
##'   stop when the branch-and-bound gap is at most \code{tol * max(1, scale)},
##'   where \code{scale} is the larger bounding-box diagonal of the two units
##'   in the plane or the larger bounding-cap diameter in meters on the sphere.
##'
##'   On the sphere (\code{use_sphere = TRUE}), edges are great-circle arcs, as
##'   in \code{s2}. For \code{method = "polygon"}, each polygon unit must fit
##'   inside a spherical cap of radius below \eqn{0.45 \pi R}{0.45 pi R}
##'   (about 81 degrees); each ring encloses its smaller region (under
##'   \eqn{2\pi}{2 pi} steradians), matching \code{sf} and \code{s2}, and holes
##'   and islands compose by even-odd parity. Polygon units that exceed this
##'   cap limit return \code{NA} with a warning (use \code{method = "discrete"}
##'   or \code{"boundary"} for global-scale polygons). Continuous spherical
##'   methods also cap branch-and-bound iterations per pair so degenerate
##'   inputs cannot stall; if any pair does not reach \code{tol}, the function
##'   returns the best lower bound found and issues a warning.
##'
##' @param x an \code{sf} or \code{sfc} object.
##' @param y an optional \code{sf} or \code{sfc} object with the same CRS as
##'   \code{x}. If \code{NULL} (the default), the function computes the
##'   distances between the units in \code{x}.
##' @param use_sphere a logical. If \code{FALSE} (the default), the
##'   coordinates are planar. If \code{TRUE}, the distances are on a sphere of
##'   radius \code{s2::s2_earth_radius_meters()}. Projected input is
##'   transformed to EPSG:4326. Input without a CRS is read as longitude and
##'   latitude, with a warning.
##' @param method a character string selecting the Hausdorff distance
##'   definition: \code{"discrete"} (the default), \code{"boundary"}, or
##'   \code{"polygon"}.
##' @param tol a positive numeric scalar giving the convergence tolerance for
##'   the continuous methods (\code{"boundary"} and \code{"polygon"}). The
##'   effective tolerance for each pair is \code{tol * max(1, scale)}, floored
##'   above floating-point rounding limits. Ignored when
##'   \code{method = "discrete"}.
##'
##' @return A numeric matrix. Without \code{y}, it is a symmetric
##'   \eqn{n \times n}{n x n} matrix, where \eqn{n} is the number of units in
##'   \code{x}. With \code{y}, it is an \eqn{n \times m}{n x m} matrix, where
##'   \eqn{m} is the number of units in \code{y}. In the plane, the distances
##'   are in the units of the CRS. On the sphere, they are in meters. An empty
##'   unit, or a unit with a non-finite coordinate (for example, a
##'   \code{MULTIPOINT} with an \code{EMPTY} member), gives \code{NA} in its row
##'   and column.
##'
##' @references Taha, A. A., & Hanbury, A. (2015). An efficient algorithm for
##'   calculating the exact Hausdorff distance. \emph{IEEE Transactions on
##'   Pattern Analysis and Machine Intelligence}, 37(11), 2153--2163.
##'   \doi{10.1109/TPAMI.2015.2408351}
##'
##' @seealso \code{\link{dist_cr}}
##'
##' @examples
##' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"),
##'                   quiet = TRUE)
##' nc <- sf::st_transform(nc[1:5, ], 32119)
##'
##' ## Planar discrete vertex-to-edge distances, in meters
##' dist_hausdorff(nc)
##' dist_hausdorff(nc[1:2, ], nc[3:5, ])
##'
##' ## Continuous boundary and filled-polygon distances, in meters
##' dist_hausdorff(nc[1:3, ], method = "boundary")
##' dist_hausdorff(nc[1:3, ], method = "polygon")
##'
##' ## Distances on the sphere, in meters
##' dist_hausdorff(nc, use_sphere = TRUE)
##' dist_hausdorff(nc[1:3, ], use_sphere = TRUE, method = "polygon")
##' @export
dist_hausdorff <- function(x, y = NULL, use_sphere = FALSE,
                           method = c("discrete", "boundary", "polygon"),
                           tol = 1e-8) {
  method <- match.arg(method)
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
      !is.finite(tol) || tol <= 0) {
    stop("tol must be a single finite positive number.", call. = FALSE)
  }
  if (!is.null(y)) .check_same_crs(x, y)
  x <- .prep_geom(x, use_sphere = use_sphere, arg = "x")
  if (!is.null(y)) y <- .prep_geom(y, use_sphere = use_sphere, arg = "y")
  to_mats <- if (method == "polygon") .to_polygon_list else .to_list_mat
  if (use_sphere) {
    radius <- s2::s2_earth_radius_meters()
    if (method == "discrete") {
      if (is.null(y)) {
        return(dist_haus_s2(to_mats(x), earth_radius = radius))
      }
      return(cdist_haus_s2_lists(to_mats(x), to_mats(y),
                                 earth_radius = radius))
    }
    if (is.null(y)) {
      return(dist_haus_s2_method(to_mats(x), method = method, tol = tol,
                                 earth_radius = radius))
    }
    return(cdist_haus_s2_lists_method(to_mats(x), to_mats(y),
                                      method = method, tol = tol,
                                      earth_radius = radius))
  }
  if (method == "discrete") {
    if (is.null(y)) {
      return(dist_haus_r2(to_mats(x)))
    }
    return(cdist_haus_r2_lists(to_mats(x), to_mats(y)))
  }
  if (is.null(y)) {
    return(dist_haus_r2_method(to_mats(x), method = method, tol = tol))
  }
  cdist_haus_r2_lists_method(to_mats(x), to_mats(y),
                             method = method, tol = tol)
}

## ---- Chebyshev centers and radii ---------------------------------------------

##' @title Chebyshev centers and radii of spatial units
##'
##' @description Computes the center and radius of the minimum enclosing ball
##'   (the Chebyshev center and radius) of each spatial unit in \code{x}.
##'
##' @details The function uses Welzl's algorithm (Welzl, 1991) to compute the
##'   exact minimum enclosing ball of the vertices of each unit. On the sphere,
##'   edges are great-circle arcs, as in \code{s2}, so this cap is also the
##'   minimum enclosing cap of the unit. Do not use \code{s2::s2_bounds_cap()}
##'   for this purpose: it bounds the latitude-longitude rectangle of the
##'   unit, and its cap is not minimal.
##'
##' @inheritParams dist_hausdorff
##'
##' @return A list with two elements:
##'   \item{center}{an \code{sfc} object with one point for each unit. In the
##'     plane, its CRS is the CRS of \code{x}. On the sphere, it is the
##'     longitude and latitude CRS of \code{x}, or EPSG:4326 if \code{x} has a
##'     projected CRS or no CRS.}
##'   \item{radius}{a numeric vector. In the plane, the radii are in the units
##'     of the CRS. On the sphere, they are in meters.}
##'   The center and radius are \code{NA} for an empty unit, or for a unit with
##'   a non-finite coordinate. On the sphere, they are also \code{NA}, with a
##'   warning, for a unit that is not inside an open hemisphere. The Chebyshev
##'   radius of such a unit is at least \eqn{\pi R / 2}{pi R / 2}, and the
##'   algorithm needs a smaller radius.
##'
##' @references Welzl, E. (1991). Smallest enclosing disks (balls and
##'   ellipsoids). In H. Maurer (Ed.), \emph{New Results and New Trends in
##'   Computer Science} (pp. 359--370). Springer. \doi{10.1007/BFb0038202}
##'
##' @seealso \code{\link{dist_cr}}
##'
##' @examples
##' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"),
##'                   quiet = TRUE)
##' nc <- sf::st_transform(nc[1:5, ], 32119)
##' get_mbc_params(nc)
##' get_mbc_params(nc, use_sphere = TRUE)
##' @export
get_mbc_params <- function(x, use_sphere = FALSE) {
  x <- .prep_geom(x, use_sphere = use_sphere, arg = "x")
  n <- length(x)
  crd <- wk::wk_coords(x)
  if (use_sphere) {
    res <- mbc_sphere_cpp(crd$x, crd$y, crd$feature_id, n,
                          earth_radius = s2::s2_earth_radius_meters())
    .warn_non_hemisphere(which(is.na(res$radius) & res$non_empty), "x")
    cx <- res$lon
    cy <- res$lat
  } else {
    res <- mbc_planar_cpp(crd$x, crd$y, crd$feature_id, n)
    cx <- res$cx
    cy <- res$cy
  }
  centers <- sf::st_sfc(
    lapply(seq_len(n), function(i) sf::st_point(c(cx[i], cy[i]))),
    crs = sf::st_crs(x)
  )
  list(center = centers, radius = res$radius)
}

## ---- CR distance ------------------------------------------------------------

##' @title Center-radii (CR) distance between spatial units
##'
##' @description Computes the center-radii (CR) distance between all pairs of
##'   spatial units in \code{x}, or between each unit in \code{x} and each unit
##'   in \code{y}.
##'
##' @details The CR distance between units \eqn{i} and \eqn{j} is
##'   \deqn{d(c_i, c_j) + |r_i - r_j|,}{d(c_i, c_j) + |r_i - r_j|,} where
##'   \eqn{c} and \eqn{r} are the
##'   Chebyshev center and radius of a unit (see \code{\link{get_mbc_params}}).
##'   In the plane, \eqn{d} is the Euclidean distance. On the sphere, \eqn{d}
##'   is the great-circle distance.
##'
##' @inheritParams dist_hausdorff
##'
##' @return A numeric matrix. Without \code{y}, it is a symmetric
##'   \eqn{n \times n}{n x n} matrix with a zero diagonal, where \eqn{n} is the
##'   number of units in \code{x}. With \code{y}, it is an
##'   \eqn{n \times m}{n x m} matrix, where \eqn{m} is the number of units in
##'   \code{y}. In the plane, the distances are in the units of the CRS. On the
##'   sphere, they are in meters. An empty unit, or a unit with a non-finite
##'   coordinate, gives \code{NA} off the diagonal. On the sphere, a unit that
##'   is not inside an open hemisphere also gives \code{NA}, with a warning.
##'
##' @seealso \code{\link{get_mbc_params}}, \code{\link{dist_cr_sphere}},
##'   \code{\link{dist_hausdorff}}
##'
##' @examples
##' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"),
##'                   quiet = TRUE)
##' nc <- sf::st_transform(nc[1:5, ], 32119)
##'
##' ## Planar distances, in meters (the unit of EPSG:32119)
##' dist_cr(nc)
##' dist_cr(nc[1:2, ], nc[3:5, ])
##'
##' ## Distances on the sphere, in meters
##' dist_cr(nc, use_sphere = TRUE)
##' @export
dist_cr <- function(x, y = NULL, use_sphere = FALSE) {
  if (!is.null(y)) .check_same_crs(x, y)
  x <- .prep_geom(x, use_sphere = use_sphere, arg = "x")
  n1 <- length(x)
  crd1 <- wk::wk_coords(x)
  if (is.null(y)) {
    if (!use_sphere) {
      return(cr_planar_cpp(crd1$x, crd1$y, crd1$feature_id, n1))
    }
    out <- cr_sphere_cpp(crd1$x, crd1$y, crd1$feature_id, n1,
                         earth_radius = s2::s2_earth_radius_meters())
    .warn_cr_sphere_if_needed(crd1, n1, anyNA(out) || n1 == 1L)
    return(out)
  }
  y <- .prep_geom(y, use_sphere = use_sphere, arg = "y")
  n2 <- length(y)
  crd2 <- wk::wk_coords(y)
  if (!use_sphere) {
    return(cr_cross_planar_cpp(crd1$x, crd1$y, crd1$feature_id, n1,
                               crd2$x, crd2$y, crd2$feature_id, n2))
  }
  out <- cr_cross_sphere_cpp(crd1$x, crd1$y, crd1$feature_id, n1,
                             crd2$x, crd2$y, crd2$feature_id, n2,
                             earth_radius = s2::s2_earth_radius_meters())
  check <- anyNA(out) || n1 == 0L || n2 == 0L
  .warn_cr_sphere_if_needed(crd1, n1, check, "x")
  .warn_cr_sphere_if_needed(crd2, n2, check, "y")
  out
}

##' @title Spherical CR distance with the exceedance bound
##'
##' @description Computes the center-radii (CR) distance on the sphere,
##'   bounded by the exceedance bound, for all pairs of spatial units in
##'   \code{x}.
##'
##' @details On a sphere of radius \eqn{R}, the exceedance bound for units
##'   \eqn{i} and \eqn{j} is
##'   \deqn{e_{ij} = \pi R - \min(r_i, r_j),}{e_ij = pi R - min(r_i, r_j),}
##'   where \eqn{r} is the Chebyshev radius of a unit (see
##'   \code{\link{get_mbc_params}}). Let \eqn{u_{ij}}{u_ij} be the spherical
##'   CR distance from \code{\link{dist_cr}}. Then:
##'   \itemize{
##'     \item \code{dist_cr_sphere(x)} returns
##'       \eqn{\min(e_{ij}, u_{ij})}{min(e_ij, u_ij)}.
##'     \item \code{dist_cr_sphere(x, exceedance = TRUE)} returns
##'       \eqn{e_{ij}}{e_ij}.
##'     \item \code{dist_cr_all(x)} returns the three matrices from one
##'       computation of the centers and radii.
##'   }
##'   The coordinates are always on the sphere, with the same rules as
##'   \code{use_sphere = TRUE} in \code{\link{dist_cr}}.
##'
##' @param x an \code{sf} or \code{sfc} object.
##' @param exceedance a logical. If \code{TRUE}, the function returns only the
##'   exceedance bound. The default is \code{FALSE}.
##'
##' @return \code{dist_cr_sphere()} returns a symmetric
##'   \eqn{n \times n}{n x n} numeric matrix, in meters, where \eqn{n} is the
##'   number of units in \code{x}. With \code{exceedance = TRUE}, the diagonal
##'   is \eqn{\pi R - r_i}{pi R - r_i}. Otherwise, the diagonal is zero.
##'
##'   \code{dist_cr_all()} returns a list with three such matrices:
##'   \item{upper}{the spherical CR distance \eqn{u_{ij}}{u_ij}.}
##'   \item{exceedance}{the exceedance bound \eqn{e_{ij}}{e_ij}.}
##'   \item{cr}{\eqn{\min(e_{ij}, u_{ij})}{min(e_ij, u_ij)}.}
##'
##'   An empty unit, a unit with a non-finite coordinate, or a unit that is not
##'   inside an open hemisphere gives \code{NA}. The last case also gives a
##'   warning.
##'
##' @seealso \code{\link{dist_cr}}, \code{\link{get_mbc_params}}
##'
##' @examples
##' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"),
##'                   quiet = TRUE)
##' nc <- nc[1:5, ]
##' dist_cr_sphere(nc)
##' dist_cr_sphere(nc, exceedance = TRUE)
##' str(dist_cr_all(nc))
##' @export
dist_cr_sphere <- function(x, exceedance = FALSE) {
  x <- .prep_geom(x, use_sphere = TRUE, arg = "x")
  n <- length(x)
  crd <- wk::wk_coords(x)
  out <- cr_sphere_meb_cpp(crd$x, crd$y, crd$feature_id, n,
                           exceedance = exceedance,
                           earth_radius = s2::s2_earth_radius_meters())
  .warn_cr_sphere_if_needed(crd, n, anyNA(out))
  out
}

##' @rdname dist_cr_sphere
##' @export
dist_cr_all <- function(x) {
  x <- .prep_geom(x, use_sphere = TRUE, arg = "x")
  n <- length(x)
  crd <- wk::wk_coords(x)
  out <- cr_all_cpp(crd$x, crd$y, crd$feature_id, n,
                    earth_radius = s2::s2_earth_radius_meters())
  .warn_cr_sphere_if_needed(crd, n, anyNA(out$exceedance))
  out
}
