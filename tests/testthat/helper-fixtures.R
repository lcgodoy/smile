## Shared fixtures and reference implementations for the tests of
## dist_hausdorff(), dist_cr(), dist_cr_sphere(), dist_cr_all(), and
## get_mbc_params().
##
## The references use sf wherever sf has an equivalent:
## - Planar Hausdorff: sf::st_distance(which = "Hausdorff") (GEOS).
## - Spherical Hausdorff: sf::st_distance() (s2) from each vertex to the
##   linework of the other unit.
## - Planar centers and radii: sf::st_minimum_bounding_circle() (GEOS).
## - Distances between centers: sf::st_distance() (GEOS or s2).
## sf has no minimum enclosing cap on the sphere, so the spherical centers and
## radii come from a pure-R port of Welzl's algorithm.

## ---- Utilities --------------------------------------------------------------

## Removes the units class from a sf::st_distance() matrix.
drop_u <- function(m) {
  m <- unclass(m)
  attr(m, "units") <- NULL
  m
}

skip_if_no_mbc <- function() {
  skip_if_not(exists("st_minimum_bounding_circle", envir = asNamespace("sf"),
                     inherits = FALSE),
              "sf::st_minimum_bounding_circle() needs sf >= 1.0-20")
  geos <- sf::sf_extSoftVersion()[["GEOS"]]
  skip_if(package_version(geos) < "3.8.0",
          "sf::st_minimum_bounding_circle() needs GEOS >= 3.8.0")
}

skip_if_no_s2 <- function() {
  skip_if_not(isTRUE(sf::sf_use_s2()), "sf does not use s2")
}

read_nc <- function() {
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  sf::st_geometry(nc)
}

## One POINT for each vertex of the single unit in `g` (an sfc of length 1).
vertices_sfc <- function(g) {
  crd <- wk::wk_coords(g)
  sf::st_sfc(
    lapply(seq_len(nrow(crd)), function(k) sf::st_point(c(crd$x[k], crd$y[k]))),
    crs = sf::st_crs(g)
  )
}

## ---- Fixtures (from setgp/tests/testthat/helper-r-reference.R) --------------

make_planar_fixtures <- function() {
  p_sq <- sf::st_polygon(list(rbind(
    c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0)
  )))
  p_notch <- sf::st_polygon(list(rbind(
    c(0, 0), c(10, 0), c(10, 10), c(5, 5), c(0, 10), c(0, 0)
  )))
  p_hole <- sf::st_polygon(list(
    rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0)),
    rbind(c(4, 4), c(4, 6), c(6, 6), c(6, 4), c(4, 4))
  ))
  p_collinear <- sf::st_polygon(list(rbind(
    c(0, 0), c(3, 0), c(6, 0), c(6, 4), c(0, 0)
  )))
  p_multi <- sf::st_multipolygon(list(
    list(rbind(c(12, 1), c(14, 1), c(14, 3), c(12, 3), c(12, 1))),
    list(rbind(c(15, 2), c(16, 2), c(16, 4), c(15, 4), c(15, 2)))
  ))
  pt <- sf::st_point(c(-4, 7))
  mpt <- sf::st_multipoint(rbind(c(-2, -1), c(1, -3), c(-1, 2)))
  ls <- sf::st_linestring(rbind(c(3, -2), c(7, -1), c(9, 1)))
  gc <- sf::st_geometrycollection(list(
    sf::st_point(c(20, 20)),
    sf::st_polygon(list(rbind(c(22, 20), c(24, 20), c(24, 22), c(22, 22),
                              c(22, 20))))
  ))

  ## Large projected coordinates (EPSG:3857 scale) test the precision.
  shift <- c(3400000, 5200000)
  p_big1 <- sf::st_polygon(list(rbind(
    shift + c(0, 0), shift + c(1200, 0), shift + c(900, 800), shift + c(0, 0)
  )))
  p_big2 <- sf::st_polygon(list(rbind(
    shift + c(2500, -1000), shift + c(3500, -1000),
    shift + c(3500, 500), shift + c(2500, 500), shift + c(2500, -1000)
  )))

  list(
    mixed = sf::st_sfc(p_sq, p_notch, p_hole, p_collinear, p_multi, pt, mpt,
                       ls, gc, crs = "EPSG:3857"),
    polygons = sf::st_sfc(p_sq, p_notch, p_hole, p_multi, crs = "EPSG:3857"),
    targets = sf::st_sfc(
      sf::st_polygon(list(rbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1),
                                c(-1, -1)))),
      sf::st_point(c(5, 5)),
      p_collinear,
      crs = "EPSG:3857"
    ),
    points = sf::st_sfc(
      sf::st_point(c(0, 0)),
      sf::st_point(c(3, 4)),
      sf::st_point(c(-5, 12)),
      crs = "EPSG:3857"
    ),
    large_proj = sf::st_sfc(p_big1, p_big2, crs = "EPSG:3857"),
    with_empty = sf::st_sfc(p_sq, sf::st_polygon(), pt, crs = "EPSG:3857")
  )
}

make_sphere_fixtures <- function() {
  poly1 <- sf::st_polygon(list(rbind(
    c(-46.8, -23.7), c(-46.4, -23.7), c(-46.4, -23.4), c(-46.8, -23.4),
    c(-46.8, -23.7)
  )))
  poly2 <- sf::st_polygon(list(rbind(
    c(-74.1, 40.6), c(-73.7, 40.6), c(-73.7, 40.9), c(-74.1, 40.9),
    c(-74.1, 40.6)
  )))
  poly3 <- sf::st_polygon(list(
    rbind(c(2.2, 48.7), c(2.5, 48.7), c(2.5, 49.0), c(2.2, 49.0), c(2.2, 48.7)),
    rbind(c(2.3, 48.8), c(2.3, 48.9), c(2.4, 48.9), c(2.4, 48.8), c(2.3, 48.8))
  ))
  multi_poly <- sf::st_multipolygon(list(
    list(rbind(c(139.6, 35.6), c(139.8, 35.6), c(139.8, 35.8), c(139.6, 35.8),
               c(139.6, 35.6))),
    list(rbind(c(140.0, 35.7), c(140.2, 35.7), c(140.2, 35.9), c(140.0, 35.9),
               c(140.0, 35.7)))
  ))
  pt <- sf::st_point(c(-0.1276, 51.5072))
  mpt <- sf::st_multipoint(rbind(c(-122.4, 37.7), c(-122.2, 37.8),
                                 c(-122.3, 37.9)))
  ls <- sf::st_linestring(rbind(c(151.1, -33.9), c(151.3, -33.8),
                                c(151.4, -33.7)))

  ## Points on a small circle with a radius just under pi / 2.
  near_hemi <- sf::st_multipoint(rbind(
    c(-80, 10), c(80, 10), c(0, 75), c(0, -5)
  ))

  ## Points that are not inside an open hemisphere (radius >= pi / 2).
  over_hemi <- sf::st_multipoint(rbind(
    c(0, 0), c(180, 0), c(0, 80), c(0, -80)
  ))

  list(
    mixed = sf::st_sfc(poly1, poly2, poly3, multi_poly, pt, mpt, ls,
                       near_hemi, crs = "EPSG:4326"),
    targets = sf::st_sfc(poly1, pt, ls, sf::st_point(c(12.49, 41.90)),
                         crs = "EPSG:4326"),
    with_empty = sf::st_sfc(poly1, sf::st_polygon(), pt, crs = "EPSG:4326"),
    over_hemi = sf::st_sfc(poly1, over_hemi, pt, crs = "EPSG:4326"),
    no_crs = sf::st_sfc(poly1, poly2, pt)
  )
}

## ---- Spherical Hausdorff reference (sf + s2) --------------------------------

## The edges of a unit: the rings of a polygon become lines. Points and lines
## stay as they are.
linework_sfc <- function(g) {
  if (any(sf::st_is(g, c("POLYGON", "MULTIPOLYGON")))) {
    return(sf::st_cast(g, "MULTILINESTRING"))
  }
  g
}

## Largest s2 distance from a vertex of `a` to the edges of `b`.
directed_haus_s2 <- function(a, b) {
  max(drop_u(sf::st_distance(vertices_sfc(a), linework_sfc(b))))
}

## Discrete Hausdorff distance on the sphere, in meters, as GEOS defines it in
## the plane: the larger of the two directed distances. An empty unit gives NA.
haus_ref_s2 <- function(x, y = NULL) {
  x <- sf::st_geometry(x)
  yy <- if (is.null(y)) x else sf::st_geometry(y)
  out <- matrix(NA_real_, length(x), length(yy))
  for (i in seq_along(x)) {
    for (j in seq_along(yy)) {
      if (sf::st_is_empty(x[i]) || sf::st_is_empty(yy[j])) next
      out[i, j] <- max(directed_haus_s2(x[i], yy[j]),
                       directed_haus_s2(yy[j], x[i]))
    }
  }
  out
}

## ---- Planar centers and radii reference (sf + GEOS) -------------------------

## sf::st_minimum_bounding_circle() returns a regular polygon, with its
## vertices on the circle. The center is the mean of the vertices. The radius
## is the largest distance from the center to a vertex. The centers have no
## CRS, so sf::st_distance() between them is Euclidean.
mbc_ref_planar <- function(x) {
  x <- sf::st_geometry(x)
  circ <- suppressWarnings(sf::st_minimum_bounding_circle(x))
  n <- length(x)
  cx <- cy <- radius <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    if (sf::st_is_empty(circ[i])) next
    xy <- sf::st_coordinates(circ[i])[, c("X", "Y"), drop = FALSE]
    if (nrow(xy) == 1L) {
      cx[i] <- xy[1L, 1L]
      cy[i] <- xy[1L, 2L]
      radius[i] <- 0
      next
    }
    ring <- xy[-nrow(xy), , drop = FALSE]  # drop the closing vertex
    ctr <- colMeans(ring)
    cx[i] <- ctr[1L]
    cy[i] <- ctr[2L]
    radius[i] <- max(sqrt((ring[, 1L] - ctr[1L])^2 + (ring[, 2L] - ctr[2L])^2))
  }
  centers <- sf::st_sfc(
    lapply(seq_len(n), function(i) sf::st_point(c(cx[i], cy[i])))
  )
  list(center = centers, radius = radius)
}

## ---- Spherical centers and radii reference (pure R) -------------------------
## Port of .welzl_cap_r() and .min_cap_r() from
## setgp/tests/testthat/helper-r-reference.R.

.ll_to_xyz_r <- function(lon, lat) {
  lon <- lon * pi / 180
  lat <- lat * pi / 180
  cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
}

.xyz_to_ll_r <- function(u) {
  c(atan2(u[2], u[1]), atan2(u[3], sqrt(u[1]^2 + u[2]^2))) * 180 / pi
}

.cross3_r <- function(a, b) {
  c(a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])
}

.chord_r <- function(P, u) {
  sqrt((P[, 1] - u[1])^2 + (P[, 2] - u[2])^2 + (P[, 3] - u[3])^2)
}

.cap2_r <- function(a, b) {
  u <- a + b
  u <- u / sqrt(sum(u^2))
  list(center = u, r = max(.chord_r(rbind(a, b), u)))
}

.cap3_r <- function(a, b, c) {
  m <- a + b + c
  m <- m / sqrt(sum(m^2))
  e1 <- .cross3_r(if (abs(m[3]) < 0.9) c(0, 0, 1) else c(1, 0, 0), m)
  e1 <- e1 / sqrt(sum(e1^2))
  e2 <- .cross3_r(m, e1)
  ba <- b - a
  ca <- c - a
  a11 <- sum(ba * e1)
  a12 <- sum(ba * e2)
  a21 <- sum(ca * e1)
  a22 <- sum(ca * e2)
  da <- sum((a - m)^2)
  r1 <- (sum((b - m)^2) - da) / 2
  r2 <- (sum((c - m)^2) - da) / 2
  det <- a11 * a22 - a12 * a21
  if (det != 0) {
    u <- m + ((r1 * a22 - a12 * r2) * e1 + (a11 * r2 - r1 * a21) * e2) / det
  } else {
    u <- .cross3_r(ba, ca)
    if (sum(u * m) < 0) u <- -u
  }
  u <- u / sqrt(sum(u^2))
  list(center = u, r = max(.chord_r(rbind(a, b, c), u)))
}

.first_outside_r <- function(P, cap, from, to, tol) {
  if (from > to) return(NA_integer_)
  out <- .chord_r(P[from:to, , drop = FALSE], cap$center) > cap$r + tol
  k <- match(TRUE, out)
  if (is.na(k)) NA_integer_ else from + k - 1L
}

.welzl_cap_r <- function(P, tol = 1e-13) {
  with_2 <- function(m, p, q) {
    cap <- .cap2_r(p, q)
    k <- .first_outside_r(P, cap, 1L, m, tol)
    while (!is.na(k)) {
      cap <- .cap3_r(p, q, P[k, ])
      k <- .first_outside_r(P, cap, k + 1L, m, tol)
    }
    cap
  }
  with_1 <- function(m, p) {
    cap <- list(center = p, r = 0)
    j <- .first_outside_r(P, cap, 1L, m, tol)
    while (!is.na(j)) {
      cap <- with_2(j - 1L, p, P[j, ])
      j <- .first_outside_r(P, cap, j + 1L, m, tol)
    }
    cap
  }
  n <- nrow(P)
  cap <- list(center = P[1L, ], r = 0)
  i <- .first_outside_r(P, cap, 2L, n, tol)
  while (!is.na(i)) {
    cap <- with_1(i - 1L, P[i, ])
    i <- .first_outside_r(P, cap, i + 1L, n, tol)
  }
  cap
}

.min_cap_r <- function(P, tol = 1e-13) {
  na <- list(center = rep(NA_real_, 3), angle = NA_real_)
  if (nrow(P) == 0L) return(na)
  P <- P[order((seq_len(nrow(P)) * 0.6180339887498949) %% 1), , drop = FALSE]
  u <- .welzl_cap_r(P, tol)$center
  cr <- cbind(P[, 2] * u[3] - P[, 3] * u[2],
              P[, 3] * u[1] - P[, 1] * u[3],
              P[, 1] * u[2] - P[, 2] * u[1])
  angle <- max(atan2(sqrt(rowSums(cr^2)), drop(P %*% u)))
  if (!is.finite(angle) || angle >= pi / 2) return(na)
  list(center = u, angle = angle)
}

## Spherical centers (EPSG:4326) and radii (meters). `x` must have a
## longitude and latitude CRS, or a projected CRS. It gives no warnings.
mbc_ref_welzl_s2 <- function(x) {
  x <- sf::st_geometry(x)
  if (!isTRUE(sf::st_is_longlat(x))) x <- sf::st_transform(x, "EPSG:4326")
  crd <- wk::wk_coords(x)
  rows <- split(seq_len(nrow(crd)),
                factor(crd$feature_id, levels = seq_along(x)))
  caps <- lapply(rows, function(i) {
    .min_cap_r(.ll_to_xyz_r(crd$x[i], crd$y[i]))
  })
  angle <- unname(vapply(caps, `[[`, numeric(1), "angle"))
  centers <- sf::st_sfc(
    lapply(caps, function(cp) sf::st_point(.xyz_to_ll_r(cp$center))),
    crs = sf::st_crs(x)
  )
  list(center = centers, radius = angle * s2::s2_earth_radius_meters())
}

## ---- CR distance reference ---------------------------------------------------

## CR distance from the reference centers and radii, with sf::st_distance()
## between the centers (Euclidean in the plane, s2 on the sphere). Without
## `y`, the diagonal is zero.
cr_ref <- function(x, y = NULL, use_sphere = FALSE) {
  mbc <- if (use_sphere) mbc_ref_welzl_s2 else mbc_ref_planar
  m1 <- mbc(x)
  m2 <- if (is.null(y)) m1 else mbc(y)
  out <- drop_u(sf::st_distance(m1$center, m2$center)) +
    abs(outer(m1$radius, m2$radius, "-"))
  if (is.null(y)) diag(out) <- 0
  out
}

## Spherical CR distance with the exceedance bound, from the references.
cr_all_ref <- function(x) {
  r <- mbc_ref_welzl_s2(x)$radius
  upper <- cr_ref(x, use_sphere = TRUE)
  exceedance <- pi * s2::s2_earth_radius_meters() - outer(r, r, pmin)
  list(upper = upper, exceedance = exceedance, cr = pmin(exceedance, upper))
}
