## ---- Plane: compare with sf::st_distance(which = "Hausdorff") ----------------

test_that("planar dist_hausdorff() matches sf on the fixtures", {
  fix <- make_planar_fixtures()

  for (nm in c("mixed", "polygons", "targets", "points", "large_proj")) {
    geom <- fix[[nm]]
    got <- dist_hausdorff(geom)
    ref <- drop_u(sf::st_distance(geom, which = "Hausdorff"))

    expect_equal(got, ref, tolerance = 1e-8, info = nm)
    expect_equal(diag(got), rep(0, length(geom)), info = nm)
    expect_equal(got, t(got), info = nm)
  }

  got <- dist_hausdorff(fix$polygons, fix$targets)
  ref <- drop_u(sf::st_distance(fix$polygons, fix$targets, which = "Hausdorff"))
  expect_equal(got, ref, tolerance = 1e-8)
})

test_that("planar dist_hausdorff() matches sf on North Carolina counties", {
  nc_proj <- sf::st_transform(read_nc()[1:25], 32119)

  expect_equal(dist_hausdorff(nc_proj),
               drop_u(sf::st_distance(nc_proj, which = "Hausdorff")),
               tolerance = 1e-8)
  expect_equal(dist_hausdorff(nc_proj[1:12], nc_proj[13:25]),
               drop_u(sf::st_distance(nc_proj[1:12], nc_proj[13:25],
                                      which = "Hausdorff")),
               tolerance = 1e-8)
})

test_that("dist_hausdorff() accepts sf objects", {
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  nc_proj <- sf::st_transform(nc[1:5, ], 32119)
  expect_equal(dist_hausdorff(nc_proj),
               dist_hausdorff(sf::st_geometry(nc_proj)))
})

## ---- Sphere: compare with sf::st_distance() (s2) -----------------------------

test_that("spherical dist_hausdorff() matches sf/s2 on the fixtures", {
  skip_if_no_s2()
  sph <- make_sphere_fixtures()

  got <- dist_hausdorff(sph$mixed, use_sphere = TRUE)
  expect_equal(got, haus_ref_s2(sph$mixed), tolerance = 1e-7)
  expect_equal(diag(got), rep(0, length(sph$mixed)))
  expect_equal(got, t(got))

  got_x <- dist_hausdorff(sph$mixed, sph$targets, use_sphere = TRUE)
  expect_equal(got_x, haus_ref_s2(sph$mixed, sph$targets), tolerance = 1e-7)
})

test_that("spherical dist_hausdorff() matches sf/s2 on North Carolina counties", {
  skip_if_no_s2()
  nc_ll <- sf::st_transform(read_nc()[1:10], "EPSG:4326")

  expect_equal(dist_hausdorff(nc_ll, use_sphere = TRUE),
               haus_ref_s2(nc_ll), tolerance = 1e-7)
  expect_equal(dist_hausdorff(nc_ll[1:4], nc_ll[5:10], use_sphere = TRUE),
               haus_ref_s2(nc_ll[1:4], nc_ll[5:10]), tolerance = 1e-7)
})

test_that("spherical dist_hausdorff() between points is the great-circle distance", {
  skip_if_no_s2()
  pts <- sf::st_sfc(
    sf::st_point(c(-46.63, -23.55)),
    sf::st_point(c(-46.50, -23.50)),
    sf::st_point(c(-74.00, 40.71)),
    sf::st_point(c(139.69, 35.69)),
    crs = "EPSG:4326"
  )
  ## A tight tolerance also checks the radius of the sphere.
  expect_equal(dist_hausdorff(pts, use_sphere = TRUE),
               drop_u(sf::st_distance(pts)), tolerance = 1e-10)
  expect_equal(dist_hausdorff(pts[1:2], pts[3:4], use_sphere = TRUE),
               drop_u(sf::st_distance(pts[1:2], pts[3:4])), tolerance = 1e-10)
})

test_that("dist_hausdorff() transforms projected input when use_sphere = TRUE", {
  nc_proj <- sf::st_transform(read_nc()[1:6], 32119)
  nc_ll <- sf::st_transform(nc_proj, "EPSG:4326")

  expect_equal(dist_hausdorff(nc_proj, use_sphere = TRUE),
               dist_hausdorff(nc_ll, use_sphere = TRUE))
  expect_equal(dist_hausdorff(nc_proj[1:2], nc_proj[3:6], use_sphere = TRUE),
               dist_hausdorff(nc_ll[1:2], nc_ll[3:6], use_sphere = TRUE))
})

## ---- Empty units, zero-length input, CRS, and warnings ----------------------

test_that("dist_hausdorff() gives NA for empty units in both spaces", {
  fix <- make_planar_fixtures()
  got <- dist_hausdorff(fix$with_empty)
  expect_equal(got,
               drop_u(sf::st_distance(fix$with_empty, which = "Hausdorff")),
               tolerance = 1e-8)
  expect_true(all(is.na(got[2L, ])))
  expect_true(all(is.na(got[, 2L])))
  expect_false(anyNA(got[-2L, -2L]))

  got_x <- dist_hausdorff(fix$with_empty, fix$targets)
  expect_true(all(is.na(got_x[2L, ])))
  expect_false(anyNA(got_x[-2L, ]))

  sph <- make_sphere_fixtures()
  got_s <- dist_hausdorff(sph$with_empty, use_sphere = TRUE)
  expect_true(all(is.na(got_s[2L, ])))
  expect_true(all(is.na(got_s[, 2L])))
  expect_false(anyNA(got_s[-2L, -2L]))

  got_sx <- dist_hausdorff(sph$targets, sph$with_empty, use_sphere = TRUE)
  expect_true(all(is.na(got_sx[, 2L])))
  expect_false(anyNA(got_sx[, -2L]))

  all_empty <- sf::st_sfc(sf::st_polygon(), sf::st_point(), crs = "EPSG:3857")
  expect_true(all(is.na(dist_hausdorff(all_empty))))
  expect_true(all(is.na(dist_hausdorff(all_empty, use_sphere = TRUE))))
})

test_that("dist_hausdorff() gives NA for units with a non-finite coordinate", {
  ## wk::wk_coords() gives a NaN row for the EMPTY member of the MULTIPOINT.
  wkt <- c("MULTIPOINT ((0 0), EMPTY, (10 0))", "POINT (50 0)", "POINT (60 0)")
  g_ll <- sf::st_as_sfc(wkt, crs = "EPSG:4326")
  g_pl <- sf::st_as_sfc(wkt, crs = "EPSG:3857")

  for (got in list(dist_hausdorff(g_pl),
                   dist_hausdorff(g_ll, use_sphere = TRUE))) {
    expect_true(all(is.na(got[1L, ])))
    expect_true(all(is.na(got[, 1L])))
    expect_false(anyNA(got[-1L, -1L]))
  }
  expect_true(all(is.na(dist_hausdorff(g_ll[1L], g_ll[2:3], use_sphere = TRUE))))

  ## The same rule as the CR distance.
  expect_true(is.na(dist_cr(g_pl)[1L, 2L]))
  expect_true(is.na(dist_cr(g_ll, use_sphere = TRUE)[1L, 2L]))
})

test_that("dist_hausdorff() handles zero-length input", {
  fix <- make_planar_fixtures()
  zero <- sf::st_sfc(crs = "EPSG:3857")

  expect_equal(dim(dist_hausdorff(zero)), c(0L, 0L))
  expect_equal(dim(dist_hausdorff(zero, fix$targets)), c(0L, 3L))
  expect_equal(dim(dist_hausdorff(fix$targets, zero)), c(3L, 0L))
  expect_equal(dim(dist_hausdorff(zero, use_sphere = TRUE)), c(0L, 0L))
})

test_that("dist_hausdorff() stops when x and y have different CRS", {
  fix <- make_planar_fixtures()
  expect_error(
    dist_hausdorff(fix$polygons, sf::st_transform(fix$targets, 4326)),
    "same CRS"
  )
})

test_that("dist_hausdorff() warns about longitude/latitude or missing CRS", {
  sph <- make_sphere_fixtures()

  ## In the plane, longitude and latitude are read as planar coordinates, as
  ## sf and GEOS do.
  expect_warning(got <- dist_hausdorff(sph$mixed[1:3]),
                 "longitude and latitude coordinates")
  expect_equal(got,
               drop_u(sf::st_distance(sph$mixed[1:3], which = "Hausdorff")),
               tolerance = 1e-8)

  expect_warning(got_nc <- dist_hausdorff(sph$no_crs, use_sphere = TRUE),
                 "x has no CRS")
  expect_equal(got_nc,
               dist_hausdorff(sf::st_set_crs(sph$no_crs, "EPSG:4326"),
                              use_sphere = TRUE))
  expect_warning(
    expect_warning(dist_hausdorff(sph$no_crs, sph$no_crs, use_sphere = TRUE),
                   "y has no CRS"),
    "x has no CRS"
  )
})

## ---- Internal C++ functions ---------------------------------------------------

test_that("the internal C++ functions match dist_hausdorff() and close raw rings", {
  fix <- make_planar_fixtures()
  mats <- .to_list_mat(fix$polygons)
  pts_m <- .to_list_mat(fix$points)

  expect_equal(dist_haus_r2(mats), dist_hausdorff(fix$polygons))
  expect_equal(cdist_haus_r2(mats, pts_m),
               dist_hausdorff(fix$polygons, fix$points))

  ## A 2-column matrix without a part id is closed automatically.
  raw_open <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  raw_closed <- rbind(raw_open, c(0, 0))
  expect_equal(dist_haus_r2(list(raw_open, raw_closed), method = "boundary")[1L, 2L],
               0)

  skip_if_no_s2()
  sph <- make_sphere_fixtures()
  radius <- s2::s2_earth_radius_meters()
  pts <- sph$targets[c(2L, 4L)]
  mats_s <- .to_list_mat(sph$mixed)

  expect_equal(dist_haus_s2(mats_s, earth_radius = radius),
               dist_hausdorff(sph$mixed, use_sphere = TRUE))
  expect_equal(cdist_haus_s2(mats_s, .to_list_mat(pts), earth_radius = radius),
               dist_hausdorff(sph$mixed, pts, use_sphere = TRUE))
  expect_equal(dist_haus_s2(list(raw_open, raw_closed), method = "boundary",
                            earth_radius = radius)[1L, 2L],
               0)
})

test_that("the C++ functions check and convert their matrix input", {
  sq <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  pts_d <- rbind(c(2, 0), c(3, 0))
  pts_i <- matrix(c(2L, 3L, 0L, 0L), 2L)
  radius <- s2::s2_earth_radius_meters()

  ## Integer storage gives the same result as double storage.
  expect_equal(dist_haus_r2(list(sq, pts_i)), dist_haus_r2(list(sq, pts_d)))
  expect_equal(cdist_haus_r2(list(sq), list(pts_i)),
               cdist_haus_r2(list(sq), list(pts_d)))
  expect_equal(dist_haus_s2(list(sq, pts_i), earth_radius = radius),
               dist_haus_s2(list(sq, pts_d), earth_radius = radius))
  expect_equal(cdist_haus_s2(list(sq), list(pts_i), earth_radius = radius),
               cdist_haus_s2(list(sq), list(pts_d), earth_radius = radius))

  ## Invalid matrix inputs, methods, or tolerances stop with clear errors.
  expect_error(dist_haus_r2(list(sq, c(0, 0, 1))), "matrix")
  expect_error(dist_haus_s2(list(sq, c(0, 0, 1))), "matrix")
  expect_error(cdist_haus_r2(list(sq), list(matrix(c(5, 6), ncol = 1L))),
               "at least 2 columns")
  expect_error(cdist_haus_s2(list(sq), list(matrix(c(5, 6), ncol = 1L))),
               "at least 2 columns")
  expect_error(dist_haus_r2(list(sq), method = "bad"), "Unknown method")
  expect_error(dist_haus_s2(list(sq), method = "bad"), "Unknown method")
  expect_error(dist_haus_r2(list(sq), method = "boundary", tol = 0), "tol")
  expect_error(dist_haus_s2(list(sq), method = "boundary", tol = 0), "tol")
})

test_that(".to_list_mat() keeps each part and handles empty units", {
  fix <- make_planar_fixtures()

  mats <- .to_list_mat(fix$polygons)
  expect_length(mats, length(fix$polygons))
  ## p_hole has two rings, and p_multi has two polygons.
  expect_equal(length(unique(mats[[3L]][, 3L])), 2L)
  expect_equal(length(unique(mats[[4L]][, 3L])), 2L)

  all_empty <- sf::st_sfc(sf::st_polygon(), sf::st_point(), crs = "EPSG:3857")
  mats_e <- .to_list_mat(all_empty)
  expect_length(mats_e, 2L)
  expect_equal(dim(mats_e[[1L]]), c(0L, 3L))
  expect_length(.to_list_mat(sf::st_sfc(crs = "EPSG:3857")), 0L)

  mats_p <- .to_polygon_list(fix$mixed)
  expect_length(mats_p, length(fix$mixed))
  expect_equal(ncol(mats_p[[1L]]), 4L)
  ## Features 1:5 are polygons/multipolygons; 6:8 are points/linestring.
  for (k in 1:5) expect_true(all(mats_p[[k]][, 4L] == 1))
  for (k in 6:8) expect_true(all(mats_p[[k]][, 4L] == 0))
  ## Feature 9 is a GEOMETRYCOLLECTION of 1 point followed by 1 polygon.
  expect_equal(mats_p[[9L]][1L, 4L], 0)
  expect_true(all(mats_p[[9L]][-1L, 4L] == 1))

  mats_pe <- .to_polygon_list(all_empty)
  expect_length(mats_pe, 2L)
  expect_equal(dim(mats_pe[[1L]]), c(0L, 4L))
})

## ---- Planar continuous methods ("boundary" and "polygon") --------------------

test_that("planar continuous methods match exact analytical truth", {
  ## 1. Bridge multipolygon vs. spanning rectangle:
  ##    discrete = 1, boundary = 5, polygon = 5.
  rect <- sf::st_polygon(list(rbind(
    c(0, 0), c(10, 0), c(10, 1), c(0, 1), c(0, 0)
  )))
  bridge <- sf::st_multipolygon(list(
    list(rbind(c(-1, 0), c(0, 0), c(0, 1), c(-1, 1), c(-1, 0))),
    list(rbind(c(10, 0), c(11, 0), c(11, 1), c(10, 1), c(10, 0)))
  ))
  g_bridge <- sf::st_sfc(rect, bridge, crs = "EPSG:3857")
  expect_equal(dist_hausdorff(g_bridge, method = "discrete")[1L, 2L], 1)
  expect_equal(dist_hausdorff(g_bridge, method = "boundary")[1L, 2L], 5,
               tolerance = 1e-8)
  expect_equal(dist_hausdorff(g_bridge, method = "polygon")[1L, 2L], 5,
               tolerance = 1e-8)

  ## 2. Filled square vs. square with hole:
  ##    discrete = 1, boundary = 1, polygon = 4 (at interior center (5, 5)).
  sq10 <- sf::st_polygon(list(rbind(
    c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0)
  )))
  sq_hole <- sf::st_polygon(list(
    rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0)),
    rbind(c(1, 1), c(1, 9), c(9, 9), c(9, 1), c(1, 1))
  ))
  g_hole <- sf::st_sfc(sq10, sq_hole, crs = "EPSG:3857")
  expect_equal(dist_hausdorff(g_hole, method = "discrete")[1L, 2L], 1)
  expect_equal(dist_hausdorff(g_hole, method = "boundary")[1L, 2L], 1,
               tolerance = 1e-8)
  expect_equal(dist_hausdorff(g_hole, method = "polygon")[1L, 2L], 4,
               tolerance = 1e-7)

  ## 3. Square vs. concave X-shape touching all 4 corners:
  ##    discrete = 4.9, boundary = 4.9, polygon = 24.5 / sqrt(49.01) (< 4.9).
  x_shape <- sf::st_polygon(list(rbind(
    c(0, 0), c(5, 4.9), c(10, 0), c(5.1, 5),
    c(10, 10), c(5, 5.1), c(0, 10), c(4.9, 5), c(0, 0)
  )))
  g_x <- sf::st_sfc(sq10, x_shape, crs = "EPSG:3857")
  expect_equal(dist_hausdorff(g_x, method = "discrete")[1L, 2L], 4.9,
               tolerance = 1e-8)
  expect_equal(dist_hausdorff(g_x, method = "boundary")[1L, 2L], 4.9,
               tolerance = 1e-8)
  expect_equal(dist_hausdorff(g_x, method = "polygon")[1L, 2L],
               24.5 / sqrt(49.01), tolerance = 1e-7)

  ## 4. MULTIPOINT at corners of unit square vs. unit square:
  ##    discrete = 0, boundary = 0.5, polygon = sqrt(0.5).
  sq1 <- sf::st_polygon(list(rbind(
    c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
  )))
  corners <- sf::st_multipoint(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)))
  g_mp <- sf::st_sfc(sq1, corners, crs = "EPSG:3857")
  expect_equal(dist_hausdorff(g_mp, method = "discrete")[1L, 2L], 0)
  expect_equal(dist_hausdorff(g_mp, method = "boundary")[1L, 2L], 0.5,
               tolerance = 1e-8)
  expect_equal(dist_hausdorff(g_mp, method = "polygon")[1L, 2L], sqrt(0.5),
               tolerance = 1e-8)

  ## 5. Closed LINESTRING must NOT be treated as a filled polygon:
  ring_line <- sf::st_linestring(rbind(
    c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0)
  ))
  pt_center <- sf::st_point(c(5, 5))
  g_line <- sf::st_sfc(sq10, ring_line, pt_center, crs = "EPSG:3857")
  ## ring_line vs sq10: center (5, 5) is in sq10 and distance 5 from ring_line.
  expect_equal(dist_hausdorff(g_line, method = "discrete")[1L, 2L], 0)
  expect_equal(dist_hausdorff(g_line, method = "boundary")[1L, 2L], 0)
  expect_equal(dist_hausdorff(g_line, method = "polygon")[1L, 2L], 5,
               tolerance = 1e-7)

  ## 6. Identical polygons with collinear and shifted start vertices give 0.
  sq1_shifted <- sf::st_polygon(list(rbind(
    c(1, 0), c(1, 0.5), c(1, 1), c(0, 1), c(0, 0), c(0.5, 0), c(1, 0)
  )))
  g_eq <- sf::st_sfc(sq1, sq1_shifted, crs = "EPSG:3857")
  expect_equal(dist_hausdorff(g_eq, method = "boundary")[1L, 2L], 0,
               tolerance = 1e-8)
  expect_equal(dist_hausdorff(g_eq, method = "polygon")[1L, 2L], 0,
               tolerance = 1e-8)

  ## 7. Polygon with hole + island inside the hole vs. full square:
  ##    max distance occurs at the 4 diagonal equidistant points where
  ##    t - 1 = sqrt(2) * (4 - t), giving exact value 6 - 3 * sqrt(2) (~1.7574).
  sq_island <- sf::st_multipolygon(list(
    list(
      rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0)),
      rbind(c(1, 1), c(1, 9), c(9, 9), c(9, 1), c(1, 1))
    ),
    list(rbind(c(4, 4), c(6, 4), c(6, 6), c(4, 6), c(4, 4)))
  ))
  g_isl <- sf::st_sfc(sq10, sq_island, crs = "EPSG:3857")
  expect_equal(dist_hausdorff(g_isl, method = "polygon")[1L, 2L],
               6 - 3 * sqrt(2), tolerance = 1e-7)
})

test_that("all three methods agree on point sets and single point vs. polygon", {
  fix <- make_planar_fixtures()

  for (m in c("boundary", "polygon")) {
    expect_equal(dist_hausdorff(fix$points, method = m),
                 dist_hausdorff(fix$points, method = "discrete"),
                 tolerance = 1e-8, info = m)
    expect_equal(dist_hausdorff(fix$points, fix$polygons, method = m),
                 dist_hausdorff(fix$points, fix$polygons, method = "discrete"),
                 tolerance = 1e-8, info = m)
    expect_equal(dist_hausdorff(fix$polygons, fix$points, method = m),
                 dist_hausdorff(fix$polygons, fix$points, method = "discrete"),
                 tolerance = 1e-8, info = m)
  }

  ## Single point inside a polygon or inside a hole also agrees across methods.
  pt_in <- sf::st_sfc(sf::st_point(c(5, 5)), sf::st_point(c(1, 1)),
                      crs = "EPSG:3857")
  for (m in c("boundary", "polygon")) {
    expect_equal(dist_hausdorff(pt_in, fix$polygons, method = m),
                 dist_hausdorff(pt_in, fix$polygons, method = "discrete"),
                 tolerance = 1e-8, info = m)
  }
})

test_that("planar continuous methods match densified/grid sf references on fixtures", {
  fix <- make_planar_fixtures()

  for (nm in c("polygons", "targets", "mixed")) {
    geom <- fix[[nm]]
    got_b <- dist_hausdorff(geom, method = "boundary")
    ref_b <- haus_ref_boundary_r2(geom, step = 0.05)
    expect_lt(max(abs(got_b - ref_b)), 0.05)
    expect_true(all(got_b >= dist_hausdorff(geom, method = "discrete") - 1e-12),
                info = nm)
    expect_equal(diag(got_b), rep(0, length(geom)), info = nm)
    expect_equal(got_b, t(got_b), info = nm)

    got_p <- dist_hausdorff(geom, method = "polygon")
    ref_p <- haus_ref_polygon_r2(geom, step = 0.1)
    expect_lt(max(abs(got_p - ref_p)), 0.1)
    expect_equal(diag(got_p), rep(0, length(geom)), info = nm)
    expect_equal(got_p, t(got_p), info = nm)
  }

  ## Cross-distance symmetry bit-for-bit between (x, y) and t(y, x).
  for (m in c("boundary", "polygon")) {
    ab <- dist_hausdorff(fix$polygons, fix$targets, method = m)
    ba <- dist_hausdorff(fix$targets, fix$polygons, method = m)
    expect_equal(ab, t(ba), tolerance = 0, info = m)
  }

  ## Large projected coordinates (EPSG:3857 scale) stay accurate.
  for (m in c("boundary", "polygon")) {
    expect_equal(dist_hausdorff(fix$large_proj, method = m),
                 dist_hausdorff(fix$large_proj, method = "discrete"),
                 tolerance = 1e-5, info = m)
  }
})

test_that("planar continuous methods match densified references on North Carolina", {
  nc_proj <- sf::st_transform(read_nc()[1:10], 32119)

  got_d <- dist_hausdorff(nc_proj, method = "discrete")
  got_b <- dist_hausdorff(nc_proj, method = "boundary")
  got_p <- dist_hausdorff(nc_proj, method = "polygon")

  ref_b <- haus_ref_boundary_r2(nc_proj, step = 250)
  ref_p <- haus_ref_polygon_r2(nc_proj, step = 1000)

  expect_lt(max(abs(got_b - ref_b)), 250)
  expect_lt(max(abs(got_p - ref_p)), 1000)
  expect_true(all(got_b >= got_d - 1e-8))
  expect_equal(diag(got_b), rep(0, 10L))
  expect_equal(diag(got_p), rep(0, 10L))
  expect_equal(got_b, t(got_b), tolerance = 0)
  expect_equal(got_p, t(got_p), tolerance = 0)
  expect_equal(diag(dist_hausdorff(nc_proj[1:3], nc_proj[1:3], method = "polygon")),
               rep(0, 3L), tolerance = 0)

  ## Triangle inequality across sampled triples.
  for (m_mat in list(got_b, got_p)) {
    for (i in 1:5) {
      for (j in 6:8) {
        for (k in 9:10) {
          expect_true(m_mat[i, k] <= m_mat[i, j] + m_mat[j, k] + 1e-6)
        }
      }
    }
  }
})

test_that("continuous methods handle empty, non-finite, zero-length, and invalid args", {
  fix <- make_planar_fixtures()
  zero <- sf::st_sfc(crs = "EPSG:3857")
  wkt <- c("MULTIPOINT ((0 0), EMPTY, (10 0))", "POINT (50 0)")
  g_nan <- sf::st_as_sfc(wkt, crs = "EPSG:3857")

  for (m in c("boundary", "polygon")) {
    got_e <- dist_hausdorff(fix$with_empty, method = m)
    expect_true(all(is.na(got_e[2L, ])))
    expect_true(all(is.na(got_e[, 2L])))
    expect_false(anyNA(got_e[-2L, -2L]))

    got_n <- dist_hausdorff(g_nan, method = m)
    expect_true(all(is.na(got_n[1L, ])))
    expect_true(all(is.na(got_n[, 1L])))
    expect_false(is.na(got_n[2L, 2L]))

    expect_equal(dim(dist_hausdorff(zero, method = m)), c(0L, 0L))
    expect_equal(dim(dist_hausdorff(zero, fix$targets, method = m)), c(0L, 3L))
  }

  expect_error(dist_hausdorff(fix$polygons, method = "unknown"))
  expect_error(dist_hausdorff(fix$polygons, method = "boundary", tol = 0), "tol")
  expect_error(dist_hausdorff(fix$polygons, method = "polygon", tol = -1), "tol")
  expect_error(dist_hausdorff(fix$polygons, method = "boundary", tol = NA), "tol")
})

## ---- Spherical continuous methods ("boundary" and "polygon") ----------------

test_that("spherical continuous methods match exact analytical truth on S^2", {
  skip_if_no_s2()
  R <- s2::s2_earth_radius_meters()

  ## 1. Great-circle arch test: distance to the equator along an arc from
  ##    (0, 45) to (90, 45) bulges to latitude atan(sqrt(2)) (~54.7356 deg),
  ##    exceeding both endpoints' 45-deg distance. This checks that S^2 does
  ##    NOT wrongly assume line-arc convexity.
  eq <- sf::st_linestring(rbind(c(0, 0), c(90, 0)))
  arc45 <- sf::st_linestring(rbind(c(0, 45), c(90, 45)))
  g_arch <- sf::st_sfc(eq, arc45, crs = "EPSG:4326")
  expect_equal(dist_hausdorff(g_arch, use_sphere = TRUE, method = "discrete")[1L, 2L],
               (pi / 4) * R, tolerance = 1e-8)
  expect_equal(dist_hausdorff(g_arch, use_sphere = TRUE, method = "boundary")[1L, 2L],
               atan(sqrt(2)) * R, tolerance = 1e-7)
  expect_equal(dist_hausdorff(g_arch, use_sphere = TRUE, method = "polygon")[1L, 2L],
               atan(sqrt(2)) * R, tolerance = 1e-7)

  ## 2. Equator bridge multipolygon vs. spanning rectangle:
  ##    discrete = 1 deg * R, boundary = polygon = 5 deg * R.
  rect_eq <- sf::st_polygon(list(rbind(
    c(0, 0), c(10, 0), c(10, 1), c(0, 1), c(0, 0)
  )))
  bridge_eq <- sf::st_multipolygon(list(
    list(rbind(c(-1, 0), c(0, 0), c(0, 1), c(-1, 1), c(-1, 0))),
    list(rbind(c(10, 0), c(11, 0), c(11, 1), c(10, 1), c(10, 0)))
  ))
  g_bridge <- sf::st_sfc(rect_eq, bridge_eq, crs = "EPSG:4326")
  deg_m <- (pi / 180) * R
  expect_equal(dist_hausdorff(g_bridge, use_sphere = TRUE, method = "discrete")[1L, 2L],
               1 * deg_m, tolerance = 1e-7)
  expect_equal(dist_hausdorff(g_bridge, use_sphere = TRUE, method = "boundary")[1L, 2L],
               5 * deg_m, tolerance = 1e-7)
  expect_equal(dist_hausdorff(g_bridge, use_sphere = TRUE, method = "polygon",
                              tol = 1e-6)[1L, 2L],
               5 * deg_m, tolerance = 1e-6)

  ## 3. Spherical square vs. square with hole centered on equator:
  ##    discrete and boundary ~ 1 deg * R; polygon = 4 deg * R at center (0, 0).
  sq_eq <- sf::st_polygon(list(rbind(
    c(-5, -5), c(5, -5), c(5, 5), c(-5, 5), c(-5, -5)
  )))
  sq_eq_hole <- sf::st_polygon(list(
    rbind(c(-5, -5), c(5, -5), c(5, 5), c(-5, 5), c(-5, -5)),
    rbind(c(-4, -4), c(-4, 4), c(4, 4), c(4, -4), c(-4, -4))
  ))
  g_hole <- sf::st_sfc(sq_eq, sq_eq_hole, crs = "EPSG:4326")
  expect_lt(dist_hausdorff(g_hole, use_sphere = TRUE, method = "boundary")[1L, 2L],
            1.1 * deg_m)
  expect_equal(dist_hausdorff(g_hole, use_sphere = TRUE, method = "polygon",
                              tol = 1e-6)[1L, 2L],
               4 * deg_m, tolerance = 1e-6)

  ## 4. Antipodal POINT vs. POLYGON on S^2:
  ##    Antipode (180, 0) of POINT (0, 0) lies inside box [178, 182] x [-1, 1],
  ##    so polygon method gives exact pi * R (> discrete).
  pt0 <- sf::st_point(c(0, 0))
  anti_box <- sf::st_polygon(list(rbind(
    c(178, -1), c(-178, -1), c(-178, 1), c(178, 1), c(178, -1)
  )))
  g_anti <- sf::st_sfc(pt0, anti_box, crs = "EPSG:4326")
  d_anti_d <- dist_hausdorff(g_anti, use_sphere = TRUE, method = "discrete")[1L, 2L]
  d_anti_b <- dist_hausdorff(g_anti, use_sphere = TRUE, method = "boundary")[1L, 2L]
  d_anti_p <- dist_hausdorff(g_anti, use_sphere = TRUE, method = "polygon")[1L, 2L]
  expect_lt(d_anti_d, pi * R - 1000)
  expect_lt(d_anti_b, pi * R - 1000)
  expect_equal(d_anti_p, pi * R, tolerance = 1e-7)

  ## 5. Clockwise vs. counter-clockwise ring equivalence (matches sf / s2).
  ccw <- sf::st_polygon(list(rbind(
    c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0)
  )))
  cw <- sf::st_polygon(list(rbind(
    c(0, 0), c(0, 10), c(10, 10), c(10, 0), c(0, 0)
  )))
  g_cw <- sf::st_sfc(ccw, cw, sf::st_point(c(5, 5)), crs = "EPSG:4326")
  m_cw <- dist_hausdorff(g_cw, use_sphere = TRUE, method = "polygon")
  expect_equal(m_cw[1L, 2L], 0, tolerance = 1e-8)
  expect_equal(m_cw[1L, 3L], m_cw[2L, 3L], tolerance = 1e-8)

  ## 6. North Pole cap vs. North Pole point (no pole singularity).
  pole_cap <- sf::st_polygon(list(rbind(
    c(0, 80), c(90, 80), c(180, 80), c(-90, 80), c(0, 80)
  )))
  g_pole <- sf::st_sfc(pole_cap, sf::st_point(c(0, 90)), crs = "EPSG:4326")
  expect_equal(dist_hausdorff(g_pole, use_sphere = TRUE, method = "polygon")[1L, 2L],
               10 * deg_m, tolerance = 1e-7)
})

test_that("spherical continuous methods warn on cap-limit and slow-convergence cases", {
  skip_if_no_s2()

  ## 1. Global-scale polygon exceeding the 0.45 * pi cap limit returns NA with
  ##    a warning for method = "polygon", while "discrete" and "boundary" run.
  global_mp <- sf::st_multipolygon(list(
    list(rbind(c(0, 0), c(5, 0), c(5, 5), c(0, 5), c(0, 0))),
    list(rbind(c(175, -5), c(180, -5), c(180, 0), c(175, 0), c(175, -5)))
  ))
  small_sq <- sf::st_polygon(list(rbind(
    c(20, 0), c(25, 0), c(25, 5), c(20, 5), c(20, 0)
  )))
  g_cap <- sf::st_sfc(global_mp, small_sq, crs = "EPSG:4326")

  expect_warning(
    got_p <- dist_hausdorff(g_cap, use_sphere = TRUE, method = "polygon"),
    "spherical cap limit"
  )
  expect_true(all(is.na(got_p[1L, ])))
  expect_true(all(is.na(got_p[, 1L])))
  expect_equal(got_p[2L, 2L], 0)

  expect_warning(
    got_xy <- dist_hausdorff(g_cap[2L], g_cap[1L], use_sphere = TRUE, method = "polygon"),
    "spherical cap limit"
  )
  expect_true(is.na(got_xy[1L, 1L]))

  expect_no_warning(got_d <- dist_hausdorff(g_cap, use_sphere = TRUE, method = "discrete"))
  expect_no_warning(got_b <- dist_hausdorff(g_cap, use_sphere = TRUE, method = "boundary"))
  expect_false(anyNA(got_d))
  expect_false(anyNA(got_b))

  ## 2. Identical polygons with mismatched edge discretization hit the bounded
  ##    iteration budget cleanly and warn instead of stalling.
  pole_cap <- sf::st_polygon(list(rbind(
    c(0, 80), c(90, 80), c(180, 80), c(-90, 80), c(0, 80)
  )))
  u0 <- .ll_to_xyz_r(0, 80)
  u1 <- .ll_to_xyz_r(90, 80)
  mid_ll <- .xyz_to_ll_r((u0 + u1) / sqrt(sum((u0 + u1)^2)))
  pole_mid <- sf::st_polygon(list(rbind(
    c(0, 80), mid_ll, c(90, 80), c(180, 80), c(-90, 80), c(0, 80)
  )))
  g_slow <- sf::st_sfc(pole_cap, pole_mid, crs = "EPSG:4326")
  expect_warning(
    got_s <- dist_hausdorff(g_slow, use_sphere = TRUE, method = "polygon"),
    "iteration limit"
  )
  expect_true(is.finite(got_s[1L, 2L]))
  expect_lt(got_s[1L, 2L], 1e-3)
})

test_that("spherical continuous methods match s2 references on fixtures and NC", {
  skip_if_no_s2()
  sph <- make_sphere_fixtures()

  for (m in c("boundary", "polygon")) {
    got_m <- dist_hausdorff(sph$mixed, use_sphere = TRUE, method = m)
    expect_equal(diag(got_m), rep(0, length(sph$mixed)))
    expect_equal(got_m, t(got_m), tolerance = 0)

    ab <- dist_hausdorff(sph$mixed, sph$targets, use_sphere = TRUE, method = m)
    ba <- dist_hausdorff(sph$targets, sph$mixed, use_sphere = TRUE, method = m)
    expect_equal(ab, t(ba), tolerance = 0, info = m)

    ## Empty units and non-finite coordinates return NA.
    got_e <- dist_hausdorff(sph$with_empty, use_sphere = TRUE, method = m)
    expect_true(all(is.na(got_e[2L, ])))
    expect_true(all(is.na(got_e[, 2L])))
    expect_false(anyNA(got_e[-2L, -2L]))
  }

  ## Cross-check North Carolina counties against densified/grid s2 references.
  nc_ll <- sf::st_transform(read_nc()[1:8], "EPSG:4326")
  got_d <- dist_hausdorff(nc_ll, use_sphere = TRUE, method = "discrete")
  got_b <- dist_hausdorff(nc_ll, use_sphere = TRUE, method = "boundary")
  got_p <- dist_hausdorff(nc_ll, use_sphere = TRUE, method = "polygon")

  ref_b <- haus_ref_boundary_s2(nc_ll, max_deg = 0.005)
  ref_p <- haus_ref_polygon_s2(nc_ll, max_deg = 0.01, grid_deg = 0.015)

  expect_lt(max(abs(got_b - ref_b)), 600)
  expect_lt(max(abs(got_p - ref_p)), 1500)
  expect_true(all(got_b >= got_d - 1e-6))
  expect_equal(diag(got_b), rep(0, 8L))
  expect_equal(diag(got_p), rep(0, 8L))
  expect_equal(got_b, t(got_b), tolerance = 0)
  expect_equal(got_p, t(got_p), tolerance = 0)

  ## Point sets agree across all three methods on the sphere.
  pts <- sf::st_sfc(
    sf::st_point(c(-46.63, -23.55)),
    sf::st_point(c(-74.00, 40.71)),
    sf::st_point(c(139.69, 35.69)),
    crs = "EPSG:4326"
  )
  for (m in c("boundary", "polygon")) {
    expect_equal(dist_hausdorff(pts, use_sphere = TRUE, method = m),
                 drop_u(sf::st_distance(pts)), tolerance = 1e-10, info = m)
  }

  ## Triangle inequality on spherical continuous matrices.
  for (m_mat in list(got_b, got_p)) {
    for (i in 1:4) {
      for (j in 5:6) {
        for (k in 7:8) {
          expect_true(m_mat[i, k] <= m_mat[i, j] + m_mat[j, k] + 1e-4)
        }
      }
    }
  }
})
