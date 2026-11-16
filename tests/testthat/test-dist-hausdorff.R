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

test_that("the point and single-pair C++ functions agree with the list functions", {
  fix <- make_planar_fixtures()
  mats <- .to_list_mat(fix$polygons)

  expect_equal(cdist_haus_r2_mp(mats, sf::st_coordinates(fix$points)),
               dist_hausdorff(fix$polygons, fix$points))
  expect_equal(haus_r2_aux_continuous(mats[[1L]], mats[[2L]]),
               dist_hausdorff(fix$polygons)[1L, 2L])

  sph <- make_sphere_fixtures()
  radius <- s2::s2_earth_radius_meters()
  pts <- sph$targets[c(2L, 4L)]
  mats_s <- .to_list_mat(sph$mixed)

  expect_equal(cdist_haus_s2_mp(mats_s, sf::st_coordinates(pts),
                                earth_radius = radius),
               dist_hausdorff(sph$mixed, pts, use_sphere = TRUE))
  expect_equal(haus_s2_aux_continuous(mats_s[[1L]], mats_s[[3L]],
                                      earth_radius = radius),
               dist_hausdorff(sph$mixed, use_sphere = TRUE)[1L, 3L])
})

test_that("the C++ functions check and convert their matrix input", {
  sq <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  pts_d <- rbind(c(2, 0), c(3, 0))
  pts_i <- matrix(c(2L, 3L, 0L, 0L), 2L)
  radius <- s2::s2_earth_radius_meters()

  ## Integer storage gives the same result as double storage.
  expect_equal(cdist_haus_r2_mp(list(sq), pts_i),
               cdist_haus_r2_mp(list(sq), pts_d))
  expect_equal(cdist_haus_s2_mp(list(sq), pts_i, earth_radius = radius),
               cdist_haus_s2_mp(list(sq), pts_d, earth_radius = radius))
  expect_equal(dist_haus_r2(list(sq, pts_i)), dist_haus_r2(list(sq, pts_d)))

  ## A vector, or a matrix with one column, gives an error.
  expect_error(haus_r2_aux_continuous(c(0, 0), c(3, 4)), "matrix")
  expect_error(haus_s2_aux_continuous(c(0, 0), c(3, 4)), "matrix")
  expect_error(cdist_haus_r2_mp(list(sq), matrix(c(5, 6), ncol = 1L)),
               "at least 2 columns")
  expect_error(cdist_haus_s2_mp(list(sq), matrix(c(5, 6), ncol = 1L)),
               "at least 2 columns")
  expect_error(dist_haus_r2(list(sq, c(0, 0, 1))), "matrix")
  expect_error(dist_haus_s2(list(sq, c(0, 0, 1))), "matrix")
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
})
