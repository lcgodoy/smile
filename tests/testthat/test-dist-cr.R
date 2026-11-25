## ---- Planar centers and radii: compare with sf ---------------------------------

test_that("planar get_mbc_params() matches sf::st_minimum_bounding_circle()", {
  skip_if_no_mbc()
  fix <- make_planar_fixtures()

  for (nm in c("mixed", "polygons", "targets", "points", "large_proj")) {
    geom <- fix[[nm]]
    got <- get_mbc_params(geom)
    ref <- mbc_ref_planar(geom)

    expect_equal(got$radius, ref$radius, tolerance = 1e-8, info = nm)
    expect_equal(unname(sf::st_coordinates(got$center)),
                 unname(sf::st_coordinates(ref$center)),
                 tolerance = 1e-8, info = nm)
    expect_equal(sf::st_crs(got$center), sf::st_crs(geom), info = nm)
  }

  nc_proj <- sf::st_transform(read_nc()[1:25], 32119)
  got <- get_mbc_params(nc_proj)
  ref <- mbc_ref_planar(nc_proj)
  expect_equal(got$radius, ref$radius, tolerance = 1e-8)
  expect_equal(unname(sf::st_coordinates(got$center)),
               unname(sf::st_coordinates(ref$center)),
               tolerance = 1e-8)
})

test_that("each planar circle holds all vertices of its unit", {
  fix <- make_planar_fixtures()
  geom <- fix$mixed
  got <- get_mbc_params(geom)
  cxy <- sf::st_coordinates(got$center)

  for (i in seq_along(geom)) {
    pts <- wk::wk_coords(geom[i])
    d <- sqrt((pts$x - cxy[i, "X"])^2 + (pts$y - cxy[i, "Y"])^2)
    expect_true(all(d <= got$radius[i] * (1 + 1e-12) + 1e-12), info = i)
    expect_equal(max(d), got$radius[i], tolerance = 1e-10, info = i)
  }
})

## ---- Planar CR distance: compare with sf ---------------------------------------

test_that("planar dist_cr() matches the sf reference", {
  skip_if_no_mbc()
  fix <- make_planar_fixtures()

  for (nm in c("mixed", "polygons", "targets", "points", "large_proj")) {
    geom <- fix[[nm]]
    got <- dist_cr(geom)
    expect_equal(got, cr_ref(geom), tolerance = 1e-8, info = nm)
    expect_equal(diag(got), rep(0, length(geom)), info = nm)
    expect_equal(got, t(got), info = nm)
  }

  expect_equal(dist_cr(fix$polygons, fix$targets),
               cr_ref(fix$polygons, fix$targets), tolerance = 1e-8)

  nc_proj <- sf::st_transform(read_nc()[1:25], 32119)
  expect_equal(dist_cr(nc_proj), cr_ref(nc_proj), tolerance = 1e-8)
  expect_equal(dist_cr(nc_proj[1:12], nc_proj[13:25]),
               cr_ref(nc_proj[1:12], nc_proj[13:25]), tolerance = 1e-8)
})

## ---- Spherical centers and radii ---------------------------------------------

test_that("spherical get_mbc_params() matches the pure-R Welzl reference", {
  sph <- make_sphere_fixtures()

  got <- get_mbc_params(sph$mixed, use_sphere = TRUE)
  ref <- mbc_ref_welzl_s2(sph$mixed)
  expect_equal(got$radius, ref$radius, tolerance = 1e-6)
  expect_equal(sf::st_coordinates(got$center), sf::st_coordinates(ref$center),
               tolerance = 1e-9)
  expect_equal(sf::st_crs(got$center), sf::st_crs("EPSG:4326"))

  ## A cap near the pole uses the other tangent-plane axis in cap3().
  polar <- sf::st_sfc(
    sf::st_multipoint(rbind(c(-10, 89.2), c(110, 89.2), c(230, 89.2))),
    crs = "EPSG:4326"
  )
  got_p <- get_mbc_params(polar, use_sphere = TRUE)
  ref_p <- mbc_ref_welzl_s2(polar)
  expect_equal(got_p$radius, ref_p$radius, tolerance = 1e-6)
  expect_equal(sf::st_coordinates(got_p$center),
               sf::st_coordinates(ref_p$center), tolerance = 1e-9)

  nc_ll <- sf::st_transform(read_nc()[1:25], "EPSG:4326")
  got_nc <- get_mbc_params(nc_ll, use_sphere = TRUE)
  ref_nc <- mbc_ref_welzl_s2(nc_ll)
  expect_equal(got_nc$radius, ref_nc$radius, tolerance = 1e-6)
  expect_equal(sf::st_coordinates(got_nc$center),
               sf::st_coordinates(ref_nc$center), tolerance = 1e-9)

  ## Projected input is transformed to EPSG:4326.
  fix <- make_planar_fixtures()
  got_tr <- get_mbc_params(fix$large_proj, use_sphere = TRUE)
  ref_tr <- mbc_ref_welzl_s2(fix$large_proj)
  expect_equal(got_tr$radius, ref_tr$radius, tolerance = 1e-6)
  expect_equal(sf::st_coordinates(got_tr$center),
               sf::st_coordinates(ref_tr$center), tolerance = 1e-9)
  expect_equal(sf::st_crs(got_tr$center), sf::st_crs("EPSG:4326"))

  ## Longitude and latitude input keeps its CRS (NAD27 for nc.shp).
  nc <- read_nc()[1:3]
  expect_equal(sf::st_crs(get_mbc_params(nc, use_sphere = TRUE)$center),
               sf::st_crs(nc))
})

test_that("each spherical cap holds all vertices of its unit (sf/s2 check)", {
  skip_if_no_s2()
  sph <- make_sphere_fixtures()
  nc_ll <- sf::st_transform(read_nc()[1:10], "EPSG:4326")

  for (geom in list(sph$mixed, nc_ll)) {
    got <- get_mbc_params(geom, use_sphere = TRUE)
    for (i in seq_along(geom)) {
      d <- drop_u(sf::st_distance(got$center[i], vertices_sfc(geom[i])))
      expect_true(all(d <= got$radius[i] * (1 + 1e-9)), info = i)
      expect_equal(max(d), got$radius[i], tolerance = 1e-7, info = i)
    }
  }
})

## ---- Spherical CR distance: compare with sf/s2 -------------------------------

test_that("spherical dist_cr() matches the sf/s2 reference", {
  skip_if_no_s2()
  sph <- make_sphere_fixtures()

  got <- dist_cr(sph$mixed, use_sphere = TRUE)
  expect_equal(got, cr_ref(sph$mixed, use_sphere = TRUE), tolerance = 1e-6)
  expect_equal(diag(got), rep(0, length(sph$mixed)))
  expect_equal(got, t(got))

  expect_equal(dist_cr(sph$mixed, sph$targets, use_sphere = TRUE),
               cr_ref(sph$mixed, sph$targets, use_sphere = TRUE),
               tolerance = 1e-6)

  nc_ll <- sf::st_transform(read_nc()[1:20], "EPSG:4326")
  expect_equal(dist_cr(nc_ll, use_sphere = TRUE),
               cr_ref(nc_ll, use_sphere = TRUE), tolerance = 1e-6)
  expect_equal(dist_cr(nc_ll[1:8], nc_ll[9:20], use_sphere = TRUE),
               cr_ref(nc_ll[1:8], nc_ll[9:20], use_sphere = TRUE),
               tolerance = 1e-6)
})

test_that("spherical dist_cr() between points is the great-circle distance", {
  skip_if_no_s2()
  ## A point has radius 0, so the CR distance is the distance between the
  ## points. A tight tolerance also checks the radius of the sphere.
  pts <- sf::st_sfc(
    sf::st_point(c(-46.63, -23.55)),
    sf::st_point(c(-74.00, 40.71)),
    sf::st_point(c(139.69, 35.69)),
    crs = "EPSG:4326"
  )
  expect_equal(dist_cr(pts, use_sphere = TRUE),
               drop_u(sf::st_distance(pts)), tolerance = 1e-10)
  expect_equal(dist_cr(pts[1L], pts[2:3], use_sphere = TRUE),
               drop_u(sf::st_distance(pts[1L], pts[2:3])), tolerance = 1e-10)
})

test_that("spherical dist_cr() transforms projected input", {
  nc_proj <- sf::st_transform(read_nc()[1:6], 32119)
  nc_ll <- sf::st_transform(nc_proj, "EPSG:4326")

  expect_equal(dist_cr(nc_proj, use_sphere = TRUE),
               dist_cr(nc_ll, use_sphere = TRUE))
  expect_equal(dist_cr(nc_proj[1:2], nc_proj[3:6], use_sphere = TRUE),
               dist_cr(nc_ll[1:2], nc_ll[3:6], use_sphere = TRUE))
})

## ---- Empty units, zero-length input, CRS, and warnings ----------------------

test_that("planar CR functions keep the NA rules for empty units", {
  fix <- make_planar_fixtures()

  got <- dist_cr(fix$with_empty)
  expect_equal(diag(got), c(0, NA_real_, 0))
  expect_true(all(is.na(got[2L, ])))
  expect_true(all(is.na(got[, 2L])))
  expect_false(is.na(got[1L, 3L]))

  got_x <- dist_cr(fix$with_empty, fix$targets)
  expect_true(all(is.na(got_x[2L, ])))
  expect_false(anyNA(got_x[c(1L, 3L), ]))

  mbc <- get_mbc_params(fix$with_empty)
  expect_true(is.na(mbc$radius[2L]))
  expect_true(sf::st_is_empty(mbc$center[2L]))
  expect_equal(mbc$radius[3L], 0)

  skip_if_no_mbc()
  expect_equal(got, cr_ref(fix$with_empty), tolerance = 1e-8)
})

test_that("spherical CR functions keep the NA rules for empty units", {
  skip_if_no_s2()
  sph <- make_sphere_fixtures()

  expect_no_warning({
    got_cr <- dist_cr(sph$with_empty, use_sphere = TRUE)
    got_mbc <- get_mbc_params(sph$with_empty, use_sphere = TRUE)
  })
  expect_equal(got_cr, cr_ref(sph$with_empty, use_sphere = TRUE),
               tolerance = 1e-6)
  expect_true(all(is.na(got_cr[2L, ])))
  expect_true(all(is.na(got_cr[, 2L])))
  expect_true(is.na(got_mbc$radius[2L]))
})

test_that("CR functions handle zero-length input and CRS mismatch", {
  fix <- make_planar_fixtures()
  zero <- sf::st_sfc(crs = "EPSG:3857")

  expect_equal(dim(dist_cr(zero)), c(0L, 0L))
  expect_equal(dim(dist_cr(zero, fix$targets)), c(0L, 3L))
  expect_equal(dim(dist_cr(zero, use_sphere = TRUE)), c(0L, 0L))
  expect_length(get_mbc_params(zero)$center, 0L)
  expect_length(get_mbc_params(zero)$radius, 0L)

  expect_error(dist_cr(fix$polygons, sf::st_transform(fix$targets, 4326)),
               "same CRS")
})

test_that("CR functions warn about units outside an open hemisphere", {
  sph <- make_sphere_fixtures()

  expect_warning(got_mbc <- get_mbc_params(sph$over_hemi, use_sphere = TRUE),
                 "not inside an open hemisphere")
  expect_true(is.na(got_mbc$radius[2L]))
  expect_true(sf::st_is_empty(got_mbc$center[2L]))

  expect_warning(dist_cr(sph$over_hemi, use_sphere = TRUE),
                 "not inside an open hemisphere")
  expect_warning(dist_cr(sph$over_hemi, sph$targets, use_sphere = TRUE),
                 "not inside an open hemisphere")
  expect_warning(dist_cr(sph$targets, sph$over_hemi, use_sphere = TRUE),
                 "not inside an open hemisphere")

  ## One invalid unit only: the 1 x 1 matrix is NA, and the warning occurs.
  expect_warning(got_1 <- dist_cr(sph$over_hemi[2L], use_sphere = TRUE),
                 "not inside an open hemisphere")
  expect_true(is.na(got_1[1L, 1L]))

  ## Zero-length x: the matrix is empty, but the warning about y occurs.
  expect_warning(dist_cr(sph$over_hemi[0L], sph$over_hemi, use_sphere = TRUE),
                 "units of y are not inside an open hemisphere: 2")

  ## The warning names the argument that holds the unit.
  expect_warning(dist_cr(sph$targets, sph$over_hemi, use_sphere = TRUE),
                 "units of y are not inside an open hemisphere: 2")
  expect_warning(dist_cr(sph$over_hemi, sph$targets, use_sphere = TRUE),
                 "units of x are not inside an open hemisphere: 2")

  skip_if_no_s2()
  expect_equal(suppressWarnings(dist_cr(sph$over_hemi, use_sphere = TRUE)),
               cr_ref(sph$over_hemi, use_sphere = TRUE), tolerance = 1e-6)
})

test_that("CR functions warn about longitude/latitude or missing CRS", {
  sph <- make_sphere_fixtures()

  expect_warning(got_ll <- dist_cr(sph$mixed[1:3]),
                 "longitude and latitude coordinates")
  expect_warning(get_mbc_params(sph$mixed[1:3]),
                 "longitude and latitude coordinates")

  expect_warning(dist_cr(sph$no_crs, use_sphere = TRUE), "x has no CRS")
  expect_warning(get_mbc_params(sph$no_crs, use_sphere = TRUE), "x has no CRS")

  skip_if_no_mbc()
  expect_equal(got_ll, cr_ref(sph$mixed[1:3]), tolerance = 1e-8)
})
