##' @name weight_mat
w_col <- function(source_unit, target) {
    source_unit <- sf::st_sfc(source_unit,
                              crs = sf::st_crs(target))
    ints <- sf::st_intersects(target, source_unit,
                          sparse = FALSE)
    out <- vector(mode = "numeric",
                  length = length(target))
    out[ints] <- as.numeric(
        sf::st_area(sf::st_intersection(target[ints], source_unit))
    )
    return(out)
}

##' @name weight_mat
##'
##' @title Building weight matrix \strong{W} for Areal Interpolation
##' 
##' @description internal use. \eqn{W_{ij} = | A_i \, cap \, B_j |}.
##' @param source a \code{sf} object - source spatial data.
##' @param source_unit a single \code{geometry} from the source dataset.
##' @param source_dt a \code{data.frame} object representing the source dataset
##'     but excludying the \code{geometry}, i.e. the spatial information,
##'     column.
##' @param target a \code{sf} object - target spatial data.
##' @param var_vec a \code{numeric} vector with variances observed at the source
##'     data.
##' @param W the weight matrix.
##' @param method a \code{character} representing the method to approximate the
##'     variance of the AI estimates. Possible values are "CS"
##'     (Cauchy-Schwartz) or "MI" (Moran's I).
##' @param rho_mi \code{numeric} calcuated Moran's I.
##' @return A \eqn{n \times m} \code{numeric} matrix. Where \eqn{n} is the
##'     number of objservations in the target and \eqn{m} is the sample size in
##'     the source dataset.
##' @keywords internal
build_w <- function(source, target) {
  source <- sf::st_geometry(source)
  target <- sf::st_geometry(target)
  proj_source <- sf::st_crs(source)
  proj_target <- sf::st_crs(target)
  if( proj_source != proj_target ) {
    warning("reprojecting target.")
    target <- sf::st_transform(target, proj_source)
  }
  w_aux <- lapply(source, w_col,
                  target = target)
  w_aux <- do.call("cbind", w_aux)
  div_target <- (1 / as.numeric(sf::st_area(target))) |>
      diag()

  return(div_target %*% w_aux)
}

##' @name weight_mat
est_w <- function(W, source_dt, target) {
    estimates <- (W %*% as.matrix(source_dt)) |>
        as.data.frame()
    cbind(target, estimates)
}

##' @title Calculates the (global) Moran's I
##' @param sf_dt a \code{sf} (with POLYGON \code{geometry}) dataset.
##' @param variable a \code{character} representing one of the variables from
##'     \code{sf_dt}.
##' @return a \code{numeric} scalar.
##' @keywords internal
morans_i <- function(sf_dt, variable) {
    stopifnot(inherits(sf_dt, "sf"))
    stopifnot(inherits(variable, "character"))
    stopifnot(length(variable) == 1)

    adj_mat <- sf::st_intersects(x = sf::st_geometry(sf_dt),
                        sparse = FALSE) |>
        as.numeric() |>
        matrix(ncol = nrow(sf_dt))
    
    diag(adj_mat) <- 0
    
    y <- scale(sf_dt[[variable]])
    .n <- length(y)
    ones <- matrix(rep(1, .n), ncol = 1)

    val1 <- as.numeric(crossprod(ones))
    val2 <- tcrossprod(ones)
    val3 <- diag(nrow = .n, ncol = .n)
    val4 <- val3 - (val2 / val1)
    val4 <- val4 %*% y

    as.numeric(
        .n / (crossprod(ones, adj_mat) %*% ones) *
        ((crossprod(val4, adj_mat) %*% val4) /
         crossprod(val4))
    )
}

##' @name weight_mat
var_w <- function(W, var_vec, target,
                  method = "CS", rho_mi) {
    if(method == "MI")
        stopifnot(!missing(rho_mi))

    stopifnot(NCOL(var_vec) == 1)
    
    .m <- length(var_vec)
    
    cov_mat <- tcrossprod(sqrt(var_vec))

    if(method == "MI") {
        mat_aux <- rho_mi * (tcrossprod(rep(1, .m)) - diag(.m))
        cov_mat <- cov_mat * mat_aux
    }
    
    se_est <- (W %*% tcrossprod(cov_mat, W)) |>
        diag() |>
        sqrt()

    return(transform(target, se_est = se_est))
}

##' @name AI
##'
##' @title Areal Interpolation
##' @param source a \code{sf} object - source spatial data.
##' @param target a \code{sf} object - target spatial data.
##' @param vars a \code{character} representing the variables (observed at the
##'     source) to be estimated at the target data.
##' @param vars_var a scalar of type \code{character} representing the name of
##'     the variable in the source dataset that stores the variances of the
##'     variable to be estimated at the target data.
##' @param var_method a \code{character} representing the method to approximate
##'     the variance of the AI estimates. Possible values are "CS"
##'     (Cauchy-Schwartz) or "MI" (Moran's I).
##' 
##' @return the target (of type \code{sf}) with estimates of the variables
##'     observed at the source data.
##' @export
ai <- function(source, target,
               vars) {
    stopifnot(inherits(source, "sf"))
    stopifnot(inherits(target, "sf"))
    stopifnot(inherits(vars, "character"))
    source_dt <- sf::st_drop_geometry(source)
    stopifnot(all(!is.na(c(source_dt[vars]))))
    W <- build_w(source, target)
    if( any(sapply(source_dt[vars], Negate(is.numeric))) )
        stop("only supports numeric variables.")
    out <- est_w(W, source_dt[vars], target)
    return(out)
}

##' @name AI
##' @export
ai_var <- function(source, target,
                   vars, vars_var,
                   var_method = "CS") {
    stopifnot(inherits(source, "sf"))
    stopifnot(inherits(target, "sf"))
    stopifnot(inherits(vars, "character"))
    stopifnot(inherits(vars_var, "character"))
    stopifnot(inherits(var_method, "character"))
    stopifnot(length(vars) == 1 & length(vars_var) == 1)
    stopifnot(var_method %in% c("CS", "MI"))
    
    source_dt <- sf::st_drop_geometry(source)
    W <- build_w(source, target)
    if( any(sapply(source_dt[vars], Negate(is.numeric))) )
        stop("only supports numeric variables.")
    
    estimates <- W %*% as.matrix(source_dt[vars])
    target <- transform(target, est = as.numeric(estimates))
    if( var_method == "MI" ) {
        rho <- morans_i(source, vars_var)
        target <- var_w(W, source_dt[[vars_var]],
                        target, method = var_method,
                        rho_mi = rho)
    } else {
        target <- var_w(W, source_dt[[vars_var]],
                        target, method = var_method)
    }
    
}
