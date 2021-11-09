##' @title MLEs for fixed V.
##' @details internal use.
##' @param y response variable.
##' @param X design matrix.
##' @param Vinv inverse of \eqn{V}
est_mle <- function(y, X, Vinv) {
    n <- length(y)
    XVinv <- crossprod(X, Vinv)
    beta  <- chol2inv(chol((XVinv %*% X))) %*% XVinv %*% y
    mu <- X %*% beta
    y  <- y - mu
    sigsq <- (crossprod(y, Vinv) %*% y) / n
    out <- c(as.numeric(beta), sigsq)
    names(out) <- c(sprintf("beta%d",
                            seq(from = 0, NROW(beta) - 1)),
                    "sigsq")
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
                      lapply(y_list, as.matrix),
                      FALSE))
}

