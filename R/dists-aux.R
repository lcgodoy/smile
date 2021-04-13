##' @name aux_mat
get_grid_list <- function(x_to_list, by) {
    x_list <- sf::st_coordinates(x_to_list)
    
    x_list <- split(x = as.data.frame(x_list),
                    f = as.character(x_to_list[[by]]))
    
    return(lapply(x_list, as.matrix))
}

##' @name aux_mat
dist_from_grids <- function(y_grid,  by) {
    out_list <- split(x = sf::st_set_geometry(y_grid, NULL),
                      f = as.character(y_grid[[by]]))
    return(
        single_dists(
            lapply(out_list,
                   function(x) as.matrix(x[ ,c(2, 3)]))
        )
    )
}

##' @name aux_mat
mult_dist_from_grids <- function(y_grid, x_grid, by) {
    y_list <- sf::st_coordinates(y_grid)
    
    y_list <- split(x = as.data.frame(y_list),
                    f = as.character(y_grid[[by[[1]]]]))

    x_list <- sf::st_coordinates(x_grid)
    
    x_list <- split(x = as.data.frame(x_list),
                    f = as.character(x_grid[[by[[2]]]]))
    
    return(mult_dists(lapply(y_list, as.matrix),
                      lapply(y_list, as.matrix),
                      FALSE))
}

