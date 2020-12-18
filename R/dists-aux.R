##' @title Internal use only
##' @param mat_list internal use
##' @param mat_list1 internal use
##' @param mat_list2 internal use
##' @param return_single internal use
##' @param pred_mat internal use
##' @param x_to_list internal use
##' @param by internal use
##' @name aux_mat
single_dists <- function(mat_list) {
    n_out <- length(mat_list)

    out_dists <- vector(mode = "list", length = length(mat_list))
    
    k <- 0
    for(i in seq_len(n_out)) {
        for(j in i:n_out) {
            k <- k + 1
            ## include parameter for parallelization
            out_dists[[k]] <- crossdist(mat_list[[i]], mat_list[[j]])
        }
    }

    return(out_dists)
}

##' @name aux_mat
mult_dists <- function(mat_list1, mat_list2,
                       return_single = FALSE) {

    n_1 <- length(mat_list1)
    n_2 <- length(mat_list2)

    out_cross <- vector(mode = "list",
                        length = n_1 * n_2)
    k <- 0

    for(i in 1:n_1) {
        for(j in 1:n_2) {
            k <- k + 1
            out_cross[[k]] <- crossdist(mat_list1[[i]], mat_list2[[j]])
        }
    }

    if(return_single) {
        out_1 <- single_dists(mat_list1)
        out_2 <- single_dists(mat_list2)
        return(
            list(
                dists_1 = out_1,
                dists_2 = out_2,
                cross   = out_cross
            )
        )
    } else
        return(out_cross)
}

##' @name aux_mat
pred_cdist <- function(mat_list, pred_mat) {
    n_1  <- length(mat_list)
    n_pd <- NROW(pred_mat)

    out <- vector(mode = "list",
                  length = n_1 * n_pd)

    k <- 0
    for(i in 1:n_1) {
        for(j in 1:n_pd) {
            k <- k + 1
            out[[k]] <- crossdist(mat_list[[i]], pred_mat[j, , drop = FALSE])
        }
    }

    return(out)
}

##' @name aux_mat
get_grid_list <- function(x_to_list, by) {
    x_list <- sf::st_coordinates(x_to_list)
    
    x_list <- split(x = as.data.frame(x_list),
                    f = as.character(x_to_list[[by]]))
    
    return(lapply(x_list, as.matrix))
}
