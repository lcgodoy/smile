##' @title MLEs for fixed V.
##' @details internal use.
##' @param y response variable.
##' @param X design matrix.
##' @param Vinv inverse of \eqn{V}
est_mle <- function(y, X, Vinv) {
    XVinv <- crossprod(X, Vinv)
    beta  <- chol2inv(chol((XVinv %*% X))) %*% XVinv %*% y
    mu <- X %*% beta
    y  <- y - mu
    sigsq <- (crossprod(y, Vinv) %*% y) / n
    out <- c(as.numeric(beta), sigma2)
    names(out) <- c(sprintf("beta%d", seq(from = 0, NROW(beta) - 1)),
                    "sigsq" = sigsq)
    return(out)
}
