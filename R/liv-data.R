##' @title Liverpool Middle Super Output Area.
##'
##' @description A dataset containing containing the MSOA's for Liverpool along
##'     with estimates for Life Expectancy at Birth. Data taken from
##'     \href{https://ij-healthgeographics.biomedcentral.com/articles/10.1186/s12942-020-00200-w}{Johnson et al. 2020}
##' 
##' @details The data was projected to EPSG 27700 and units changed to km
##' 
##' @format A \code{sf} data frame with 61 rows and 4 variables:
##' \describe{
##'     \item{msoa11cd}{MSOA code} \item{msoa11cd}{MSOA name}
##'     \item{lev_est}{Estimated life expectancy at birth, in years}
##'     \item{area}{MSOA area, in \eqn{km^2}}
##'   }
##' 
##' @source \url{https://ij-healthgeographics.biomedcentral.com/articles/10.1186/s12942-020-00200-w}
"liv_msoa"

##' @title Liverpool Lower Super Output Area.
##'
##' @description A dataset containing containing the LSOA's for Liverpool along
##'     with estimates for Index of Multiple Deprivation. Data taken from
##'     \href{https://ij-healthgeographics.biomedcentral.com/articles/10.1186/s12942-020-00200-w}{Johnson et al. 2020}
##' 
##' @details The data was projected to EPSG 27700 and units changed to km
##' 
##' @format A \code{sf} data frame with 298 rows and 6 variables:
##' \describe{
##'     \item{lsoa11cd}{LSOA code}
##'     \item{lsoa11cd}{LSOA name} \item{male}{Male population}
##'     \item{female}{Female population}
##'     \item{imdscore}{Index of Multiple Deprivation}
##'     \item{area}{LMSOA area, in \eqn{km^2}}
##'   }
##' 
##' @source \url{https://ij-healthgeographics.biomedcentral.com/articles/10.1186/s12942-020-00200-w}
"liv_lsoa"
