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

##' @title New York City survey data.
##'
##' @description A dataset containing containing the census tracts for New York
##'     City along with estimates for median income and a margin of error for
##'     this estimates.
##' 
##' @details The data is project using EPSG 4326.
##' 
##' @format A \code{sf} data frame with 2128 rows and 5 variables:
##' \describe{
##'     \item{GEOID}{unique identifier}
##'     \item{NAME}{census tract name}
##'     \item{variable}{variable estimated}
##'     \item{estimate}{median income estimate}
##'     \item{moe}{median income estimate margin of error}
##'   }
"nyc_surv"

##' @title New York City community districts spatial geometries
##'
##' @description A dataset containing containing the CD's for New York City.
##' 
##' @details The data is project using EPSG 4326.
##' 
##' @format A \code{sf} data frame with 71 rows and 3 variables:
##' \describe{
##'     \item{boro_cd}{unique identifier}
##'     \item{shape_area}{Shape Area}
##'     \item{shape_length}{Shape Length}
##'     \item{est}{median income estimated using areal interpolation}
##'     \item{se_est}{standard error associated with the estimates}
##'   }
"nyc_comd"

##' @title Nova Lima census tracts
##'
##' @description A dataset containing containing the census tracts for the city
##'     of Nova Lima in Minas Gerais - Brazil.
##'
##' @details The data is project using the SIRGAS 2000.
##'
##' @format A \code{sf} data frame with 113 rows and 14 variables:
##' \describe{
##'     \item{cd_setor}{unique identifier}
##'     \item{hh_density}{average household density}
##'     \item{var_hhd}{variance of the household density}
##'     \item{avg_income}{average income per household}
##'     \item{var_income}{variance of the income per household}
##'     \item{pop}{population in the census tract}
##'     \item{avg_age}{average age of the inhabitants in the census tract}
##'     \item{var_age}{variance of the variable age in the census tract}
##'     \item{prop_women}{proportion of women}
##'     \item{prop_elder}{proportion of people with 55 years of age or older}
##'     \item{illit_rate}{illiteracy rate}
##'     \item{prop_white}{proportion of self-declared white people}
##'     \item{prop_black}{proportion of self-declared black people}
##'     \item{prop_native}{proportion of self-declared native people}
##'   }
"nl_ct"
